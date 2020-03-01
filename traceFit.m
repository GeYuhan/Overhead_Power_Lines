function [F] = traceFit(frequency,tau,Data,option)
%[F] = traceFit(frequency,tau,Data,option)
%   输入：
%   frequency
%   tau由延时提取而来,且已分好组
%   Data是输入的迹，那么是一个行向量（不是的话，给它转化成行向量）
%   option是控制结构体，含有以下字段
%       Np          指定极点数目,必须小于采样点数（frequency的长度）
%       passivity   是否进行无源性检查，并强迫无源（1=YES，0=No）
%       weight      最小二乘拟合权重(0=不设任何权重，1=提高基频处响应精度)
%       realFit     是否采用实极点拟合（1=Yes，0=No）
%       relax       是否采用宽松方法（1=yes，0=no）
%   输出：
%   F是一结构体
%       F.T最终的延时
%       F.P最终的极点
%       F.C留数
%       F.D常数项
%       F.MFE最大拟合误差
%       F.RMSE均方根误差

%% 数据处理
[Nc Nc Ns]=size(Data);
Ns=length(frequency);
%%保证frequency是列向量
[row column]=size(frequency);
if row <column
    frequency=frequency.';
end
Response=reshape(Data,length(Data),1);%重新构造的响应矩阵，保证为列向量
s=2i*pi*frequency;%列向量，Ns by 1
%% 确定极点
Np=option.Np;
% disp(['极点数目为 ' num2str(Np)])
poles=initialPoles(frequency,Np,1);%2表示生成实极点，1复极点
%% 形成计算用的延时向量
Ng=length(tau);
t=[];
for i=1:Ng
    t=[t,repmat(tau(i),1,Np)];
end
%% 普通最小二乘
%% 确定权重
switch option.weight
    case 0
        Weight=ones(Ns,1);%有其他权重再加
    case 1
        for i=1:length(frequency)%pscad权重，给50Hz加1000，保留上一句是为了生成列向量
            if frequency(i)>50&frequency(i-1)<=50
                Weight(i-1)=1000;
                Weight(i)=1;
            else
                Weight(i)=1;
            end
        end
        Weight=reshape(Weight,Ns,1);%上面生成的是行向量，这里需要列向量来填充
end
% disp(['最小二乘权重为1'])
%% 大系数矩阵定义 AAx=bb
if option.relax
   AA=zeros(Np+1);
   bb=zeros(Np+1,1);
else
    AA=zeros(Np);
    bb=zeros(Np,1);
end
b=Weight.*Response;%给响应加权
b=[real(b);imag(b)];
scale=(norm(b))^2;%每个导体一个权重向量，relax参数 
scale=sqrt(scale)./Ns; %对所有元素的平方和进行开方，求平均
%% 迭代求解
Nimax=10;
tolerance=0.5;%用于判断极点是否收敛
x=0;
for i=1:Nimax
    disp(['迭代第 ' num2str(i) '次'])
        %% 标记复数极点
    cindex=zeros(size(poles));%实极点记为零，复极点记为1，其共轭（紧接着下一位）记为2
    %例子：         poles=（[1 0 i 1 0 i]）,
    %结果cindex长度为7，    [0 0 1 2 0 1 2]
    for m=1:length(poles)
        if imag(poles(m))~=0
            if m==1
                cindex(m)=1;
            else
                if cindex(m-1)==0 || cindex(m-1)==2
                    cindex(m)=1; cindex(m+1)=2;
                else
                    cindex(m)=2;
                end
            end
        end
    end
    %%构造系数矩阵
    Ncol=2*Np;%小系数矩阵A的列数
    Ncol1=Np;%小系数矩阵A前面与σ无关的部分
    
    H=ones(Ns,Ncol);%中间矩阵H，用来存放h（s）
    P=ones(Ns,Ncol);%中间矩阵P，用来存放1/（s-p）
    A=zeros(Ns,Ncol);%小系数矩阵A
    
    P=repmat(s,1,2*Np);
    A=repmat(poles,Ns,2);%借用一下A矩阵
    P=1./(P-A);
    if option.relax
        P(:,end+1)=ones(Ns,1);
    end
    %%对复极点做修改
    for m=1:length(poles)
        if cindex(m)==1
            P1=P(:,m);
            P2=P(:,m+1);
            P(:,m)=P1+P2;%
            P(:,m+1)=(P1-P2)*1i;
            P(:,m+length(poles))=P1+P2;
            P(:,m+1+length(poles))=(P1-P2)*1i;
        end
    end%相应的求新极点也要修改
    H(:,Ncol1+1:end)=repmat(-1.*Response,1,Np);
    if option.relax
        H(:,end+1)=-1.*Response;
    end
    A=H.*P;%%得到了小系数矩阵
    %将A扩充为迹拟合的形式
    delayFactorMatrix=s*t;%这个得是一个矩阵
    Atemp=repmat(A(:,1:Ncol1),1,Ng).*exp(-delayFactorMatrix);
    A=[Atemp,A(:,Ncol1+1:end)];%此时A为迹拟合的正确系数矩阵
    %%对系数矩阵A加权
    [~,colA]=size(A);
    W=repmat(Weight,1,colA);
    A=W.*A;%加权
    %虚实分解以保证共轭（b已在循环外分解过）
    A=[real(A);imag(A)];
    Ncol2=Np*Ng;
    Ncol3=Np*(Ng+1);
    if option.relax
        A(end+1,:)=[zeros(1,Ncol2), real(sum(P (:,Ncol1+1:end)))*scale];%Ncol1=Np
        [Q R]=qr(A,0);%QR分解，提取出与σ相关的部分
        R22=R(Ncol2+1:end,Ncol2+1:end);
        AA(1:Np+1,:)=R22;
        %             size(bb((n-1)*(Np+1)+1:n*(Np+1),1))
        %             size(Q(end,Ncol1+1:end)'*Ns*scale)
        bb(1:Np+1,1)=Q(end,Ncol2+1:end)'*Ns*scale;
    else
        [Q R]=qr(A,0);%QR分解，提取出与σ相关的部分
        R22=R(Ncol2+1:Ncol2+Np,Ncol2+1:Ncol2+Np);
        AA(1:Np,:)=R22;
        bb(1:Np,1)=Q(:,Ncol2+1:Ncol2+Np)'*b;
    end
    %如果是relax方法，这里要进行修改
    %% 改善AA矩阵条件数
    [rowOfAA colOfAA]=size(AA);
    col=1:colOfAA;%AA矩阵的列数
    rnormArray=arrayfun(@(col) norm(AA(:,col),2),col);%%行向量，暂用rnormArray存放AA每一列向量的2范数
    rnormArray=1./(rnormArray);%得到所需的倒数
    CondScale=repmat(rnormArray,rowOfAA,1);%构造出缩放矩阵
    AA=AA.*CondScale;%%缩放完毕
    %% 求解最小二乘方程
    x=AA\bb;%求出了σ的留数,[c1~,c2~,...]
    %% 还原AA改善条件数之前的解
    x=x.*rnormArray.';
    if sum(abs(x)==inf)>=1
        x((abs(x)==inf))=sign(x((abs(x)==inf)))*realmax;%如果某一数超过了机器最大值，就将其职位最大值。否则inf会在后面的求新极点eig中报错
    end
    if option.relax
        if x(end)==0
            x(end)=1;
            %         elseif x(end)<Tol
        end
    else
        x(end+1)=1;%应该为Dnew
    end
    %% 计算新极点
    X=repmat(x(1:end-1),1,Np);%这里只需要留数部分的大小
    lambda=diag(poles);
    %对lambda中的复极点进行修改
    for m=1:length(poles)
        if cindex(m)==1
            lambda(m:m+1,m:m+1)=[real(poles(m)),imag(poles(m));-imag(poles(m)),real(poles(m))];
            X(m:m+1,m:m+1)=[2*x(m),2*x(m+1);0,0];
        end
    end
    D=lambda-X./x(end);%公式中的d~^-1
    newPoles=eig(D).';
    %%翻转不稳定极点
    unstables=real(newPoles)>0;%存储不稳定极点的坐标
    if sum(unstables)>0%说明有不稳定极点
        disp(['存在不稳定极点，将翻转'])
        newPoles(unstables)=newPoles(unstables)-2*real(newPoles(unstables));%翻转不稳定极点
    end
    %%判断收敛
    if norm(newPoles-poles)/norm(poles)<=tolerance
        %         disp(['极点已收敛 tolerance=' num2str(tolerance)])
        break
    else
        poles=newPoles;
    end
end%迭代求解
if i==Nimax
    disp(['矢量拟合达到最大迭代次数',num2str(Nimax)])
    %Np=Np+2;
    %这里应该返回重新计算，暂先不写
end
%% 输出计算
Ne=1;%为了简便，后面的Ne我不删了，在迹拟合中它就是1
option.pattern=0;%同样为了简便，只是对传播函数的迹拟合，无常数项
%%储存频率
F.freq=frequency;%列向量
F.T=t;
%%矩阵定义
F.C=zeros(Ne,1,Np*Ng);%存放留数，每一页对应一个极点
F.D=zeros(Ne,1);%存放常数项
%%极点
F.P=repmat(poles,1,Ng);%行向量
%%留数计算,每个元素的系数矩阵都一样
P=ones(Ns,Ncol1);%中间矩阵P，用来存放1/（s-p）
A=zeros(Ns,Ncol1);
P(:,1:end)=repmat(s,1,Np);
A(:,1:end)=repmat(poles,Ns,1);%借用一下A矩阵
P=1./(P-A);
for m=1:length(poles)
    if cindex(m)==1
        P1=P(:,m);
        P2=P(:,m+1);
        P(:,m)=P1+P2;%
        P(:,m+1)=(P1-P2)*1i;
    end
end%使得留数共轭，先变成实数
A=P;%得到系数矩阵，Ns by Np
%将A扩充为迹拟合的形式
A=repmat(A(:,1:end),1,Ng).*exp(-delayFactorMatrix);%此时A为迹拟合的正确系数矩阵
[~,colA]=size(A);
W=repmat(Weight,1,colA);
A=W.*A;%加权
A=[real(A);imag(A)];
%b=Weight.*Response;前面已响应加权
%% 改善A矩阵条件数
[rowOfA colOfA]=size(A);
col=1:colOfA;%A矩阵的列数
rnormArray=arrayfun(@(col) norm(A(:,col),2),col);%%行向量，暂用rnormArray存放A每一列向量的2范数
rnormArray=1./(rnormArray);%得到所需的倒数
CondScale=repmat(rnormArray,rowOfA,1);%构造出缩放矩阵
A=A.*CondScale;%%缩放完毕
%% MOR降阶,筛选出子极点集
[Srow,Scol]=size(A);
Sindex=min(Srow,Scol);
[U,S,V]=svd(A);
save traceFit S A
%         save svdResidue U S V
%         save datA 
while S(Sindex,Sindex)<1e-14
    S(Sindex,Sindex)=0;
    Sindex=Sindex-1;
end
A=U*S*V';
%%
for n=1:Ne %分别计算每个元素的留数
    x=A\b(:,n);%求得留数,列向量,[d,c1,c2,...]
    x=x.*rnormArray.';%恢复条件数优化之前的解
    for m=1:length(poles)
        if cindex(m)==1
            x1=x(m+option.pattern);
            x2=x(m+option.pattern+1);
            x(m+option.pattern)=x1+x2*1i;
            x(m+option.pattern+1)=x1-x2*1i;
        end
    end%将共轭留数变为复数
    if option.pattern==1
        F.D(n,1)=x(1);%对Yc来说，第一项是d
        F.C(n,1,:)=x(2:end);
    else
        F.D(n,1)=0;%对Hmps来说，无d
        F.C(n,1,:)=x(1:end);
    end
%     F.D(n,1)=abs(F.D(n,1))*real(F.D(n,1))/abs(real(F.D(n,1)));%看来不虚实分解，只是简单的取实部，会丢失相位的信息
end
%% 提取子极点
index=(F.C~=0);
F.C=F.C(index);
F.T=F.T(index);
F.P=F.P(index);
% size(F.C)
%% 拟合曲线数据、RMSE与MFE（maxFittingError）计算，代码移植于dataPlot
%生成拟合曲线数据
f=zeros(size(Data));
for i=1:Ns
    for j=1:length(F.P)
        f(i)=f(i)+F.C(:,:,j)./(s(i)-F.P(j));
    end
    f(i)=f(i)+F.D; 
end
%画幅值的函数
% figure;
% for i=1:Nc
%     for j=1:Nc
%         semilogx(F.freq,abs(squeeze(Data(i,j,:))),'b'); 
%         hold on
%         semilogx(F.freq,abs(squeeze(f(i,j,:))),'r.');
%     end
% end
% hold off
% xlabel('频率 Hz')
% ylabel('幅值 ')
% legend('原数据','拟合结果')
%% 画误差的函数
% figure;
%这里要计算RMSE
residualError=abs(Data-f);%每个点的残差
F.MFE=max(residualError)/max(abs(Data));%找出最大的拟合误差并储存,依据pscad定义
SSR=sum(residualError.^2);%残差平方和
F.RMSE=sqrt(SSR/length(residualError));%均方根误差，并储存
%         semilogx(F.freq,residualError,'r'); %%残差
%这里要计算最大拟合误差MFE
%         hold on


% title('幅值拟合残差')
% xlabel('频率 Hz')
% ylabel('幅值 ')
% hold off
 %% 画相角的函数
% figure;
% for i=1:Nc
%     for j=1:Nc
%         semilogx(F.freq,180*unwrap(angle(squeeze(Data(i,j,:))))/pi,'b'); 
%         hold on
%         semilogx(F.freq,180*unwrap(angle(squeeze(f(i,j,:))))/pi,'r.');
%     end
% end
% hold off
% xlabel('频率 Hz')
% ylabel('相角 度')
% legend('原数据','拟合结果')
end



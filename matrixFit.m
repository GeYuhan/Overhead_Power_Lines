function [F] = matrixFit(frequency,Data,option)
%   F = matrixFit(frequency,Data)
%要求输入frequency是列向量，Ns by 1，单位Hz。如果不是列向量，会自动转置为列向量
%Data是三维矩阵，Nc by Nc by Ns，可以是复矩阵
%option是控制结构体，含有以下字段
%	pattern     选择拟合形式（Yc/HmMPS）
%	Np          指定极点数目,必须小于采样点数（frequency的长度）
%	passivity   是否进行无源性检查，并强迫无源（1=YES，0=No）
%	weight      最小二乘拟合权重(0=不设任何权重，1=提高基频处响应精度)
%	realFit     是否采用实极点拟合（1=Yes，0=No）
%   relax       是否采用宽松方法（1=yes，0=no）
%
%下面是详细解释：
%option.pattern用来选择是Yc数据，还是Hmps数据
%当option.pattern==1,是Yc数据
%       F按如下形式拟合：
%              N     Cn
%       F(s)= SUM  -------  + D 
%             n=1   s-pn
%      
%当option.pattern==0，是Hmps数据
%       F按如下形式拟合：
%              N     Cn
%       F(s)= SUM  -------
%             n=1   s-pn
%      
%F是结构体
%   F.freq==frequency，频率列向量 1 Ns
%   F.P是极点行向量，1 by Np
%   F.C是留数矩阵，三维Nc Nc Np，每一页对应一个极点
%   F.D是常数项矩阵，Nc Nc，H为零
%   F.T是延时向量，Nk*Np，Yc为零，H为非零。由于本程序是对Yc与最小相位HmMPS拟合，故均为零。
%   F.MFE是最大拟合误差（maxFittingError），标量。幅值与相角不区分？现在只有幅值
%   F.RMSE是最大的均方根误差，标量。现在只算幅值
%
%下面是拟合中的设置：
%极点数目当前人为指定 Np
%最大迭代次数人为设定 Nimax，收敛判据 tolerance=0.01
%没有写拟合精度判断，用于增加极点数目。【函数外实现】
%使用relax方法
%使用条件数优化，使用虚实分解
%最小二乘权重还未指定，当前为1 weight【已改为pscad权重】
%
%% 数据处理
[Nc Nc Ns]=size(Data);
Ns=length(frequency);
%%保证frequency是列向量
[row column]=size(frequency);
if row <column
    frequency=frequency.';
end
%%三维矩阵变二维
symmetry=0;%默认为不对称
if isequal(Data(:,:,fix(Ns/2)),Data(:,:,fix(Ns/2)).')
    symmetry=1;
    Ne=(Nc+1)*Nc/2;
    Data1=zeros(1,Ne,Ns);
    map=logical(triu(ones(Nc)));
    for n=1:Ns
        a=Data(:,:,n);
        Data1(:,:,n)=reshape(a(map),1,Ne);
    end  %Data1页频率，列元素，一行
else
    Ne=Nc*Nc;
    for n=1:Ns
        Data1(:,:,n)=reshape(Data(:,:,n),1,Ne);
    end  %Data1页频率，列元素，一行
end
Data2=shiftdim(Data1,1);%Data2 Ne行元素，Ns列频率

Response=Data2.';%重新构造的响应矩阵，Ns by Ne
s=2i*pi*frequency;%列向量，Ns by 1
%% 确定极点
Np=option.Np;
% disp(['极点数目为 ' num2str(Np)])
poles=initialPoles(frequency,Np,1);%2表示生成实极点，1复极点
%% 普通最小二乘
%% 确定权重
switch option.weight
    case 0
        Weight=ones(Ns,Ne);%有其他权重再加
    case 1
        for i=1:length(frequency)%pscad权重，给50Hz加1000，保留上一句是为了生成列向量
            if frequency(i)>50&frequency(i-1)<=50
                Weight(i-1)=1000;
                Weight(i)=1;
            else
                Weight(i)=1;
            end
        end
        Weight=repmat(Weight.',1,Ne);%上面生成的是行向量，这里需要列向量来填充
end
% disp(['最小二乘权重为1'])
%% 大系数矩阵定义 AAx=bb
if option.relax
   AA=zeros(Ne*(Np+1),Np+1);
   bb=zeros(Ne*(Np+1),1);
else
    AA=zeros(Ne*Np,Np);
    bb=zeros(Ne*Np,1);
end
b=Weight.*Response;%给响应加权
b=[real(b);imag(b)];
scale=0;%relax参数
for m=1:Ne
    scale=scale+(norm(b(:,m)))^2;%每个导体一个权重向量 
end
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
    %% 大系数矩阵AA与bb的构造
    for n=1:Ne
        Ncol=2*Np+option.pattern;%小系数矩阵A的列数
        Ncol1=Np+option.pattern;%小系数矩阵A前面与σ无关的部分
    
        H=ones(Ns,Ncol);%中间矩阵H，用来存放h（s）
        P=ones(Ns,Ncol);%中间矩阵P，用来存放1/（s-p）
        A=zeros(Ns,Ncol);%小系数矩阵A

        P(:,1+option.pattern:end)=repmat(s,1,2*Np);
        A(:,1+option.pattern:end)=repmat(poles,Ns,2);%借用一下A矩阵
        P=1./(P-A);
        if option.relax
            P(:,end+1)=ones(Ns,1);
        end
        %%对复极点做修改
        for m=1:length(poles)
            if cindex(m)==1
                P1=P(:,m+option.pattern);
                P2=P(:,m+option.pattern+1);
                P(:,m+option.pattern)=P1+P2;%
                P(:,m+option.pattern+1)=(P1-P2)*1i;
                P(:,m+option.pattern+length(poles))=P1+P2;
                P(:,m+option.pattern+1+length(poles))=(P1-P2)*1i;
            end
        end%相应的求新极点也要修改
        H(:,Ncol1+1:end)=repmat(-1.*Response(:,n),1,Np);
        if option.relax
            H(:,end+1)=-1.*Response(:,n);
        end
        A=H.*P;%%得到了小系数矩阵
        [~,colA]=size(A);
        W=repmat(Weight(:,n),1,colA);
        A=W.*A;%加权
        %虚实分解以保证共轭（b已在循环外分解过）
        A=[real(A);imag(A)];
        if option.relax
            A(end+1,:)=[zeros(1,Ncol1), real(sum(P (:,Ncol1+1:end)))*scale];
            [Q R]=qr(A,0);%QR分解，提取出与σ相关的部分
            R22=R(Ncol1+1:end,Ncol1+1:end);
            AA((n-1)*(Np+1)+1:n*(Np+1),:)=R22;
%             size(bb((n-1)*(Np+1)+1:n*(Np+1),1))
%             size(Q(end,Ncol1+1:end)'*Ns*scale)
            bb((n-1)*(Np+1)+1:n*(Np+1),1)=Q(end,Ncol1+1:end)'*Ns*scale;
        else
            [Q R]=qr(A,0);%QR分解，提取出与σ相关的部分
            R22=R(Ncol1+1:Ncol1+Np,Ncol1+1:Ncol1+Np);
            AA((n-1)*Np+1:n*Np,:)=R22;
            bb((n-1)*Np+1:n*Np,1)=Q(:,Ncol1+1:Ncol1+Np)'*b(:,n);
        end
        %如果是relax方法，这里要进行修改
    end
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
%%储存频率
F.freq=frequency;%列向量
F.T=0;%将延时项置为零，后面不予理会。
%%矩阵定义
F.C=zeros(Ne,1,Np);%存放留数，每一页对应一个极点
F.D=zeros(Ne,1);%存放常数项
%%极点
F.P=poles;%行向量
%%留数计算,每个元素的系数矩阵都一样
P=ones(Ns,Ncol1);%中间矩阵P，用来存放1/（s-p）
A=zeros(Ns,Ncol1);
P(:,1+option.pattern:end)=repmat(s,1,Np);
A(:,1+option.pattern:end)=repmat(poles,Ns,1);%借用一下A矩阵
P=1./(P-A);
for m=1:length(poles)
    if cindex(m)==1
        P1=P(:,m+option.pattern);
        P2=P(:,m+option.pattern+1);
        P(:,m+option.pattern)=P1+P2;%
        P(:,m+option.pattern+1)=(P1-P2)*1i;
    end
end%使得留数共轭，先变成实数
A=P;%得到系数矩阵，Ns by Np+1
W=repmat(Weight(:,n),1,Ncol1);%这里n=Ne，是上面的某个循环结束留下的结果。本来应该是不同的元素不同的权重，可在实际操作中每个元素权重都一样，也就无所谓差别了
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
%% 将向量恢复成原有矩阵形式
if symmetry%如果是对称
    index=1;%用来存放每一段向量的起始位置
    C=zeros(Nc,Nc,Np);%用来存放上三角阵
    D=zeros(Nc,Nc);
    %%生成上三角阵
    for n=1:Nc
        len=Nc-n+1;%length，与matlab已有函数名冲突，故改为len
        C(n,n:end,:)=F.C(index:index+len-1,1,:);
        D(n,n:end)=F.D(index:index+len-1,1);
        index=index+len;      
    end
    %%由上三角阵得原矩阵
    i=1:Np;
    cellC=arrayfun(@(a) triu(C(:,:,a),1).'+C(:,:,a),i,'UniformOutput',false);%cellC是一个元胞数组，每一个元胞存储一个极点对应的留数矩阵
    F.C=cat(3,cellC{:});%将元胞数组转换为三维矩阵
    F.D=triu(D,1).'+D;%将D矩阵的行列转置并相加得到最后的常数项矩阵；因为D是二维矩阵，故没有C那么麻烦   
else%如果是非对称
    F.D=reshape(F.D,Nc,Nc);
    F.C=reshape(F.C,Nc,Nc,Np);
end

%% 拟合曲线数据、RMSE与MFE（maxFittingError）计算，代码移植于dataPlot
%生成拟合曲线数据
f=zeros(size(Data));
for i=1:Ns
    for j=1:Np
        f(:,:,i)=f(:,:,i)+F.C(:,:,j)./(s(i)-F.P(j));
    end
    f(:,:,i)=f(:,:,i)+F.D; 
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
for i=1:Nc
    for j=1:Nc
        %这里要计算RMSE
        residualError=abs(squeeze(Data(i,j,:)-f(i,j,:)));%每个点的残差
        MFE(i,j)=max(residualError)/max(abs(Data(i,j,:)));%找出最大的拟合误差并储存,依据pscad定义
        SSR=sum(residualError.^2);%残差平方和
        RMSE(i,j)=sqrt(SSR/length(residualError));%均方根误差，并储存
%         semilogx(F.freq,residualError,'r'); %%残差
        %这里要计算最大拟合误差MFE
%         hold on
    end
end
%找出最大的最大拟合误差并输出到拟合结构体
F.MFE=max(max(MFE));
F.RMSE=max(max(RMSE));%所有元素最大的均方误差输出到结构体
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
    

    











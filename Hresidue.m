function F = Hresidue(freq,tau,Poles,Data,doMOR)
%X = Hresidue(freq,tau,Poles,data)
%   计算留数的函数
%-------------------输出-------------------
%   F是结构体：
%         F.freq 存储频率
%         F.P 存储极点向量
%         F.C 存储留数矩阵
%         F.T
%         存储延时τ向量，已对相量做了处理，与F.C三维长度相同，使得F.C（:,:,i）对应于F.T(i).留数程序中F.T不为标量零（matrixFit中为标量零）
%         F.MFE是最大拟合误差（maxFittingError），标量。幅值与相角不区分？现在只有幅值
%         F.RMSE是最大的均方根误差，标量。现在只算幅值
%         F.MRPR是最大留数极点比，标量。复数之比取绝对值
%   ----------------输入---------------------------      
%   freq是频率，Ns长度，需要列向量，如果不是，会自动转为列向量来生成s
%   tau是延时向量，行列无所谓，Nk长度
%   Poles是=元胞=极点矩阵，要求每一行对应一个τ,时延组数为 numel（Poles），即行数目
%   Data是数据,三维矩阵 Nc*Nc*Ns
%   doMOR:是否进行降阶？0——不降阶；1——降阶

%% 输入处理
%%对极点处理
% [Nk Np]=size(Poles);
Nk=numel(Poles);%读取延时组数目
% poles=reshape(Poles.',1,Nk*Np);%将所有极点组合成一个行向量，按延时顺序排列[p11,p12,p21,p22,p31,p32]
poles=[];
for i=1:Nk
    Np(i)=length(Poles{i});%延时组对应的极点数目
    poles=[poles,reshape(Poles{i},1,Np(i))];%将所有极点组合成一个行向量，按延时顺序排列[p11,p12,p21,p22,p31,p32]
end
Npp=length(poles);%极点总数目
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
%%对频率处理
[Nc,Nc,Ns]=size(Data);
Ns=length(freq);
freq=reshape(freq,Ns,1);%将频率转为列向量
s=2i*pi*freq;%得到列向量的s
%%对Data处理，三维矩阵变二维
%%
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
%% 加权
% Weight=sqrt((Response));
Weight=ones(Ns,Ne);
% Weight=repmat(Weight,1,Ne);
% Response=Response.*Weight;%响应加权
%% 留数计算
%%矩阵定义
F.freq=freq;%列向量
F.P=poles;%极点行向量
F.C=zeros(Ne,1,Npp);%存放留数，每一页对应一个极点
F.D=zeros(Nc,Nc);%由于本程序专用来计算H的留数，故此项为零矩阵
F.T=[];%我的设想就是每个留数矩阵对应一个τ，这样他就是Npp=Np(1)+Np(2)+Np(3)长度

t=[];%t矩阵存放着τ
for i=1:Nk%延时τ的数目,%得知道每个模量的极点数目-----------------------------------
    t=[t,repmat(tau(i),1,Np(i))];%T矩阵存放着exp（-sτ）
end
F.T=t;
%2019/3/4
if 0
%% 大系数矩阵定义 AAx=bb
rowNum=min(2*Ns+1-Npp,Npp+1);%考虑到欠定方程组情况，取二者最小
AA=zeros(Ne*rowNum,Npp+1);
bb=zeros(Ne*rowNum,1);
scale=0;%relax参数
for m=1:Ne
    scale=scale+(norm(Response(:,m)))^2;%每个导体一个权重向量 
end
scale=sqrt(scale)./Ns; %对所有元素的平方和进行开方，求平均
Nimax=4;%最大迭代次数
tolerance=0.5;%用于判断极点是否收敛
x=0;
for i=1:Nimax
    %%构造系数矩阵,每个元素都相同
    H=ones(Ns,2*Npp+1);%中间矩阵H，用来存放h（s）
    P=repmat(s,1,2*Npp);%中间矩阵P，用来存放1/（s-p）
    A=repmat(poles,Ns,2);%借用一下A矩阵
    P=1./(P-A);%存放了 1/（s-p）
    P(:,end+1)=ones(Ns,1);%relax方法所需
    for m=1:length(poles)
        if cindex(m)==1
            P1=P(:,m);
            P2=P(:,m+1);
            P(:,m)=P1+P2;%
            P(:,m+1)=(P1-P2)*1i;
            P(:,m+length(poles))=P1+P2;
            P(:,m+1+length(poles))=(P1-P2)*1i;
        end
    end%使得留数共轭，先变成实数
    for n=1:Ne%构造大矩阵AA
        H(:,Npp+1:end)=repmat(-1.*Response(:,n),1,Npp+1);
        A=H.*P;
        A(:,1:Npp)=A(:,1:Npp).*exp(-s*t);%得到系数矩阵，Ns by Npp, 存放了 exp(-sτ)/（s-p）
        [~,Acol]=size(A);
        % W=repmat(Weight(:,1),1,Acol);
        % A=A.*W;%如果是每个元素权重不同的话，没计算一个元素，就得改权重
        A=[real(A);imag(A)];
        A(end+1,:)=[zeros(1,Npp), real(sum(P (:,Npp+1:end)))*scale];
        [Q R]=qr(A,0);%QR分解，提取出与σ相关的部分
        R22=R(Npp+1:end,Npp+1:end);
        AA((n-1)*rowNum+1:n*rowNum,:)=R22;
        bb((n-1)*rowNum+1:n*rowNum,1)=Q(end,Npp+1:end)'*Ns*scale;
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
    if x(end)==0
        x(end)=1;
    end
    %% 计算新极点
    X=repmat(x(1:end-1),1,Npp);%这里只需要留数部分的大小
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
end
end
%% 构造小矩阵A
%%构造系数矩阵,每个元素都相同
H=ones(Ns,Npp);%中间矩阵H，用来存放h（s）
P=repmat(s,1,Npp);%中间矩阵P，用来存放1/（s-p）
A=repmat(poles,Ns,1);%借用一下A矩阵
P=1./(P-A);%存放了 1/（s-p）
% A=P.*exp(-s*t);
for m=1:length(poles)
    if cindex(m)==1
        P1=P(:,m);
        P2=P(:,m+1);
        P(:,m)=P1+P2;%
        P(:,m+1)=(P1-P2)*1i;
    end
end%使得留数共轭，先变成实数
%%2019/3/7
A=P.*exp(-s*t);%得到系数矩阵，Ns by Npp, 存放了 exp(-sτ)/（s-p）
save math Response
% W=repmat(Weight(:,i),1,Acol);
% A=A.*W;%如果是每个元素权重不同的话，没计算一个元素，就得改权重
A=[real(A);imag(A)];
%% 改善A矩阵条件数
[rowOfA colOfA]=size(A);
col=1:colOfA;%A矩阵的列数
rnormArray=arrayfun(@(col) norm(A(:,col),2),col);%%行向量，暂用rnormArray存放A每一列向量的2范数
rnormArray=1./(rnormArray);%得到所需的倒数
CondScale=repmat(rnormArray,rowOfA,1);%构造出缩放矩阵
A=A.*CondScale;%%缩放完毕
%% 分别计算每个元素的留数
for n=1:Ne
    b=Response(:,n).*Weight(:,n);%第n个元素的响应数据
    b=[real(b);imag(b)];
    %%求解
    x=A\b;%求得留数,列向量,[c11,c12,...c21,c22,...]
    x=x.*rnormArray.';%%还原A改善条件数之前的解  
    for m=1:length(poles)
        if cindex(m)==1
            x1=x(m);
            x2=x(m+1);
            x(m)=x1+x2*1i;
            x(m+1)=x1-x2*1i;
        end
    end%将共轭留数变为复数
    F.C(n,1,:)=x(1:end);
    MRPRatio(n)=max(abs(x)./abs(poles).');%最大留数极点比(Max Residue/Pole Ratio),x是列向量，poles是行向量，最终MRPRatio也是列向量
end
F.MRPR=max(MRPRatio);%变为标量输出
%% 判断是否进行降阶。MOR降阶,可防止过大的留数极点比
if doMOR
    disp(['留数极点比超过阈值，进行降阶']);
    [Srow,Scol]=size(A);
    Sindex=min(Srow,Scol);
    global MRPRtolerance %MRPRtolerance在fittingSet定义
    while F.MRPR>MRPRtolerance
        [U,S,V]=svd(A);
        save svdResidue U S V
        save datA         
        S(Sindex,Sindex)=0;
        Sindex=Sindex-1;
        A=U*S*V';
        %% 分别计算每个元素的留数
        for n=1:Ne
            b=Response(:,n);%第n个元素的响应数据
            b=[real(b);imag(b)];
            %%求解
            x=A\b;%求得留数,列向量,[c11,c12,...c21,c22,...]
            x=x.*rnormArray.';%%还原A改善条件数之前的解
            %     AAA*x-b<=eps
            
            for m=1:length(poles)
                if cindex(m)==1
                    x1=x(m);
                    x2=x(m+1);
                    x(m)=x1+x2*1i;
                    x(m+1)=x1-x2*1i;
                end
            end%将共轭留数变为复数
            F.C(n,1,:)=x(1:end);
            MRPRatio(n)=max(abs(x)./abs(poles).');%最大留数极点比(Max Residue/Pole Ratio),x是列向量，poles是行向量，最终MRPRatio也是列向量
        end
        F.MRPR=max(MRPRatio);%变为标量输出
    end
    %% 提取剩余极点
    index=(F.C(1,1,:)~=0);%只看左上角的元素即可
    F.C=F.C(:,:,index);
    F.T=F.T(index);
    F.P=F.P(index);
end%MOR
Npp=length(F.P);
%% 将向量恢复成原有矩阵形式
if symmetry%如果是对称
    index=1;%用来存放每一段向量的起始位置
    C=zeros(Nc,Nc,Npp);%用来存放上三角阵----------------------?为什么是Np，看下面应该是总数目
    %%生成上三角阵
    for n=1:Nc
        len=n;%length，与matlab已有函数名冲突，故改为len
        C(1:n,n,:)=F.C(index:index+len-1,1,:);
        index=index+len;      
    end
    %%由上三角阵得原矩阵
    i=1:Npp;%已修改为总极点数目
    cellC=arrayfun(@(a) triu(C(:,:,a),1).'+C(:,:,a),i,'UniformOutput',false);%cellC是一个元胞数组，每一个元胞存储一个极点对应的留数矩阵
    F.C=cat(3,cellC{:});%将元胞数组转换为三维矩阵  
else%如果是非对称
    F.C=reshape(F.C,Nc,Nc,Npp);%已修改为总数目
end
%% 拟合曲线数据、RMSE与MFE（maxFittingError）计算，代码移植于dataPlot
%生成拟合曲线数据
% Np=length(F.P);%前面拟合用的Np是每个模量的，算数据的Np是一共的极点数目
f=zeros(size(Data));
for i=1:Ns
    for j=1:Npp
        f(:,:,i)=f(:,:,i)+exp(-s(i)*F.T(j)).*F.C(:,:,j)./(s(i)-F.P(j));
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
%% 求误差的函数
% figure;
% denominator=max(max(abs(Data(i,j,:)))
for i=1:Nc
    for j=1:Nc
        %这里要计算RMSE
        residualError=abs(squeeze(Data(i,j,:)-f(i,j,:)));%每个点的残差
%         residualError=abs(squeeze(abs(Data(i,j,:))-abs(f(i,j,:))));%每个点的残差，是模值的残差
        MFE(i,j)=max(residualError)/max(max(max(abs(squeeze(Data)),[],3)));%找出最大的拟合误差并储存,依据pscad定义
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
% F.MFE=MFE;
% F.RMSE=RMSE;%所有元素最大的均方误差输出到结构体
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


%%矢量拟合脚本
%% 特征导纳Yc与传播矩阵(H与Hm)计算
[H, Hm, Ti, gamma, minusGammaL,Yc]=deal(zeros(size(Y)));
traceH=zeros(1,length(freq));
omega=2*pi*freq;%计算ω
s=omega*1i;%s==jω
for i=1:length(s)
    sqrtYZ=sqrtm(Y(:,:,i)*Z(:,:,i));%计算√YZ
    Yc(:,:,i)=sqrtYZ/Z(:,:,i);%Yc=√YZ Z^-1
    exponent=-sqrtYZ*l;%exp（A）中的A
    H(:,:,i)=expm(exponent);%得到传播矩阵
    traceH(i)=trace(H(:,:,i));
    for j=1:length(traceH)
        if traceH(j)==0
            traceH(j)=exp(-745);
        end
    end  
    %%用来防止拟合算法中出现nan，将零替换为最小数
    [row,col]=size(H(:,:,i));
%     for j=1:row
%         for k=1:col
%             if H(j,k,i)==0
%                 H(j,k,i)=exp(-745);
%             end
%         end
%     end
    [Ti(:,:,i), gamma(:,:,i)]=eig(sqrtYZ);%√(YZ )的特征向量矩阵（即电流相模变换矩阵）和特征值γ
    minusGammaL(:,:,i)=-gamma(:,:,i)*l;%，即-γl，是exponet的指数，延时τ的提取要用
    Hm(:,:,i)=expm(minusGammaL(:,:,i));%得到模传播矩阵,Hm=expm(-γl),这是公式
    for j=1:row
        if Hm(j,j,i)==0
            Hm(j,j,i)=exp(-745);%避免为零，置为最小数
        end
    end
%     if i==77
%         temp=Hm(1,1,i);
%         Hm(1,1,i)=Hm(2,2,i);
%         Hm(2,2,i)=temp;
%     end
%     if i>=66
%         temp=Hm(2,2,i);
%         Hm(2,2,i)=Hm(3,3,i);
%         Hm(3,3,i)=temp;
%     end
end
%% 对H进行拟合
tic
[HFit,HmFit,tau]=optimalHFit(freq,Hm,minusGammaL,H,toleranceH,l);
if paint==1
dataPlot(H,HFit);%绘图
end
%% 对Yc进行矢量拟合
for i=1:NpMaxYc
%     option.Np=i;
%     YcFit = matrixFit(freq,Yc,option);
    opts.N=i ;%           %Order of approximation.
    opts.poletype='logcmplx'; %Mix of linearly spaced and logarithmically spaced poles
    opts.stable=1;      %翻转不稳定极点
    opts.weightparam=1; %5 --> weighting with inverse magnitude norm
    opts.Niter1=7;    %Number of iterations for fitting sum of elements (fast!)
    opts.Niter2=4;    %Number of iterations for matrix fitting
    opts.asymp=2;      %Fitting includes D
    opts.plot=0;       %不绘图
    opts.screen=0;      %拟合过程中不输出提示
    poles=[];
    [SER,rmserr,bigHfit]=VFdriver(Yc,s,poles,opts);
    YcFit.freq=freq;
    YcFit.P=SER.poles;
    YcFit.C=SER.R;
    YcFit.D=SER.D;
    YcFit.T=0;
    YcFit.MFE=[];
    YcFit.RMSE=rmserr;
    YcFit=Error(YcFit,Yc);%计算最大误差
    if YcFit.MFE<toleranceYc
        fprintf('Yc在极点数目%d时精度%f,目标精度%f\n',i,YcFit.MFE,toleranceYc);
        break
    end
end
if opts.N==NpMaxYc
    fprintf('Yc已取到最大极点数目%d，拟合精度%f，目标精度%f\n',option.Np,YcFit.MFE,toleranceH);
end
if paint==1
dataPlot(Yc,YcFit);%绘图
end
toc
%%ʸ����Ͻű�
%% ��������Yc�봫������(H��Hm)����
[H, Hm, Ti, gamma, minusGammaL,Yc]=deal(zeros(size(Y)));
traceH=zeros(1,length(freq));
omega=2*pi*freq;%�����
s=omega*1i;%s==j��
for i=1:length(s)
    sqrtYZ=sqrtm(Y(:,:,i)*Z(:,:,i));%�����YZ
    Yc(:,:,i)=sqrtYZ/Z(:,:,i);%Yc=��YZ Z^-1
    exponent=-sqrtYZ*l;%exp��A���е�A
    H(:,:,i)=expm(exponent);%�õ���������
    traceH(i)=trace(H(:,:,i));
    for j=1:length(traceH)
        if traceH(j)==0
            traceH(j)=exp(-745);
        end
    end  
    %%������ֹ����㷨�г���nan�������滻Ϊ��С��
    [row,col]=size(H(:,:,i));
%     for j=1:row
%         for k=1:col
%             if H(j,k,i)==0
%                 H(j,k,i)=exp(-745);
%             end
%         end
%     end
    [Ti(:,:,i), gamma(:,:,i)]=eig(sqrtYZ);%��(YZ )�������������󣨼�������ģ�任���󣩺�����ֵ��
    minusGammaL(:,:,i)=-gamma(:,:,i)*l;%����-��l����exponet��ָ������ʱ�ӵ���ȡҪ��
    Hm(:,:,i)=expm(minusGammaL(:,:,i));%�õ�ģ��������,Hm=expm(-��l),���ǹ�ʽ
    for j=1:row
        if Hm(j,j,i)==0
            Hm(j,j,i)=exp(-745);%����Ϊ�㣬��Ϊ��С��
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
%% ��H�������
tic
[HFit,HmFit,tau]=optimalHFit(freq,Hm,minusGammaL,H,toleranceH,l);
if paint==1
dataPlot(H,HFit);%��ͼ
end
%% ��Yc����ʸ�����
for i=1:NpMaxYc
%     option.Np=i;
%     YcFit = matrixFit(freq,Yc,option);
    opts.N=i ;%           %Order of approximation.
    opts.poletype='logcmplx'; %Mix of linearly spaced and logarithmically spaced poles
    opts.stable=1;      %��ת���ȶ�����
    opts.weightparam=1; %5 --> weighting with inverse magnitude norm
    opts.Niter1=7;    %Number of iterations for fitting sum of elements (fast!)
    opts.Niter2=4;    %Number of iterations for matrix fitting
    opts.asymp=2;      %Fitting includes D
    opts.plot=0;       %����ͼ
    opts.screen=0;      %��Ϲ����в������ʾ
    poles=[];
    [SER,rmserr,bigHfit]=VFdriver(Yc,s,poles,opts);
    YcFit.freq=freq;
    YcFit.P=SER.poles;
    YcFit.C=SER.R;
    YcFit.D=SER.D;
    YcFit.T=0;
    YcFit.MFE=[];
    YcFit.RMSE=rmserr;
    YcFit=Error(YcFit,Yc);%����������
    if YcFit.MFE<toleranceYc
        fprintf('Yc�ڼ�����Ŀ%dʱ����%f,Ŀ�꾫��%f\n',i,YcFit.MFE,toleranceYc);
        break
    end
end
if opts.N==NpMaxYc
    fprintf('Yc��ȡ����󼫵���Ŀ%d����Ͼ���%f��Ŀ�꾫��%f\n',option.Np,YcFit.MFE,toleranceH);
end
if paint==1
dataPlot(Yc,YcFit);%��ͼ
end
toc
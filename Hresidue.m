function F = Hresidue(freq,tau,Poles,Data,doMOR)
%X = Hresidue(freq,tau,Poles,data)
%   ���������ĺ���
%-------------------���-------------------
%   F�ǽṹ�壺
%         F.freq �洢Ƶ��
%         F.P �洢��������
%         F.C �洢��������
%         F.T
%         �洢��ʱ���������Ѷ��������˴�����F.C��ά������ͬ��ʹ��F.C��:,:,i����Ӧ��F.T(i).����������F.T��Ϊ�����㣨matrixFit��Ϊ�����㣩
%         F.MFE����������maxFittingError������������ֵ����ǲ����֣�����ֻ�з�ֵ
%         F.RMSE�����ľ�����������������ֻ���ֵ
%         F.MRPR�������������ȣ�����������֮��ȡ����ֵ
%   ----------------����---------------------------      
%   freq��Ƶ�ʣ�Ns���ȣ���Ҫ��������������ǣ����Զ�תΪ������������s
%   tau����ʱ��������������ν��Nk����
%   Poles��=Ԫ��=�������Ҫ��ÿһ�ж�Ӧһ����,ʱ������Ϊ numel��Poles����������Ŀ
%   Data������,��ά���� Nc*Nc*Ns
%   doMOR:�Ƿ���н��ף�0���������ף�1��������

%% ���봦��
%%�Լ��㴦��
% [Nk Np]=size(Poles);
Nk=numel(Poles);%��ȡ��ʱ����Ŀ
% poles=reshape(Poles.',1,Nk*Np);%�����м�����ϳ�һ��������������ʱ˳������[p11,p12,p21,p22,p31,p32]
poles=[];
for i=1:Nk
    Np(i)=length(Poles{i});%��ʱ���Ӧ�ļ�����Ŀ
    poles=[poles,reshape(Poles{i},1,Np(i))];%�����м�����ϳ�һ��������������ʱ˳������[p11,p12,p21,p22,p31,p32]
end
Npp=length(poles);%��������Ŀ
%% ��Ǹ�������
cindex=zeros(size(poles));%ʵ�����Ϊ�㣬�������Ϊ1���乲���������һλ����Ϊ2
%���ӣ�         poles=��[1 0 i 1 0 i]��,
%���cindex����Ϊ7��    [0 0 1 2 0 1 2]
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
%%��Ƶ�ʴ���
[Nc,Nc,Ns]=size(Data);
Ns=length(freq);
freq=reshape(freq,Ns,1);%��Ƶ��תΪ������
s=2i*pi*freq;%�õ���������s
%%��Data������ά������ά
%%
symmetry=0;%Ĭ��Ϊ���Գ�
if isequal(Data(:,:,fix(Ns/2)),Data(:,:,fix(Ns/2)).')
    symmetry=1;
    Ne=(Nc+1)*Nc/2;
    Data1=zeros(1,Ne,Ns);
    map=logical(triu(ones(Nc)));
    for n=1:Ns
        a=Data(:,:,n);
        Data1(:,:,n)=reshape(a(map),1,Ne);
    end  %Data1ҳƵ�ʣ���Ԫ�أ�һ��
else
    Ne=Nc*Nc;
    for n=1:Ns
        Data1(:,:,n)=reshape(Data(:,:,n),1,Ne);
    end  %Data1ҳƵ�ʣ���Ԫ�أ�һ��
end
Data2=shiftdim(Data1,1);%Data2 Ne��Ԫ�أ�Ns��Ƶ��
Response=Data2.';%���¹������Ӧ����Ns by Ne
%% ��Ȩ
% Weight=sqrt((Response));
Weight=ones(Ns,Ne);
% Weight=repmat(Weight,1,Ne);
% Response=Response.*Weight;%��Ӧ��Ȩ
%% ��������
%%������
F.freq=freq;%������
F.P=poles;%����������
F.C=zeros(Ne,1,Npp);%���������ÿһҳ��Ӧһ������
F.D=zeros(Nc,Nc);%���ڱ�����ר��������H���������ʴ���Ϊ�����
F.T=[];%�ҵ��������ÿ�����������Ӧһ���ӣ�����������Npp=Np(1)+Np(2)+Np(3)����

t=[];%t�������Ŧ�
for i=1:Nk%��ʱ�ӵ���Ŀ,%��֪��ÿ��ģ���ļ�����Ŀ-----------------------------------
    t=[t,repmat(tau(i),1,Np(i))];%T��������exp��-s�ӣ�
end
F.T=t;
%2019/3/4
if 0
%% ��ϵ�������� AAx=bb
rowNum=min(2*Ns+1-Npp,Npp+1);%���ǵ�Ƿ�������������ȡ������С
AA=zeros(Ne*rowNum,Npp+1);
bb=zeros(Ne*rowNum,1);
scale=0;%relax����
for m=1:Ne
    scale=scale+(norm(Response(:,m)))^2;%ÿ������һ��Ȩ������ 
end
scale=sqrt(scale)./Ns; %������Ԫ�ص�ƽ���ͽ��п�������ƽ��
Nimax=4;%����������
tolerance=0.5;%�����жϼ����Ƿ�����
x=0;
for i=1:Nimax
    %%����ϵ������,ÿ��Ԫ�ض���ͬ
    H=ones(Ns,2*Npp+1);%�м����H���������h��s��
    P=repmat(s,1,2*Npp);%�м����P���������1/��s-p��
    A=repmat(poles,Ns,2);%����һ��A����
    P=1./(P-A);%����� 1/��s-p��
    P(:,end+1)=ones(Ns,1);%relax��������
    for m=1:length(poles)
        if cindex(m)==1
            P1=P(:,m);
            P2=P(:,m+1);
            P(:,m)=P1+P2;%
            P(:,m+1)=(P1-P2)*1i;
            P(:,m+length(poles))=P1+P2;
            P(:,m+1+length(poles))=(P1-P2)*1i;
        end
    end%ʹ����������ȱ��ʵ��
    for n=1:Ne%��������AA
        H(:,Npp+1:end)=repmat(-1.*Response(:,n),1,Npp+1);
        A=H.*P;
        A(:,1:Npp)=A(:,1:Npp).*exp(-s*t);%�õ�ϵ������Ns by Npp, ����� exp(-s��)/��s-p��
        [~,Acol]=size(A);
        % W=repmat(Weight(:,1),1,Acol);
        % A=A.*W;%�����ÿ��Ԫ��Ȩ�ز�ͬ�Ļ���û����һ��Ԫ�أ��͵ø�Ȩ��
        A=[real(A);imag(A)];
        A(end+1,:)=[zeros(1,Npp), real(sum(P (:,Npp+1:end)))*scale];
        [Q R]=qr(A,0);%QR�ֽ⣬��ȡ�������صĲ���
        R22=R(Npp+1:end,Npp+1:end);
        AA((n-1)*rowNum+1:n*rowNum,:)=R22;
        bb((n-1)*rowNum+1:n*rowNum,1)=Q(end,Npp+1:end)'*Ns*scale;
    end
    %% ����AA����������
    [rowOfAA colOfAA]=size(AA);
    col=1:colOfAA;%AA���������
    rnormArray=arrayfun(@(col) norm(AA(:,col),2),col);%%������������rnormArray���AAÿһ��������2����
    rnormArray=1./(rnormArray);%�õ�����ĵ���
    CondScale=repmat(rnormArray,rowOfAA,1);%��������ž���
    AA=AA.*CondScale;%%�������
    %% �����С���˷���
    x=AA\bb;%����˦ҵ�����,[c1~,c2~,...]
    %% ��ԭAA����������֮ǰ�Ľ�
    x=x.*rnormArray.';
    if x(end)==0
        x(end)=1;
    end
    %% �����¼���
    X=repmat(x(1:end-1),1,Npp);%����ֻ��Ҫ�������ֵĴ�С
    lambda=diag(poles);
    %��lambda�еĸ���������޸�
    for m=1:length(poles)
        if cindex(m)==1
            lambda(m:m+1,m:m+1)=[real(poles(m)),imag(poles(m));-imag(poles(m)),real(poles(m))];
            X(m:m+1,m:m+1)=[2*x(m),2*x(m+1);0,0];
        end
    end
    D=lambda-X./x(end);%��ʽ�е�d~^-1
    newPoles=eig(D).';
    %%��ת���ȶ�����
    unstables=real(newPoles)>0;%�洢���ȶ����������
    if sum(unstables)>0%˵���в��ȶ�����
        disp(['���ڲ��ȶ����㣬����ת'])
        newPoles(unstables)=newPoles(unstables)-2*real(newPoles(unstables));%��ת���ȶ�����
    end
    %%�ж�����
    if norm(newPoles-poles)/norm(poles)<=tolerance
        %         disp(['���������� tolerance=' num2str(tolerance)])
        break
    else
        poles=newPoles;
    end
end
end
%% ����С����A
%%����ϵ������,ÿ��Ԫ�ض���ͬ
H=ones(Ns,Npp);%�м����H���������h��s��
P=repmat(s,1,Npp);%�м����P���������1/��s-p��
A=repmat(poles,Ns,1);%����һ��A����
P=1./(P-A);%����� 1/��s-p��
% A=P.*exp(-s*t);
for m=1:length(poles)
    if cindex(m)==1
        P1=P(:,m);
        P2=P(:,m+1);
        P(:,m)=P1+P2;%
        P(:,m+1)=(P1-P2)*1i;
    end
end%ʹ����������ȱ��ʵ��
%%2019/3/7
A=P.*exp(-s*t);%�õ�ϵ������Ns by Npp, ����� exp(-s��)/��s-p��
save math Response
% W=repmat(Weight(:,i),1,Acol);
% A=A.*W;%�����ÿ��Ԫ��Ȩ�ز�ͬ�Ļ���û����һ��Ԫ�أ��͵ø�Ȩ��
A=[real(A);imag(A)];
%% ����A����������
[rowOfA colOfA]=size(A);
col=1:colOfA;%A���������
rnormArray=arrayfun(@(col) norm(A(:,col),2),col);%%������������rnormArray���Aÿһ��������2����
rnormArray=1./(rnormArray);%�õ�����ĵ���
CondScale=repmat(rnormArray,rowOfA,1);%��������ž���
A=A.*CondScale;%%�������
%% �ֱ����ÿ��Ԫ�ص�����
for n=1:Ne
    b=Response(:,n).*Weight(:,n);%��n��Ԫ�ص���Ӧ����
    b=[real(b);imag(b)];
    %%���
    x=A\b;%�������,������,[c11,c12,...c21,c22,...]
    x=x.*rnormArray.';%%��ԭA����������֮ǰ�Ľ�  
    for m=1:length(poles)
        if cindex(m)==1
            x1=x(m);
            x2=x(m+1);
            x(m)=x1+x2*1i;
            x(m+1)=x1-x2*1i;
        end
    end%������������Ϊ����
    F.C(n,1,:)=x(1:end);
    MRPRatio(n)=max(abs(x)./abs(poles).');%������������(Max Residue/Pole Ratio),x����������poles��������������MRPRatioҲ��������
end
F.MRPR=max(MRPRatio);%��Ϊ�������
%% �ж��Ƿ���н��ס�MOR����,�ɷ�ֹ��������������
if doMOR
    disp(['��������ȳ�����ֵ�����н���']);
    [Srow,Scol]=size(A);
    Sindex=min(Srow,Scol);
    global MRPRtolerance %MRPRtolerance��fittingSet����
    while F.MRPR>MRPRtolerance
        [U,S,V]=svd(A);
        save svdResidue U S V
        save datA         
        S(Sindex,Sindex)=0;
        Sindex=Sindex-1;
        A=U*S*V';
        %% �ֱ����ÿ��Ԫ�ص�����
        for n=1:Ne
            b=Response(:,n);%��n��Ԫ�ص���Ӧ����
            b=[real(b);imag(b)];
            %%���
            x=A\b;%�������,������,[c11,c12,...c21,c22,...]
            x=x.*rnormArray.';%%��ԭA����������֮ǰ�Ľ�
            %     AAA*x-b<=eps
            
            for m=1:length(poles)
                if cindex(m)==1
                    x1=x(m);
                    x2=x(m+1);
                    x(m)=x1+x2*1i;
                    x(m+1)=x1-x2*1i;
                end
            end%������������Ϊ����
            F.C(n,1,:)=x(1:end);
            MRPRatio(n)=max(abs(x)./abs(poles).');%������������(Max Residue/Pole Ratio),x����������poles��������������MRPRatioҲ��������
        end
        F.MRPR=max(MRPRatio);%��Ϊ�������
    end
    %% ��ȡʣ�༫��
    index=(F.C(1,1,:)~=0);%ֻ�����Ͻǵ�Ԫ�ؼ���
    F.C=F.C(:,:,index);
    F.T=F.T(index);
    F.P=F.P(index);
end%MOR
Npp=length(F.P);
%% �������ָ���ԭ�о�����ʽ
if symmetry%����ǶԳ�
    index=1;%�������ÿһ����������ʼλ��
    C=zeros(Nc,Nc,Npp);%���������������----------------------?Ϊʲô��Np��������Ӧ��������Ŀ
    %%������������
    for n=1:Nc
        len=n;%length����matlab���к�������ͻ���ʸ�Ϊlen
        C(1:n,n,:)=F.C(index:index+len-1,1,:);
        index=index+len;      
    end
    %%�����������ԭ����
    i=1:Npp;%���޸�Ϊ�ܼ�����Ŀ
    cellC=arrayfun(@(a) triu(C(:,:,a),1).'+C(:,:,a),i,'UniformOutput',false);%cellC��һ��Ԫ�����飬ÿһ��Ԫ���洢һ�������Ӧ����������
    F.C=cat(3,cellC{:});%��Ԫ������ת��Ϊ��ά����  
else%����ǷǶԳ�
    F.C=reshape(F.C,Nc,Nc,Npp);%���޸�Ϊ����Ŀ
end
%% ����������ݡ�RMSE��MFE��maxFittingError�����㣬������ֲ��dataPlot
%���������������
% Np=length(F.P);%ǰ������õ�Np��ÿ��ģ���ģ������ݵ�Np��һ���ļ�����Ŀ
f=zeros(size(Data));
for i=1:Ns
    for j=1:Npp
        f(:,:,i)=f(:,:,i)+exp(-s(i)*F.T(j)).*F.C(:,:,j)./(s(i)-F.P(j));
    end
    f(:,:,i)=f(:,:,i)+F.D; 
end
%����ֵ�ĺ���
% figure;
% for i=1:Nc
%     for j=1:Nc
%         semilogx(F.freq,abs(squeeze(Data(i,j,:))),'b'); 
%         hold on
%         semilogx(F.freq,abs(squeeze(f(i,j,:))),'r.');
%     end
% end
% hold off
% xlabel('Ƶ�� Hz')
% ylabel('��ֵ ')
% legend('ԭ����','��Ͻ��')
%% �����ĺ���
% figure;
% denominator=max(max(abs(Data(i,j,:)))
for i=1:Nc
    for j=1:Nc
        %����Ҫ����RMSE
        residualError=abs(squeeze(Data(i,j,:)-f(i,j,:)));%ÿ����Ĳв�
%         residualError=abs(squeeze(abs(Data(i,j,:))-abs(f(i,j,:))));%ÿ����Ĳв��ģֵ�Ĳв�
        MFE(i,j)=max(residualError)/max(max(max(abs(squeeze(Data)),[],3)));%�ҳ��������������,����pscad����
        SSR=sum(residualError.^2);%�в�ƽ����
        RMSE(i,j)=sqrt(SSR/length(residualError));%��������������
%         semilogx(F.freq,residualError,'r'); %%�в�
        %����Ҫ�������������MFE
%         hold on
    end
end
%�ҳ��������������������Ͻṹ��
F.MFE=max(max(MFE));
F.RMSE=max(max(RMSE));%����Ԫ�����ľ������������ṹ��
% F.MFE=MFE;
% F.RMSE=RMSE;%����Ԫ�����ľ������������ṹ��
% title('��ֵ��ϲв�')
% xlabel('Ƶ�� Hz')
% ylabel('��ֵ ')
% hold off
 %% ����ǵĺ���
% figure;
% for i=1:Nc
%     for j=1:Nc
%         semilogx(F.freq,180*unwrap(angle(squeeze(Data(i,j,:))))/pi,'b'); 
%         hold on
%         semilogx(F.freq,180*unwrap(angle(squeeze(f(i,j,:))))/pi,'r.');
%     end
% end
% hold off
% xlabel('Ƶ�� Hz')
% ylabel('��� ��')
% legend('ԭ����','��Ͻ��')
end


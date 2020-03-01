function [F] = traceFit(frequency,tau,Data,option)
%[F] = traceFit(frequency,tau,Data,option)
%   ���룺
%   frequency
%   tau����ʱ��ȡ����,���ѷֺ���
%   Data������ļ�����ô��һ�������������ǵĻ�������ת������������
%   option�ǿ��ƽṹ�壬���������ֶ�
%       Np          ָ��������Ŀ,����С�ڲ���������frequency�ĳ��ȣ�
%       passivity   �Ƿ������Դ�Լ�飬��ǿ����Դ��1=YES��0=No��
%       weight      ��С�������Ȩ��(0=�����κ�Ȩ�أ�1=��߻�Ƶ����Ӧ����)
%       realFit     �Ƿ����ʵ������ϣ�1=Yes��0=No��
%       relax       �Ƿ���ÿ��ɷ�����1=yes��0=no��
%   �����
%   F��һ�ṹ��
%       F.T���յ���ʱ
%       F.P���յļ���
%       F.C����
%       F.D������
%       F.MFE���������
%       F.RMSE���������

%% ���ݴ���
[Nc Nc Ns]=size(Data);
Ns=length(frequency);
%%��֤frequency��������
[row column]=size(frequency);
if row <column
    frequency=frequency.';
end
Response=reshape(Data,length(Data),1);%���¹������Ӧ���󣬱�֤Ϊ������
s=2i*pi*frequency;%��������Ns by 1
%% ȷ������
Np=option.Np;
% disp(['������ĿΪ ' num2str(Np)])
poles=initialPoles(frequency,Np,1);%2��ʾ����ʵ���㣬1������
%% �γɼ����õ���ʱ����
Ng=length(tau);
t=[];
for i=1:Ng
    t=[t,repmat(tau(i),1,Np)];
end
%% ��ͨ��С����
%% ȷ��Ȩ��
switch option.weight
    case 0
        Weight=ones(Ns,1);%������Ȩ���ټ�
    case 1
        for i=1:length(frequency)%pscadȨ�أ���50Hz��1000��������һ����Ϊ������������
            if frequency(i)>50&frequency(i-1)<=50
                Weight(i-1)=1000;
                Weight(i)=1;
            else
                Weight(i)=1;
            end
        end
        Weight=reshape(Weight,Ns,1);%�������ɵ�����������������Ҫ�����������
end
% disp(['��С����Ȩ��Ϊ1'])
%% ��ϵ�������� AAx=bb
if option.relax
   AA=zeros(Np+1);
   bb=zeros(Np+1,1);
else
    AA=zeros(Np);
    bb=zeros(Np,1);
end
b=Weight.*Response;%����Ӧ��Ȩ
b=[real(b);imag(b)];
scale=(norm(b))^2;%ÿ������һ��Ȩ��������relax���� 
scale=sqrt(scale)./Ns; %������Ԫ�ص�ƽ���ͽ��п�������ƽ��
%% �������
Nimax=10;
tolerance=0.5;%�����жϼ����Ƿ�����
x=0;
for i=1:Nimax
    disp(['������ ' num2str(i) '��'])
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
    %%����ϵ������
    Ncol=2*Np;%Сϵ������A������
    Ncol1=Np;%Сϵ������Aǰ������޹صĲ���
    
    H=ones(Ns,Ncol);%�м����H���������h��s��
    P=ones(Ns,Ncol);%�м����P���������1/��s-p��
    A=zeros(Ns,Ncol);%Сϵ������A
    
    P=repmat(s,1,2*Np);
    A=repmat(poles,Ns,2);%����һ��A����
    P=1./(P-A);
    if option.relax
        P(:,end+1)=ones(Ns,1);
    end
    %%�Ը��������޸�
    for m=1:length(poles)
        if cindex(m)==1
            P1=P(:,m);
            P2=P(:,m+1);
            P(:,m)=P1+P2;%
            P(:,m+1)=(P1-P2)*1i;
            P(:,m+length(poles))=P1+P2;
            P(:,m+1+length(poles))=(P1-P2)*1i;
        end
    end%��Ӧ�����¼���ҲҪ�޸�
    H(:,Ncol1+1:end)=repmat(-1.*Response,1,Np);
    if option.relax
        H(:,end+1)=-1.*Response;
    end
    A=H.*P;%%�õ���Сϵ������
    %��A����Ϊ����ϵ���ʽ
    delayFactorMatrix=s*t;%�������һ������
    Atemp=repmat(A(:,1:Ncol1),1,Ng).*exp(-delayFactorMatrix);
    A=[Atemp,A(:,Ncol1+1:end)];%��ʱAΪ����ϵ���ȷϵ������
    %%��ϵ������A��Ȩ
    [~,colA]=size(A);
    W=repmat(Weight,1,colA);
    A=W.*A;%��Ȩ
    %��ʵ�ֽ��Ա�֤���b����ѭ����ֽ����
    A=[real(A);imag(A)];
    Ncol2=Np*Ng;
    Ncol3=Np*(Ng+1);
    if option.relax
        A(end+1,:)=[zeros(1,Ncol2), real(sum(P (:,Ncol1+1:end)))*scale];%Ncol1=Np
        [Q R]=qr(A,0);%QR�ֽ⣬��ȡ�������صĲ���
        R22=R(Ncol2+1:end,Ncol2+1:end);
        AA(1:Np+1,:)=R22;
        %             size(bb((n-1)*(Np+1)+1:n*(Np+1),1))
        %             size(Q(end,Ncol1+1:end)'*Ns*scale)
        bb(1:Np+1,1)=Q(end,Ncol2+1:end)'*Ns*scale;
    else
        [Q R]=qr(A,0);%QR�ֽ⣬��ȡ�������صĲ���
        R22=R(Ncol2+1:Ncol2+Np,Ncol2+1:Ncol2+Np);
        AA(1:Np,:)=R22;
        bb(1:Np,1)=Q(:,Ncol2+1:Ncol2+Np)'*b;
    end
    %�����relax����������Ҫ�����޸�
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
    if sum(abs(x)==inf)>=1
        x((abs(x)==inf))=sign(x((abs(x)==inf)))*realmax;%���ĳһ�������˻������ֵ���ͽ���ְλ���ֵ������inf���ں�������¼���eig�б���
    end
    if option.relax
        if x(end)==0
            x(end)=1;
            %         elseif x(end)<Tol
        end
    else
        x(end+1)=1;%Ӧ��ΪDnew
    end
    %% �����¼���
    X=repmat(x(1:end-1),1,Np);%����ֻ��Ҫ�������ֵĴ�С
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
end%�������
if i==Nimax
    disp(['ʸ����ϴﵽ����������',num2str(Nimax)])
    %Np=Np+2;
    %����Ӧ�÷������¼��㣬���Ȳ�д
end
%% �������
Ne=1;%Ϊ�˼�㣬�����Ne�Ҳ�ɾ�ˣ��ڼ������������1
option.pattern=0;%ͬ��Ϊ�˼�㣬ֻ�ǶԴ��������ļ���ϣ��޳�����
%%����Ƶ��
F.freq=frequency;%������
F.T=t;
%%������
F.C=zeros(Ne,1,Np*Ng);%���������ÿһҳ��Ӧһ������
F.D=zeros(Ne,1);%��ų�����
%%����
F.P=repmat(poles,1,Ng);%������
%%��������,ÿ��Ԫ�ص�ϵ������һ��
P=ones(Ns,Ncol1);%�м����P���������1/��s-p��
A=zeros(Ns,Ncol1);
P(:,1:end)=repmat(s,1,Np);
A(:,1:end)=repmat(poles,Ns,1);%����һ��A����
P=1./(P-A);
for m=1:length(poles)
    if cindex(m)==1
        P1=P(:,m);
        P2=P(:,m+1);
        P(:,m)=P1+P2;%
        P(:,m+1)=(P1-P2)*1i;
    end
end%ʹ����������ȱ��ʵ��
A=P;%�õ�ϵ������Ns by Np
%��A����Ϊ����ϵ���ʽ
A=repmat(A(:,1:end),1,Ng).*exp(-delayFactorMatrix);%��ʱAΪ����ϵ���ȷϵ������
[~,colA]=size(A);
W=repmat(Weight,1,colA);
A=W.*A;%��Ȩ
A=[real(A);imag(A)];
%b=Weight.*Response;ǰ������Ӧ��Ȩ
%% ����A����������
[rowOfA colOfA]=size(A);
col=1:colOfA;%A���������
rnormArray=arrayfun(@(col) norm(A(:,col),2),col);%%������������rnormArray���Aÿһ��������2����
rnormArray=1./(rnormArray);%�õ�����ĵ���
CondScale=repmat(rnormArray,rowOfA,1);%��������ž���
A=A.*CondScale;%%�������
%% MOR����,ɸѡ���Ӽ��㼯
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
for n=1:Ne %�ֱ����ÿ��Ԫ�ص�����
    x=A\b(:,n);%�������,������,[d,c1,c2,...]
    x=x.*rnormArray.';%�ָ��������Ż�֮ǰ�Ľ�
    for m=1:length(poles)
        if cindex(m)==1
            x1=x(m+option.pattern);
            x2=x(m+option.pattern+1);
            x(m+option.pattern)=x1+x2*1i;
            x(m+option.pattern+1)=x1-x2*1i;
        end
    end%������������Ϊ����
    if option.pattern==1
        F.D(n,1)=x(1);%��Yc��˵����һ����d
        F.C(n,1,:)=x(2:end);
    else
        F.D(n,1)=0;%��Hmps��˵����d
        F.C(n,1,:)=x(1:end);
    end
%     F.D(n,1)=abs(F.D(n,1))*real(F.D(n,1))/abs(real(F.D(n,1)));%��������ʵ�ֽ⣬ֻ�Ǽ򵥵�ȡʵ�����ᶪʧ��λ����Ϣ
end
%% ��ȡ�Ӽ���
index=(F.C~=0);
F.C=F.C(index);
F.T=F.T(index);
F.P=F.P(index);
% size(F.C)
%% ����������ݡ�RMSE��MFE��maxFittingError�����㣬������ֲ��dataPlot
%���������������
f=zeros(size(Data));
for i=1:Ns
    for j=1:length(F.P)
        f(i)=f(i)+F.C(:,:,j)./(s(i)-F.P(j));
    end
    f(i)=f(i)+F.D; 
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
%����Ҫ����RMSE
residualError=abs(Data-f);%ÿ����Ĳв�
F.MFE=max(residualError)/max(abs(Data));%�ҳ��������������,����pscad����
SSR=sum(residualError.^2);%�в�ƽ����
F.RMSE=sqrt(SSR/length(residualError));%��������������
%         semilogx(F.freq,residualError,'r'); %%�в�
%����Ҫ�������������MFE
%         hold on


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



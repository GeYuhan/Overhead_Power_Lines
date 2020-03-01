function [F] = matrixFit(frequency,Data,option)
%   F = matrixFit(frequency,Data)
%Ҫ������frequency����������Ns by 1����λHz��������������������Զ�ת��Ϊ������
%Data����ά����Nc by Nc by Ns�������Ǹ�����
%option�ǿ��ƽṹ�壬���������ֶ�
%	pattern     ѡ�������ʽ��Yc/HmMPS��
%	Np          ָ��������Ŀ,����С�ڲ���������frequency�ĳ��ȣ�
%	passivity   �Ƿ������Դ�Լ�飬��ǿ����Դ��1=YES��0=No��
%	weight      ��С�������Ȩ��(0=�����κ�Ȩ�أ�1=��߻�Ƶ����Ӧ����)
%	realFit     �Ƿ����ʵ������ϣ�1=Yes��0=No��
%   relax       �Ƿ���ÿ��ɷ�����1=yes��0=no��
%
%��������ϸ���ͣ�
%option.pattern����ѡ����Yc���ݣ�����Hmps����
%��option.pattern==1,��Yc����
%       F��������ʽ��ϣ�
%              N     Cn
%       F(s)= SUM  -------  + D 
%             n=1   s-pn
%      
%��option.pattern==0����Hmps����
%       F��������ʽ��ϣ�
%              N     Cn
%       F(s)= SUM  -------
%             n=1   s-pn
%      
%F�ǽṹ��
%   F.freq==frequency��Ƶ�������� 1 Ns
%   F.P�Ǽ�����������1 by Np
%   F.C������������άNc Nc Np��ÿһҳ��Ӧһ������
%   F.D�ǳ��������Nc Nc��HΪ��
%   F.T����ʱ������Nk*Np��YcΪ�㣬HΪ���㡣���ڱ������Ƕ�Yc����С��λHmMPS��ϣ��ʾ�Ϊ�㡣
%   F.MFE����������maxFittingError������������ֵ����ǲ����֣�����ֻ�з�ֵ
%   F.RMSE�����ľ�����������������ֻ���ֵ
%
%����������е����ã�
%������Ŀ��ǰ��Ϊָ�� Np
%������������Ϊ�趨 Nimax�������о� tolerance=0.01
%û��д��Ͼ����жϣ��������Ӽ�����Ŀ����������ʵ�֡�
%ʹ��relax����
%ʹ���������Ż���ʹ����ʵ�ֽ�
%��С����Ȩ�ػ�δָ������ǰΪ1 weight���Ѹ�ΪpscadȨ�ء�
%
%% ���ݴ���
[Nc Nc Ns]=size(Data);
Ns=length(frequency);
%%��֤frequency��������
[row column]=size(frequency);
if row <column
    frequency=frequency.';
end
%%��ά������ά
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
s=2i*pi*frequency;%��������Ns by 1
%% ȷ������
Np=option.Np;
% disp(['������ĿΪ ' num2str(Np)])
poles=initialPoles(frequency,Np,1);%2��ʾ����ʵ���㣬1������
%% ��ͨ��С����
%% ȷ��Ȩ��
switch option.weight
    case 0
        Weight=ones(Ns,Ne);%������Ȩ���ټ�
    case 1
        for i=1:length(frequency)%pscadȨ�أ���50Hz��1000��������һ����Ϊ������������
            if frequency(i)>50&frequency(i-1)<=50
                Weight(i-1)=1000;
                Weight(i)=1;
            else
                Weight(i)=1;
            end
        end
        Weight=repmat(Weight.',1,Ne);%�������ɵ�����������������Ҫ�����������
end
% disp(['��С����Ȩ��Ϊ1'])
%% ��ϵ�������� AAx=bb
if option.relax
   AA=zeros(Ne*(Np+1),Np+1);
   bb=zeros(Ne*(Np+1),1);
else
    AA=zeros(Ne*Np,Np);
    bb=zeros(Ne*Np,1);
end
b=Weight.*Response;%����Ӧ��Ȩ
b=[real(b);imag(b)];
scale=0;%relax����
for m=1:Ne
    scale=scale+(norm(b(:,m)))^2;%ÿ������һ��Ȩ������ 
end
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
    %% ��ϵ������AA��bb�Ĺ���
    for n=1:Ne
        Ncol=2*Np+option.pattern;%Сϵ������A������
        Ncol1=Np+option.pattern;%Сϵ������Aǰ������޹صĲ���
    
        H=ones(Ns,Ncol);%�м����H���������h��s��
        P=ones(Ns,Ncol);%�м����P���������1/��s-p��
        A=zeros(Ns,Ncol);%Сϵ������A

        P(:,1+option.pattern:end)=repmat(s,1,2*Np);
        A(:,1+option.pattern:end)=repmat(poles,Ns,2);%����һ��A����
        P=1./(P-A);
        if option.relax
            P(:,end+1)=ones(Ns,1);
        end
        %%�Ը��������޸�
        for m=1:length(poles)
            if cindex(m)==1
                P1=P(:,m+option.pattern);
                P2=P(:,m+option.pattern+1);
                P(:,m+option.pattern)=P1+P2;%
                P(:,m+option.pattern+1)=(P1-P2)*1i;
                P(:,m+option.pattern+length(poles))=P1+P2;
                P(:,m+option.pattern+1+length(poles))=(P1-P2)*1i;
            end
        end%��Ӧ�����¼���ҲҪ�޸�
        H(:,Ncol1+1:end)=repmat(-1.*Response(:,n),1,Np);
        if option.relax
            H(:,end+1)=-1.*Response(:,n);
        end
        A=H.*P;%%�õ���Сϵ������
        [~,colA]=size(A);
        W=repmat(Weight(:,n),1,colA);
        A=W.*A;%��Ȩ
        %��ʵ�ֽ��Ա�֤���b����ѭ����ֽ����
        A=[real(A);imag(A)];
        if option.relax
            A(end+1,:)=[zeros(1,Ncol1), real(sum(P (:,Ncol1+1:end)))*scale];
            [Q R]=qr(A,0);%QR�ֽ⣬��ȡ�������صĲ���
            R22=R(Ncol1+1:end,Ncol1+1:end);
            AA((n-1)*(Np+1)+1:n*(Np+1),:)=R22;
%             size(bb((n-1)*(Np+1)+1:n*(Np+1),1))
%             size(Q(end,Ncol1+1:end)'*Ns*scale)
            bb((n-1)*(Np+1)+1:n*(Np+1),1)=Q(end,Ncol1+1:end)'*Ns*scale;
        else
            [Q R]=qr(A,0);%QR�ֽ⣬��ȡ�������صĲ���
            R22=R(Ncol1+1:Ncol1+Np,Ncol1+1:Ncol1+Np);
            AA((n-1)*Np+1:n*Np,:)=R22;
            bb((n-1)*Np+1:n*Np,1)=Q(:,Ncol1+1:Ncol1+Np)'*b(:,n);
        end
        %�����relax����������Ҫ�����޸�
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
%%����Ƶ��
F.freq=frequency;%������
F.T=0;%����ʱ����Ϊ�㣬���治����ᡣ
%%������
F.C=zeros(Ne,1,Np);%���������ÿһҳ��Ӧһ������
F.D=zeros(Ne,1);%��ų�����
%%����
F.P=poles;%������
%%��������,ÿ��Ԫ�ص�ϵ������һ��
P=ones(Ns,Ncol1);%�м����P���������1/��s-p��
A=zeros(Ns,Ncol1);
P(:,1+option.pattern:end)=repmat(s,1,Np);
A(:,1+option.pattern:end)=repmat(poles,Ns,1);%����һ��A����
P=1./(P-A);
for m=1:length(poles)
    if cindex(m)==1
        P1=P(:,m+option.pattern);
        P2=P(:,m+option.pattern+1);
        P(:,m+option.pattern)=P1+P2;%
        P(:,m+option.pattern+1)=(P1-P2)*1i;
    end
end%ʹ����������ȱ��ʵ��
A=P;%�õ�ϵ������Ns by Np+1
W=repmat(Weight(:,n),1,Ncol1);%����n=Ne���������ĳ��ѭ���������µĽ��������Ӧ���ǲ�ͬ��Ԫ�ز�ͬ��Ȩ�أ�����ʵ�ʲ�����ÿ��Ԫ��Ȩ�ض�һ����Ҳ������ν�����
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
%% �������ָ���ԭ�о�����ʽ
if symmetry%����ǶԳ�
    index=1;%�������ÿһ����������ʼλ��
    C=zeros(Nc,Nc,Np);%���������������
    D=zeros(Nc,Nc);
    %%������������
    for n=1:Nc
        len=Nc-n+1;%length����matlab���к�������ͻ���ʸ�Ϊlen
        C(n,n:end,:)=F.C(index:index+len-1,1,:);
        D(n,n:end)=F.D(index:index+len-1,1);
        index=index+len;      
    end
    %%�����������ԭ����
    i=1:Np;
    cellC=arrayfun(@(a) triu(C(:,:,a),1).'+C(:,:,a),i,'UniformOutput',false);%cellC��һ��Ԫ�����飬ÿһ��Ԫ���洢һ�������Ӧ����������
    F.C=cat(3,cellC{:});%��Ԫ������ת��Ϊ��ά����
    F.D=triu(D,1).'+D;%��D���������ת�ò���ӵõ����ĳ����������ΪD�Ƕ�ά���󣬹�û��C��ô�鷳   
else%����ǷǶԳ�
    F.D=reshape(F.D,Nc,Nc);
    F.C=reshape(F.C,Nc,Nc,Np);
end

%% ����������ݡ�RMSE��MFE��maxFittingError�����㣬������ֲ��dataPlot
%���������������
f=zeros(size(Data));
for i=1:Ns
    for j=1:Np
        f(:,:,i)=f(:,:,i)+F.C(:,:,j)./(s(i)-F.P(j));
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
for i=1:Nc
    for j=1:Nc
        %����Ҫ����RMSE
        residualError=abs(squeeze(Data(i,j,:)-f(i,j,:)));%ÿ����Ĳв�
        MFE(i,j)=max(residualError)/max(abs(Data(i,j,:)));%�ҳ��������������,����pscad����
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
    

    











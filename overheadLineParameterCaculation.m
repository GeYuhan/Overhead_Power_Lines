function [Z,Y]=overheadLineParameterCaculation
%% ���ؼܿ�������
overheadLineParameter
global freq
%% ��������
%%����ӵ�����Ŀ�����Բ�뾶
numOfLines=N*numOfPhase;%�������Ŀ
numOfConductors=numOfLines+numOfGroundLines;%������Ŀ���������������ߣ�
radiusOfCircle=d/2/sin(pi/N);%����Բ�İ뾶����λm
%%������ѵ����ӵ���������ĵ�λ��
XrelativeToPhase=cell2mat(arrayfun(@(x) radiusOfCircle*cos(pi/4*(x-1)),1:N,'UniformOutput', false));%������
YrelativeToPhase=cell2mat(arrayfun(@(x) radiusOfCircle*sin(pi/4*(x-1)) ,1:N,'UniformOutput',false));
%%������е��壨���������ߣ���λ������
XcoordinateOfConductor=ones(N,1)*XcoordinateOfPhase+repmat(XrelativeToPhase.',1,numOfPhase);%����,numOfphase����3��
XcoordinateOfConductor=[reshape(XcoordinateOfConductor,1,numOfLines),XcoordinateOfGround];%������,������ĩβ
YcoordinateOfConductor=ones(N,1)*YcoordinateOfPhase+repmat(YrelativeToPhase.',1,numOfPhase);
YcoordinateOfConductor=[reshape(YcoordinateOfConductor,1,numOfLines),YcoordinateOfGround];%������,������ĩβ
%%��������ߣ��ɽ�����λ�ü�������ĩβ���������� numOfLines+numOfGroundLines
%%�뾶������絼���������ɣ���Z�������ʱ��
radiusOfConductor=[repmat(radiusOfLine,1,numOfLines),repmat(radiusOfGroundLine,1,numOfGroundLines)];%������е���İ뾶��������
conductivityOfConductor=[repmat(conductivityOfLine,1,numOfLines),repmat(conductivityOfGroundLine,1,numOfGroundLines)];%������е���ĵ絼�ʣ�������
%%�����ŵ��ʼ���
Epsilon0 =8.854187817e-12;% 8.854187817e-12;%��ս�糣������λF/m
miu0=4*pi*1e-7;%��մŵ��ʣ���λH/m
miu=miur*miu0;%�����ŵ���
%%��Ƶ����͸����ȼ���
omega=2*pi*freq;%�����
s=omega*1i;%s==j��
p=1./sqrt(s*miu*conductivityOfGround);%��͸�����,�õ��������絼��
%% Y��������
P=zeros(numOfConductors);%��λϵ������
for i=1:numOfConductors
    for j=1:numOfConductors
        Xdifference=(XcoordinateOfConductor(i)-XcoordinateOfConductor(j));
        Ydifference=(YcoordinateOfConductor(i)-YcoordinateOfConductor(j));
        Ysum=(YcoordinateOfConductor(i)+YcoordinateOfConductor(j));
        if i==j
            P(i,j)=1/(2*pi*Epsilon0)*log(Ysum/radiusOfConductor(i));
        else
            P(i,j)=1/(2*pi*Epsilon0)*0.5*log((Xdifference^2+Ysum^2)/(Xdifference^2+Ydifference^2));
        end
    end
end
B=inv(P);%�����Ӧϵ������
Y=arrayfun(@(x) s(x)*B ,1:length(s),'UniformOutput',false);
Y=cat(3,Y{:});
%%��Y����G
% Y=Y+1e-10;
%% Z��������
%%carson��ʽ���裺���Ѹ�Ϊ������
%%��1����������λ�Ƶ���
%%��2������Ϊ�����޴����ý��
%%��3���������໥ƽ�к����޳���
Z=zeros(numOfConductors,numOfConductors,length(s));%�迹����,��ά

for i=1:numOfConductors
    for j=1:numOfConductors
        Xdifference=(XcoordinateOfConductor(i)-XcoordinateOfConductor(j));
        Ydifference=(YcoordinateOfConductor(i)-YcoordinateOfConductor(j));
        Ysum=(YcoordinateOfConductor(i)+YcoordinateOfConductor(j));
        if i==j
            %%���迹
%             %%���ֺ�������,������ȥ��
%             f33=@(u,y,k) exp(-2*y*u)./(u+sqrt(u.^2+1/(p(k)^2)));%��������
%             f33(u,YcoordinateOfConductor(i),k)
%             f3=arrayfun(@(k)s(k)*miu0/pi*integral(@(u) f33(u,YcoordinateOfConductor(i),k),0,inf),1:length(s),'UniformOutput',false);%��ʽ�����Ԫ������
%             selfImpedance3=cat(3,f3{:});%��Ԫ������ת��Ϊ��ά����
%             figure(i);
%             plot(squeeze(selfImpedance3));
%             if i==2
%                 E1=selfImpedance3;
%             end
            %%������������
            f22=@(conductivity,radius,k) sqrt(s(k)*miu/conductivity)/(2*pi*radius)*besseli(0,radius*sqrt(s(k)*miu*conductivity))/besseli(1,radius*sqrt(s(k)*miu*conductivity));
            f2=arrayfun(@(k) f22(conductivityOfConductor(i),radiusOfConductor(i),k),1:length(s),'UniformOutput',false);%������������·�ĵĵ絼����뾶
            selfImpedance2(i,j,:)=cat(3,f2{:});%��Ԫ������ת��Ϊ��ά����
            %%������
            f1=arrayfun(@(k) s(k)*miu0/2/pi*log(2*(YcoordinateOfConductor(i)+p(k))/radiusOfConductor(i)),1:length(s),'UniformOutput',false);
            selfImpedance1(i,j,:)=cat(3,f1{:});%��Ԫ������ת��Ϊ��ά����
            %%�������
            Z(i,j,:)=selfImpedance1(i,j,:)+selfImpedance2(i,j,:);%��ά����
            %��selfImpedance1��selfImpedance2������ϣ�i��j����������Ϊȱʡ��Ϊ�˲鿴�м����
        else
            %%���迹
%             %%������,�����迹�ĵĵ������޸Ķ���
%             f33=@(u,ysum,xdifference,k) cos(u*xdifference).*exp(-ysum*u)./(u+sqrt(u.^2+1/(p(k)^2)));%��������
%             f3=arrayfun(@(k)s(k)*miu0/pi*integral(@(u) f33(u,Ysum,Xdifference,k),0,inf),1:length(s),'UniformOutput',false);
%             mutualImpedance3=cat(3,f3{:});%��Ԫ������ת��Ϊ��ά����
            %%����������迹�ĵ�һ���޸Ķ���
            f1=arrayfun(@(k) s(k)*miu0/2/pi*log(sqrt((Xdifference^2+(Ysum+2*p(k))^2)/(Xdifference^2+Ydifference^2))),1:length(s),'UniformOutput',false);
            mutualImpedance1=cat(3,f1{:});%��Ԫ������ת��Ϊ��ά����
            %%�������
            Z(i,j,:)=mutualImpedance1;%��ά���� 
            
        end
    end
end
%% ������������ѵ��ߺϲ�
%%��ȡ��ʦ�İ취����������Ӧ��ͬ��
%%����ϲ�����
[P1, P2]=deal(zeros(numOfConductors));%P1��˾���P2�ҳ˾���numOfphase��������
PP=diag(ones(1,N));%����N����N���ѵ��߻����þ���
PP(1,:)=ones(1,N);
for n=1:numOfPhase
    index1=(n-1)*N;
    index2=n*N;
    P1(index1+1:index2,index1+1:index2)=PP;
    P2(index1+1:index2,index1+1:index2)=PP.';
end
%%Y�ϲ��������ǵ���
Y=arrayfun(@(x) P1*Y(:,:,x)*P2 ,1:length(s),'UniformOutput',false);%�õ�����Ԫ������
Y=cat(3,Y{:});%ת��Ϊ��ά����
%%�����ľ���ȡ�ϲ����Ԫ�أ�
Yphase=zeros(numOfPhase,numOfPhase,length(s));
for i=1:numOfPhase
    for j=1:numOfPhase
        Yphase(i,j,:)=Y((i-1)*N+1,(j-1)*N+1,:);
    end
end

%%Z�ϲ��������ǵ���
inverseZ=arrayfun(@(x) P1*Z(:,:,x)^-1*P2 ,1:length(s),'UniformOutput',false);%�õ�����Ԫ������
inverseZ=cat(3,inverseZ{:});%ת��Ϊ��ά����
%%%�����ľ���ȡ�ϲ����Ԫ�أ�
simplifyInverseZ=zeros(numOfPhase,numOfPhase,length(s));
for i=1:numOfPhase
    for j=1:numOfPhase
        simplifyInverseZ(i,j,:)=inverseZ((i-1)*N+1,(j-1)*N+1,:);
    end
end
Zphase=arrayfun(@(x) simplifyInverseZ(:,:,x)^-1 ,1:length(s),'UniformOutput',false);
Zphase=cat(3,Zphase{:});%ת��Ϊ��ά����
%% �����
%�����ǻ�λ
if isTransposition
    %���ǻ�λ
    Zs=1/3*(Zphase(1,1,:)+Zphase(2,2,:)+Zphase(3,3,:));%�Խ���Ԫ��
    Zm=1/3*(Zphase(1,2,:)+Zphase(2,3,:)+Zphase(3,1,:));%�ǶԽ���Ԫ��
    Ys=1/3*(Yphase(1,1,:)+Yphase(2,2,:)+Yphase(3,3,:));
    Ym=1/3*(Yphase(1,2,:)+Yphase(2,3,:)+Yphase(3,1,:));
    Zpositive=Zs-Zm;
    Zzero=Zs+2*Zm;
    Ypositive=Ys-Ym;
    Yzero=Ys+2*Ym;
    for i=1:length(s)
    Zphase(:,:,i)=diag([1 1 1]).*Zs(:,:,i)+(ones(numOfPhase)-diag([1,1,1])).*Zm(:,:,i);
    Yphase(:,:,i)=diag([1,1,1]).*Ys(:,:,i)+(ones(numOfPhase)-diag([1,1,1])).*Ym(:,:,i);
    end
else
    %δ��λ
    a=exp(2j*pi/3);
    A=[1 1 1;1 a^2 a;1 a a^2];
    for i=1:length(s)        
        Zsym(:,:,i)=A\Zphase(:,:,i)*A;%�ԳƷ��������迹����
        Ysym(:,:,i)=A\Yphase(:,:,i)*A;%�ԳƷ����������ɾ���
    end
end
%% ���
Z=Zphase;
Y=Yphase;
end
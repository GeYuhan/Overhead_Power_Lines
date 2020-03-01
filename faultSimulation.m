function [t,VA,VB,VC,IA,IB,IC]=faultSimulation(HFit,YcFit,tau)
%%ʱ�����
%ʱ���������
%           Ny    G_i                 Ng              Nh(k)  R_k,i
%   Yc==G_0+sum --------          H==sum exp��-s��_k�� sum  ---------
%           i=1   s-q_i              k=1               i=1  s-p_k,i
%
%Yc״̬�ռ�������ò�����
%
%     (2+��t*q_i)              (��t*G_i)                 
%a_i==-----------       G_i~==-----------   G_i^==(a_i+1)*G_i~
%     (2-��t*q_i)             (2-��t*q_i)
%
%         Ny
%   G==G_0+sumG_i^   
%         i=1
%
%H״̬�ռ�������ò�����
%
%        (2+��t*p_k,i)                 (��t*R_k,i)
%a_k,i==---------------      R_k,i~==--------------
%        (2-��t*p_k,i)                (2-��t*p_k,i)
% 
%������parametersOfLine
% 
% 
%% ��Ļ��ʾ
disp(['���ڼ���']);
%% ����ʱ������
% clear Isr
% clc
T= 0.525;% Observation time
deltaT=50e-6;%��λ��,ͨ�������ڵ�1/20
Niteration=fix(T/deltaT);%��������
% tauInteger= fix(tau./deltaT);    % �õ���tau�Ԧ�t�ı�����tau��һ�����������洢��������ʱ
tauMax = fix(max(tau)./deltaT);    %���Ħ�
[Nc,~]=size(HFit.C(:,:,1));
%% ���ɴ�������
line1(1)=transmissionLine(T,deltaT,tau,HFit,YcFit);
line1(2)=transmissionLine(T,deltaT,tau,HFit,YcFit);
%% A��B���˵�ѹ����
%��ǰʱ�̵�ѹ����ֵ
vn=zeros(2*Nc,1);%vn=[vAn;vBn]
vAn=zeros(Nc,1);
vBn=zeros(Nc,1);
in=zeros(2*Nc,1);%in=[iAn;iBn]
iAn=zeros(Nc,1);
iBn=zeros(Nc,1);
%�����ѹ��������
VA  = zeros(Nc,Niteration+tauMax);       % V0ltage at node A 
VB  = zeros(Nc,Niteration+tauMax);       % V0ltage at node B 
IA  = zeros(Nc,Niteration+tauMax);       %���� at node A 
IB  = zeros(Nc,Niteration+tauMax);       % ����at node B 
%% A��
t= (0:deltaT:(Niteration-1)*deltaT);
Isr(3,:) = sin(100*pi*t);
Isr(1,:) = sin(100*pi*t-2*pi/3);
Isr(2,:) = sin(100*pi*t-4*pi/3);
Isr=1e3*Isr;
Ys  = diag([1; 1; 1]); 
nodeA=node(3);%3����ڵ�
%% B�� Fault��
tf=0.425;%s.���Ͽ�ʼʱ��
td=0.05;%s.���ϳ���ʱ��
kf=fix(tf/deltaT);
kclose=fix((tf+td)/deltaT);
%%�����е�k
% tag=(k>=kf&k<kclose)*1e5;
YfCanculation=@(k) diag([1/1e6; 1/1e6; 1/1e6])*(k>=kf&k<kclose)*1e8;%��ʱ��仯�ĵ�����ģ�����
If=zeros(3,1);
nodeB=node(3);
%% C�� ��������
Il=1e3*Isr;%��ֻ��Ĭ����ͬ������
Yl=diag([1e3; 1e3; 1e3]);     % Admittance of load connected at node B ��
nodeC=node(3);
%% ��������
% for k = tauMax+2:Niteration+tauMax-3  %����Ϊʲô�������,���ɲ���
for k=1:Niteration
    line1(1)=line1(1).iteration;%�ܷ�һ���԰������߶�����?���ǲ���
    line1(2)=line1(2).iteration;
    GGTilde{1}=line1(1).getGGTilde;
    GGTilde{2}=line1(2).getGGTilde;
    iHist0{1}=line1(1).getiHist0;
    iHist0{2}=line1(2).getiHist0;
    %%ȡ��ÿһ�ڵ��ϵĵ�Ч�������Ч����Դ
    %%��Ч����
    GA=GGTilde{1}(1:end/2,1:end/2);
    GB(:,:,1)=GGTilde{1}(end/2+1:end,end/2+1:end);
    GB(:,:,2)=GGTilde{2}(1:end/2,1:end/2);
    GC=GGTilde{2}(end/2+1:end,end/2+1:end);
    %%��Ч����Դ
    iHist0A=iHist0{1}(1:end/2);
    iHist0B(:,:,1)=iHist0{1}(end/2+1:end);
    iHist0B(:,:,2)=iHist0{2}(1:end/2);
    iHist0C=iHist0{2}(end/2+1:end);
    %% ��ѹ��������
    Yf=YfCanculation(k);%��������ϵ�ӵص���
    [nodeA,vAn,iAn]=nodeA.caculation(Ys,Isr(:,k),GA,iHist0A);
    [nodeB,vBn,iBn]=nodeB.caculation(Yf,If,GB,iHist0B);
    [nodeC,vCn,iCn]=nodeC.caculation(Yl,Il(:,k),GC,iHist0C);
    vn(:,1)=[vAn;vBn];
    vn(:,2)=[vBn;vCn];
    in(:,1)=[iAn;iBn(:,1)];
    in(:,2)=[iBn(:,2);iCn];
    %% ����
    line1(1)=line1(1).setvTemp(vn(:,1));
    line1(1)=line1(1).setiTemp(in(:,1));
    line1(2)=line1(2).setvTemp(vn(:,2));
    line1(2)=line1(2).setiTemp(in(:,2));
end
%%��¼��ʱ�̵�ѹ����ֵ
[VA,IA]=nodeA.getUI;%nodeA
[VB,IB]=nodeB.getUI;%nodeB
[VC,IC]=nodeC.getUI;%nodeC
%% ��ѹ��ͼ
draw=0;
if draw==1
    figure
    plot(t,VA,'--',t,VC)
    ylabel('��ֵ V')
    xlabel('ʱ�� s')
    title('ʼ�˵�ѹ���ն˵�ѹ')
    legend('ʼ��A��' , 'ʼ��B��' , 'ʼ��C��' , '�ն�A��' , '�ն�B��' , '�ն�C��')
    %% ������ͼ
%     figure
%     plot(t,IA(2,:),'--')
%     figure
%     plot(t,IC(2,:))
    figure
    plot(t,IA,'--',t,IC)
    ylabel('��ֵ A')
    xlabel('ʱ�� s')
    title('ʼ�˵������ն˵���')
    legend('ʼ��A��' , 'ʼ��B��' , 'ʼ��C��' , '�ն�A��' , '�ն�B��' , '�ն�C��')
end%if draw==1
%% ��Ļ��ʾ
disp(['�������']);

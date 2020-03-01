function [t,VA,VB,IA,IB]=doSimulation(HFit,YcFit,tau)
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
%% ����ʱ������
% Isr=[];
% clc
T= 0.5e-3;% Observation time
deltaT=0.1e-6;%��λ��,ͨ�������ڵ�1/20
Niteration=fix(T/deltaT);%��������
% tauInteger= fix(tau./deltaT);    % �õ���tau�Ԧ�t�ı�����tau��һ�����������洢��������ʱ
tauMax = fix(max(tau)./deltaT);    %���Ħ�
[Nc,~]=size(HFit.C(:,:,1));
%% ���ɴ�������
line1=transmissionLine(T,deltaT,tau,HFit,YcFit);
%% ��Դ��������
% t= (0:deltaT:(Niteration+tauMax-1)*deltaT);     % һ��Niteration+tauMax����
t= (0:deltaT:(Niteration-1)*deltaT);     % һ��Niteration+tauMax����
% Ks = menu('CHOOSE THE TYPE OF INPUT SOURCE' , '1 -unit step' , '2 -sinusoidal'); 
Ks=1;
if Ks == 1  % unit step source   %�����Ϊ�׵粨
    Isr(1,:) = 50*(exp(-20000*t)-exp(-1666666.6*t));
%     Isr(1,:) = 0;%A�಻��Դ
%     Isr(2,:) = 50*(exp(-20000*t)-exp(-1666666.6*t)); 
    Isr(2,:) = 0;%B�಻��Դ
%     Isr(3,:) = 50*(exp(-20000*t)-exp(-1666666.6*t));
    Isr(3,:) = 0;%C�಻��Դ
else                                                                                                                                                                                                                                                                                                                                                                           
    Ks == 2  % sinusoidal source
    Isr(1,:) = sin(100*pi*t);
    Isr(2,:) = sin(100*pi*t-2*pi/3);
    Isr(3,:) = sin(100*pi*t-4*pi/3);
end
%ѭ�������ɵ�Դ�����
%% A��B���˵�ѹ����
%��ǰʱ�̵�ѹ����ֵ
vn=zeros(2*Nc,1);%vn=[vAn;vBn]
vAn=zeros(Nc,1);
vBn=zeros(Nc,1);
in=zeros(2*Nc,1);%in=[iAn;iBn]
iAn=zeros(Nc,1);
iBn=zeros(Nc,1);
%�����ѹ��������
% VA  = zeros(Nc,Niteration+tauMax);       % V0ltage at node A 
% VB  = zeros(Nc,Niteration+tauMax);       % V0ltage at node B 
% IA  = zeros(Nc,Niteration+tauMax);       %���� at node A 
% IB  = zeros(Nc,Niteration+tauMax);       % ����at node B 
VA  = zeros(Nc,Niteration);       % V0ltage at node A 
VB  = zeros(Nc,Niteration);       % V0ltage at node B 
IA  = zeros(Nc,Niteration);       %���� at node A 
IB  = zeros(Nc,Niteration);       % ����at node B 
%% ��ѹ���������
%Note����Ҫ��ϳ��¾���
Ys  = diag([1/600; 1/600; 1/600]);       % Admittance of the source, connected at node A ����Դ�ڵ���
% Ys  = diag([1/600; 1e3; 1e3]);       % Admittance of the source, connected at node A ����Դ�ڵ���
% Ys  = diag([1e3; 1e30; 1e30]);       % Admittance of the source, connected at node A ����Դ�ڵ���
Yl  =diag([1/600; 1/600; 1/600]);     % Admittance of load connected at node B ��
% Yl  =diag([1/600; 1e3; 1e3]);     % Admittance of load connected at node B ��
Ysl=diag([diag(Ys);diag(Yl)]);%
%% ��������
% for k = tauMax+2:Niteration+tauMax-3  %����Ϊʲô�������,���ɲ���
% for k=1:Niteration+tauMax-3
for k=1:Niteration
    if k==3872
        nothing=1;
    end
    line1=line1.iteration;   
    GGTilde=line1.getGGTilde;
    iHist0=line1.getiHist0;
    %% ��ѹ��������
    %��ѹ����
    Gsum=GGTilde+Ysl;
    vn=Gsum\(iHist0+[Isr(:,k);zeros(3,1)]);
    line1=line1.setvTemp(vn);
    vAn=vn(1:end/2);
    vBn=vn(end/2+1:end);
    %��������
    in=GGTilde*vn-iHist0;  
    line1=line1.setiTemp(in);
    iAn=in(1:end/2);
    iBn=in(end/2+1:end);
    %%��¼��ʱ�̵�ѹ����ֵ
    VA(:,k) = vAn;%nodeA
    VB(:,k) = vBn;%nodeB
    IA(:,k) = iAn;%nodeA
    IB(:,k) = iBn;%nodeB
end
%% ��ѹ��ͼ
draw=0;
if draw==1
% % figure
% % plot(t,VA(1,:))
% % hold on
% % plot(t,VA(2,:))
% % plot(t,VA(3,:))
% % plot(t,VB(1,:))
% % plot(t,VB(2,:))
% % plot(t,VB(3,:))
% % hold off
% figure
% plot(t,VA) 
% ylabel('Amplitude in volts')
% xlabel('Time in seconds') 
% title('�Ͷ˵�ѹ')
% legend('Sending end phase A' , 'Sending end phase B' , 'Sending end phase C')
% figure
% plot(t,VB) 
% ylabel('Amplitude in volts')
% xlabel('Time in seconds')
% title('�ܶ˵�ѹ')
% legend('Receiving end phase A' , 'Receiving end phase B' , 'Receiving end phase C')
figure
plot(t,VA,'--',t,VB) 
ylabel('��ֵ V')
xlabel('ʱ�� s')
title('ʼ�˵�ѹ���ն˵�ѹ')
legend('ʼ��A��' , 'ʼ��B��' , 'ʼ��C��' , '�ն�A��' , '�ն�B��' , '�ն�C��')   
%% ������ͼ
figure
plot(t,IA,'--',t,IB) 
ylabel('��ֵ A')
xlabel('ʱ�� s')
title('ʼ�˵������ն˵���')
legend('ʼ��A��' , 'ʼ��B��' , 'ʼ��C��' , '�ն�A��' , '�ն�B��' , '�ն�C��')  
end%if draw==1

end%function
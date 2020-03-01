function [t,VA,VB,IA,IB]=doSimulation(HFit,YcFit,tau)
%%时域计算
%时域参数计算
%           Ny    G_i                 Ng              Nh(k)  R_k,i
%   Yc==G_0+sum --------          H==sum exp（-sτ_k） sum  ---------
%           i=1   s-q_i              k=1               i=1  s-p_k,i
%
%Yc状态空间迭代所用参数：
%
%     (2+Δt*q_i)              (Δt*G_i)                 
%a_i==-----------       G_i~==-----------   G_i^==(a_i+1)*G_i~
%     (2-Δt*q_i)             (2-Δt*q_i)
%
%         Ny
%   G==G_0+sumG_i^   
%         i=1
%
%H状态空间迭代所用参数：
%
%        (2+Δt*p_k,i)                 (Δt*R_k,i)
%a_k,i==---------------      R_k,i~==--------------
%        (2-Δt*p_k,i)                (2-Δt*p_k,i)
% 
%先运行parametersOfLine
% 
% 
%% 仿真时间设置
% Isr=[];
% clc
T= 0.5e-3;% Observation time
deltaT=0.1e-6;%单位秒,通常是周期的1/20
Niteration=fix(T/deltaT);%迭代次数
% tauInteger= fix(tau./deltaT);    % 得到了tau对Δt的倍数。tau是一个行向量，存储着三个延时
tauMax = fix(max(tau)./deltaT);    %最大的τ
[Nc,~]=size(HFit.C(:,:,1));
%% 生成传输线类
line1=transmissionLine(T,deltaT,tau,HFit,YcFit);
%% 电源数据生成
% t= (0:deltaT:(Niteration+tauMax-1)*deltaT);     % 一共Niteration+tauMax个点
t= (0:deltaT:(Niteration-1)*deltaT);     % 一共Niteration+tauMax个点
% Ks = menu('CHOOSE THE TYPE OF INPUT SOURCE' , '1 -unit step' , '2 -sinusoidal'); 
Ks=1;
if Ks == 1  % unit step source   %后面改为雷电波
    Isr(1,:) = 50*(exp(-20000*t)-exp(-1666666.6*t));
%     Isr(1,:) = 0;%A相不加源
%     Isr(2,:) = 50*(exp(-20000*t)-exp(-1666666.6*t)); 
    Isr(2,:) = 0;%B相不加源
%     Isr(3,:) = 50*(exp(-20000*t)-exp(-1666666.6*t));
    Isr(3,:) = 0;%C相不加源
else                                                                                                                                                                                                                                                                                                                                                                           
    Ks == 2  % sinusoidal source
    Isr(1,:) = sin(100*pi*t);
    Isr(2,:) = sin(100*pi*t-2*pi/3);
    Isr(3,:) = sin(100*pi*t-4*pi/3);
end
%循环内生成电源矩阵把
%% A，B两端电压电流
%当前时刻电压电流值
vn=zeros(2*Nc,1);%vn=[vAn;vBn]
vAn=zeros(Nc,1);
vBn=zeros(Nc,1);
in=zeros(2*Nc,1);%in=[iAn;iBn]
iAn=zeros(Nc,1);
iBn=zeros(Nc,1);
%输出电压电流向量
% VA  = zeros(Nc,Niteration+tauMax);       % V0ltage at node A 
% VB  = zeros(Nc,Niteration+tauMax);       % V0ltage at node B 
% IA  = zeros(Nc,Niteration+tauMax);       %电流 at node A 
% IB  = zeros(Nc,Niteration+tauMax);       % 电流at node B 
VA  = zeros(Nc,Niteration);       % V0ltage at node A 
VB  = zeros(Nc,Niteration);       % V0ltage at node B 
IA  = zeros(Nc,Niteration);       %电流 at node A 
IB  = zeros(Nc,Niteration);       % 电流at node B 
%% 电压计算参数，
%Note：需要组合乘新矩阵
Ys  = diag([1/600; 1/600; 1/600]);       % Admittance of the source, connected at node A ，电源内导纳
% Ys  = diag([1/600; 1e3; 1e3]);       % Admittance of the source, connected at node A ，电源内导纳
% Ys  = diag([1e3; 1e30; 1e30]);       % Admittance of the source, connected at node A ，电源内导纳
Yl  =diag([1/600; 1/600; 1/600]);     % Admittance of load connected at node B ，
% Yl  =diag([1/600; 1e3; 1e3]);     % Admittance of load connected at node B ，
Ysl=diag([diag(Ys);diag(Yl)]);%
%% 迭代计算
% for k = tauMax+2:Niteration+tauMax-3  %不懂为什么是这个数,依旧不懂
% for k=1:Niteration+tauMax-3
for k=1:Niteration
    if k==3872
        nothing=1;
    end
    line1=line1.iteration;   
    GGTilde=line1.getGGTilde;
    iHist0=line1.getiHist0;
    %% 电压电流计算
    %电压计算
    Gsum=GGTilde+Ysl;
    vn=Gsum\(iHist0+[Isr(:,k);zeros(3,1)]);
    line1=line1.setvTemp(vn);
    vAn=vn(1:end/2);
    vBn=vn(end/2+1:end);
    %电流计算
    in=GGTilde*vn-iHist0;  
    line1=line1.setiTemp(in);
    iAn=in(1:end/2);
    iBn=in(end/2+1:end);
    %%记录此时刻电压电流值
    VA(:,k) = vAn;%nodeA
    VB(:,k) = vBn;%nodeB
    IA(:,k) = iAn;%nodeA
    IB(:,k) = iBn;%nodeB
end
%% 电压绘图
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
% title('送端电压')
% legend('Sending end phase A' , 'Sending end phase B' , 'Sending end phase C')
% figure
% plot(t,VB) 
% ylabel('Amplitude in volts')
% xlabel('Time in seconds')
% title('受端电压')
% legend('Receiving end phase A' , 'Receiving end phase B' , 'Receiving end phase C')
figure
plot(t,VA,'--',t,VB) 
ylabel('幅值 V')
xlabel('时间 s')
title('始端电压与终端电压')
legend('始端A相' , '始端B相' , '始端C相' , '终端A相' , '终端B相' , '终端C相')   
%% 电流绘图
figure
plot(t,IA,'--',t,IB) 
ylabel('幅值 A')
xlabel('时间 s')
title('始端电流与终端电流')
legend('始端A相' , '始端B相' , '始端C相' , '终端A相' , '终端B相' , '终端C相')  
end%if draw==1

end%function
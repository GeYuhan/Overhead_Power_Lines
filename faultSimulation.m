function [t,VA,VB,VC,IA,IB,IC]=faultSimulation(HFit,YcFit,tau)
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
%% 屏幕提示
disp(['正在计算']);
%% 仿真时间设置
% clear Isr
% clc
T= 0.525;% Observation time
deltaT=50e-6;%单位秒,通常是周期的1/20
Niteration=fix(T/deltaT);%迭代次数
% tauInteger= fix(tau./deltaT);    % 得到了tau对Δt的倍数。tau是一个行向量，存储着三个延时
tauMax = fix(max(tau)./deltaT);    %最大的τ
[Nc,~]=size(HFit.C(:,:,1));
%% 生成传输线类
line1(1)=transmissionLine(T,deltaT,tau,HFit,YcFit);
line1(2)=transmissionLine(T,deltaT,tau,HFit,YcFit);
%% A，B两端电压电流
%当前时刻电压电流值
vn=zeros(2*Nc,1);%vn=[vAn;vBn]
vAn=zeros(Nc,1);
vBn=zeros(Nc,1);
in=zeros(2*Nc,1);%in=[iAn;iBn]
iAn=zeros(Nc,1);
iBn=zeros(Nc,1);
%输出电压电流向量
VA  = zeros(Nc,Niteration+tauMax);       % V0ltage at node A 
VB  = zeros(Nc,Niteration+tauMax);       % V0ltage at node B 
IA  = zeros(Nc,Niteration+tauMax);       %电流 at node A 
IB  = zeros(Nc,Niteration+tauMax);       % 电流at node B 
%% A端
t= (0:deltaT:(Niteration-1)*deltaT);
Isr(3,:) = sin(100*pi*t);
Isr(1,:) = sin(100*pi*t-2*pi/3);
Isr(2,:) = sin(100*pi*t-4*pi/3);
Isr=1e3*Isr;
Ys  = diag([1; 1; 1]); 
nodeA=node(3);%3导体节点
%% B端 Fault点
tf=0.425;%s.故障开始时间
td=0.05;%s.故障持续时间
kf=fix(tf/deltaT);
kclose=fix((tf+td)/deltaT);
%%迭代中的k
% tag=(k>=kf&k<kclose)*1e5;
YfCanculation=@(k) diag([1/1e6; 1/1e6; 1/1e6])*(k>=kf&k<kclose)*1e8;%随时间变化的导纳来模拟故障
If=zeros(3,1);
nodeB=node(3);
%% C端 无穷大电网
Il=1e3*Isr;%我只能默认是同步的了
Yl=diag([1e3; 1e3; 1e3]);     % Admittance of load connected at node B ，
nodeC=node(3);
%% 迭代计算
% for k = tauMax+2:Niteration+tauMax-3  %不懂为什么是这个数,依旧不懂
for k=1:Niteration
    line1(1)=line1(1).iteration;%能否一次性把两条线都算了?答案是不能
    line1(2)=line1(2).iteration;
    GGTilde{1}=line1(1).getGGTilde;
    GGTilde{2}=line1(2).getGGTilde;
    iHist0{1}=line1(1).getiHist0;
    iHist0{2}=line1(2).getiHist0;
    %%取出每一节点上的等效导纳与等效电流源
    %%等效导纳
    GA=GGTilde{1}(1:end/2,1:end/2);
    GB(:,:,1)=GGTilde{1}(end/2+1:end,end/2+1:end);
    GB(:,:,2)=GGTilde{2}(1:end/2,1:end/2);
    GC=GGTilde{2}(end/2+1:end,end/2+1:end);
    %%等效电流源
    iHist0A=iHist0{1}(1:end/2);
    iHist0B(:,:,1)=iHist0{1}(end/2+1:end);
    iHist0B(:,:,2)=iHist0{2}(1:end/2);
    iHist0C=iHist0{2}(end/2+1:end);
    %% 电压电流计算
    Yf=YfCanculation(k);%计算出故障点接地导纳
    [nodeA,vAn,iAn]=nodeA.caculation(Ys,Isr(:,k),GA,iHist0A);
    [nodeB,vBn,iBn]=nodeB.caculation(Yf,If,GB,iHist0B);
    [nodeC,vCn,iCn]=nodeC.caculation(Yl,Il(:,k),GC,iHist0C);
    vn(:,1)=[vAn;vBn];
    vn(:,2)=[vBn;vCn];
    in(:,1)=[iAn;iBn(:,1)];
    in(:,2)=[iBn(:,2);iCn];
    %% 更新
    line1(1)=line1(1).setvTemp(vn(:,1));
    line1(1)=line1(1).setiTemp(in(:,1));
    line1(2)=line1(2).setvTemp(vn(:,2));
    line1(2)=line1(2).setiTemp(in(:,2));
end
%%记录此时刻电压电流值
[VA,IA]=nodeA.getUI;%nodeA
[VB,IB]=nodeB.getUI;%nodeB
[VC,IC]=nodeC.getUI;%nodeC
%% 电压绘图
draw=0;
if draw==1
    figure
    plot(t,VA,'--',t,VC)
    ylabel('幅值 V')
    xlabel('时间 s')
    title('始端电压与终端电压')
    legend('始端A相' , '始端B相' , '始端C相' , '终端A相' , '终端B相' , '终端C相')
    %% 电流绘图
%     figure
%     plot(t,IA(2,:),'--')
%     figure
%     plot(t,IC(2,:))
    figure
    plot(t,IA,'--',t,IC)
    ylabel('幅值 A')
    xlabel('时间 s')
    title('始端电流与终端电流')
    legend('始端A相' , '始端B相' , '始端C相' , '终端A相' , '终端B相' , '终端C相')
end%if draw==1
%% 屏幕提示
disp(['计算完毕']);

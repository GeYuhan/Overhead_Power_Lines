%%%%%%%%%%%%%%%%%%%%电缆参数文件%%%%%%%%%%%%%
Case=3;%%%电缆层数%Case==1,只含导芯，Case==2,含导芯及屏蔽层，Case==3,含导芯、屏蔽层及护层
Phase=1;%%%电缆相数%Phase==1,表示单相，Phase==3,三相电缆
Erec=1;%Erec=1,表示埋地，E_rec=2,表示架空
SOL=1;%%土壤层数%SOL==1，表示单层土壤，SOL==2表示双层。

mu0=4e-7*pi;%空气磁导率，H/m
espo0=8.85e-12;%真空介电常数，F/m
gama0=0;%%%空气的电导率，S/m
%%%%%不变的%%%%%%%%%%%

%%%%%导芯参数%%%%%%%%
r1=0.6/100;%m，导芯内径
r2=1.635/100;%m，导芯外径
rou1=1.82e-8;%Ω.m，导芯电阻率；
mu1=mu0;
%%%%屏蔽层铅芯参数%%%%
r3=4.835/100;%m，铅芯内径
r4=5.235/100;%m，铅芯外径
rou2=2.8e-7;%Ω.m，铅芯电阻率；
mu2=mu0;
%%%%金属护层黄铜带层参数
r5=5.275/100;%m，黄铜带内径
r6=5.355/100;%m，黄铜带外径
rou3=7.54e-8;%Ω.m，黄铜带电阻率；
mu3=mu0;
%%%%%电缆绝缘材料相关参数%%%%%
espo1=3.5*espo0;%%%绝缘层的介电常数
espo2=3.5*espo0;%%%粘接带的介电常数
espo3=8*espo0;%%%绝缘层的介电常数
%%%%护套层外径，大地内径%%%%%
r7=5.905/100;%cm，护套层外径
%%%%大地相关参数%%%%%%%%%
%%%第一层土壤相关参数%%%
rou4=100;%Ω.m，土壤电阻率；
mu4=mu0;
espo4=espo0;%%%大地的介电常数
d1=20;%m第一层土壤的厚度
%%%第二层土壤相关参数
rou5=50;%Ω.m，土壤电阻率；
mu5=mu0;
espo5=espo0;%%%大地的介电常数

% h=0.6;%m,埋设深度%%%%
%%%%%%三相导体的情况%%%%%
ha=0.6;%m，A相导体
hb=0.6;%m,B相导体
hc=0.6;%m,C相导体
xa=5;%m,A-B之间的距离
xb=10;%m,A-C之间的距离
xc=15;%m,B-C之间的距离

habc=[ha, hb, hc];%埋地深度矩阵
xabc=[xa,xb,xc];%相间距





%%%%画出电缆结构
% alfa=0:2*pi/100:2*pi;
% figure
% plot(r1*cos(alfa),-1*h+r1*sin(alfa),'r-');
% hold on
% plot(r2*cos(alfa),-1*h+r2*sin(alfa),'k-');
% hold on
% plot(r3*cos(alfa),-1*h+r3*sin(alfa),'b-');
% hold on
% plot(r4*cos(alfa),-1*h+r4*sin(alfa),'k-');
% hold on
% plot(r5*cos(alfa),-1*h+r5*sin(alfa),'m-');
% hold on
% plot(r6*cos(alfa),-1*h+r6*sin(alfa),'k-');




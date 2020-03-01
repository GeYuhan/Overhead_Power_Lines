%% 频率设置
lowerLimit=0.5;%HZ
highLimit=1e6;%Hz
totalSolutionIncreaments=20;%采样间隔数目
global freq
% freq=50;
freq=logspace(log10(lowerLimit),log10(highLimit),totalSolutionIncreaments+1);%对数采样，增加低频精度
%% 线路几何与物理参数
transmissionLineType=1;%电力传输线种类，0表示电力电缆，1表示架空输电线
l=100e3;%线路长度
if transmissionLineType
    [Z,Y]=overheadLineParameterCaculation;
else
    [Z,Y]=cableParameterCaculation;
end
%% 矢量拟合相关参数设置
%%H拟合设置参数
global NpMaxH MRPRtolerance
NpMaxH=19;%每个模量最大极点数目
toleranceH=0.2e-2;%（0.2%参考自pscad）
MRPRtolerance=100;%最大留数极点比阈值
%%Yc拟合参数设置
NpMaxYc=50;%最大极点数目
toleranceYc=0.2e-2;
%%绘图
paint=1;%值为1时绘制，0不输出图像
%% 执行拟合脚本 
curveFitting
%% 时域计算
tic
%雷电波
[t,VA,VB,IA,IB]=doSimulation(HFit,YcFit,tau);
% [t,VA,VB,IA,IB]=doSimulation(HFitPscad,YcFitPscad,tauPscad);
%单相接地故障
% [t,VA,VB,VC,IA,IB,IC]=faultSimulation(HFit,YcFit,tau);
toc
%%simulationSet.m 是时域设置文件
%%
figure
plot(t,VA(3,:))
hold on
plot(t,VB(3,:))
hold off
%% 故障误差计算并画图
for i=1:length(t)
    if t(i)>=0.4
        tt=t(i:end);
        IIA=IA(:,i:end);
        IIC=IC(:,i:end);
        break
    end
end
for i=1:3
    figure(i)
%     subplot(2,1,1)
%     plot(Domain,1e3*ISend2(:,i),'--',Domain1,1e3*ISend3(:,i));
    plot(tt,IIA(i,:),'-.',Domain,1e3*ISend3(:,i));
    legend('自编程序','Pscad')
    xlim([0.4,0.52])
%     subplot(2,1,2)
%     Error(i,:)=abs(IIA(i,:).'-1e3*ISend2(:,i))./abs(1e3*ISend2(:,i));
%     plot(tt,Error(i,:),'r')
end
%% 
for i=1:length(tt)
    [value1, index1]=max(abs(IIA.'));
end
for i=1:length(Domain)
    [value2, index2]=max(abs(1e3*ISend3));
end
abs((value1-value2)./max(value2,value1))
    
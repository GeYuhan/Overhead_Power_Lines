%% Ƶ������
lowerLimit=0.5;%HZ
highLimit=1e6;%Hz
totalSolutionIncreaments=20;%���������Ŀ
global freq
% freq=50;
freq=logspace(log10(lowerLimit),log10(highLimit),totalSolutionIncreaments+1);%�������������ӵ�Ƶ����
%% ��·�������������
transmissionLineType=1;%�������������࣬0��ʾ�������£�1��ʾ�ܿ������
l=100e3;%��·����
if transmissionLineType
    [Z,Y]=overheadLineParameterCaculation;
else
    [Z,Y]=cableParameterCaculation;
end
%% ʸ�������ز�������
%%H������ò���
global NpMaxH MRPRtolerance
NpMaxH=19;%ÿ��ģ����󼫵���Ŀ
toleranceH=0.2e-2;%��0.2%�ο���pscad��
MRPRtolerance=100;%��������������ֵ
%%Yc��ϲ�������
NpMaxYc=50;%��󼫵���Ŀ
toleranceYc=0.2e-2;
%%��ͼ
paint=1;%ֵΪ1ʱ���ƣ�0�����ͼ��
%% ִ����Ͻű� 
curveFitting
%% ʱ�����
tic
%�׵粨
[t,VA,VB,IA,IB]=doSimulation(HFit,YcFit,tau);
% [t,VA,VB,IA,IB]=doSimulation(HFitPscad,YcFitPscad,tauPscad);
%����ӵع���
% [t,VA,VB,VC,IA,IB,IC]=faultSimulation(HFit,YcFit,tau);
toc
%%simulationSet.m ��ʱ�������ļ�
%%
figure
plot(t,VA(3,:))
hold on
plot(t,VB(3,:))
hold off
%% ���������㲢��ͼ
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
    legend('�Ա����','Pscad')
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
    
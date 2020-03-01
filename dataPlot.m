 function [ MFE,RMSE ] = dataPlot(Data,F)
%   dataPlot(Data,F)
%绘图函数，用于比较拟合结果与原数据
%Data是原数据（三维矩阵Nc Nc Ns），
%F是拟合结果（结构体）
%     F.P是极点行向量 
%     F.C是留数矩阵
%     F.D是常数项矩阵
%     F.T是延时向量,对Yc来说，此项为标量零。对H来说，长度与F.C三维长度相同，即F.T(i)对应于F.C(:,:,i)

%% 计算F曲线数据
Np=length(F.P);
Ns=length(F.freq);
[Nc Nc Ns]=size(Data);
s=2i*pi*F.freq;
%%
if F.T==0%如果识别出是Yc数据
    F.T=repmat(0,1,Np);%那么就将F.T置为零向量，便可以与H统一处理
end
f=zeros(size(Data));
for i=1:Ns
    for j=1:Np
        f(:,:,i)=f(:,:,i)+exp(-s(i)*F.T(j)).*F.C(:,:,j)./(s(i)-F.P(j));
    end
    f(:,:,i)=f(:,:,i)+F.D; 
end
%%
% close(1)
% close(2)
% close(3)
%画幅值的函数
figure;
for i=1:Nc
    for j=1:Nc
        semilogx(F.freq,abs(squeeze(Data(i,j,:))),'b'); 
        hold on
        semilogx(F.freq,abs(squeeze(f(i,j,:))),'r.');
    end
end
hold off
xlabel('频率 Hz')
ylabel('幅值 ')
legend('原数据','拟合结果')
xlim([F.freq(1) F.freq(end)])
%% 画误差的函数
figure;
for i=1:Nc
    for j=1:Nc
        residualError=abs(squeeze(Data(i,j,:)-f(i,j,:)));%每个点的残差
%         residualError=abs(squeeze(abs(Data(i,j,:))-abs(f(i,j,:))));%每个点的残差，是模值的残差
        MFE(i,j)=max(residualError)/max(max(max(abs(squeeze(Data)),[],3)));%找出最大的拟合误差并储存
%         MFE(i,j)=max(residualError)/max(abs(Data(i,j,:)));%找出最大的拟合误差并储存
        SSR=sum(residualError.^2);%残差平方和
        RMSE(i,j)=sqrt(SSR/length(residualError));%均方根误差，并储存
%         semilogx(F.freq,residualError./max(abs(squeeze(Data(i,j,:)))),'r'); %%相对误差
        semilogx(F.freq,residualError./max(max(max(abs(squeeze(Data)),[],3))),'r'); %%相对误差,分母是整个H的最大模值
        hold on
    end
end
title('相对误差')
xlabel('频率 Hz')
ylabel(' ')
xlim([F.freq(1) F.freq(end)])
hold off   
%% 画相角的函数
figure;
for i=1:Nc
    for j=1:Nc
        semilogx(F.freq,180*unwrap(angle(squeeze(Data(i,j,:))))/pi,'b'); 
        hold on
        semilogx(F.freq,180*unwrap(angle(squeeze(f(i,j,:))))/pi,'r.');
    end
end
hold off
xlabel('频率 Hz')
ylabel('相角 度')
legend('原数据','拟合结果')
xlim([F.freq(1) F.freq(end)])
end


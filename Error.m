function [ F ] = Error( F,Data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% 计算F曲线数据
[Nc,Nc,Ns]=size(Data);
Np=length(F.P);
% Ns=length(F.freq);
s=2i*pi*F.freq;
f=zeros(size(Data));
if F.T==0%如果识别出是Yc与HmMps数据
    F.T=repmat(0,1,Np);%那么就将F.T置为零向量，便可以与H统一处理
end
for i=1:Ns
    for j=1:Np
        f(:,:,i)=f(:,:,i)+exp(-s(i)*F.T(j)).*F.C(:,:,j)./(s(i)-F.P(j));
    end
    f(:,:,i)=f(:,:,i)+F.D; 
end
for i=1:Nc
    for j=1:Nc
        %这里要计算RMSE
        residualError=abs(squeeze(Data(i,j,:)-f(i,j,:)));%每个点的残差
        MFE(i,j)=max(residualError)/max(max(max(abs(squeeze(Data)),[],3)));%找出最大的拟合误差并储存,依据pscad定义
        SSR=sum(residualError.^2);%残差平方和
        RMSE(i,j)=sqrt(SSR/length(residualError));%均方根误差，并储存
%         semilogx(F.freq,residualError,'r'); %%残差
        %这里要计算最大拟合误差MFE
%         hold on
    end
end
%找出最大的最大拟合误差并输出到拟合结构体
F.MFE=max(max(MFE));
F.RMSE=max(max(RMSE));%所有元素最大的均方误差输出到结构体

% residualError=abs(squeeze(Data(1,1,:)-f(1,1,:)));%每个点的残差
% F.MFE=max(residualError)/max(abs(squeeze(Data(1,1,:))));%找出最大的拟合误差并储存
% SSR=sum(residualError.^2);%残差平方和
% F.RMSE=sqrt(SSR/length(residualError));%均方根误差，并储存
end


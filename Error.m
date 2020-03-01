function [ F ] = Error( F,Data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% ����F��������
[Nc,Nc,Ns]=size(Data);
Np=length(F.P);
% Ns=length(F.freq);
s=2i*pi*F.freq;
f=zeros(size(Data));
if F.T==0%���ʶ�����Yc��HmMps����
    F.T=repmat(0,1,Np);%��ô�ͽ�F.T��Ϊ���������������Hͳһ����
end
for i=1:Ns
    for j=1:Np
        f(:,:,i)=f(:,:,i)+exp(-s(i)*F.T(j)).*F.C(:,:,j)./(s(i)-F.P(j));
    end
    f(:,:,i)=f(:,:,i)+F.D; 
end
for i=1:Nc
    for j=1:Nc
        %����Ҫ����RMSE
        residualError=abs(squeeze(Data(i,j,:)-f(i,j,:)));%ÿ����Ĳв�
        MFE(i,j)=max(residualError)/max(max(max(abs(squeeze(Data)),[],3)));%�ҳ��������������,����pscad����
        SSR=sum(residualError.^2);%�в�ƽ����
        RMSE(i,j)=sqrt(SSR/length(residualError));%��������������
%         semilogx(F.freq,residualError,'r'); %%�в�
        %����Ҫ�������������MFE
%         hold on
    end
end
%�ҳ��������������������Ͻṹ��
F.MFE=max(max(MFE));
F.RMSE=max(max(RMSE));%����Ԫ�����ľ������������ṹ��

% residualError=abs(squeeze(Data(1,1,:)-f(1,1,:)));%ÿ����Ĳв�
% F.MFE=max(residualError)/max(abs(squeeze(Data(1,1,:))));%�ҳ��������������
% SSR=sum(residualError.^2);%�в�ƽ����
% F.RMSE=sqrt(SSR/length(residualError));%��������������
end


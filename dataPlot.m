 function [ MFE,RMSE ] = dataPlot(Data,F)
%   dataPlot(Data,F)
%��ͼ���������ڱȽ���Ͻ����ԭ����
%Data��ԭ���ݣ���ά����Nc Nc Ns����
%F����Ͻ�����ṹ�壩
%     F.P�Ǽ��������� 
%     F.C����������
%     F.D�ǳ��������
%     F.T����ʱ����,��Yc��˵������Ϊ�����㡣��H��˵��������F.C��ά������ͬ����F.T(i)��Ӧ��F.C(:,:,i)

%% ����F��������
Np=length(F.P);
Ns=length(F.freq);
[Nc Nc Ns]=size(Data);
s=2i*pi*F.freq;
%%
if F.T==0%���ʶ�����Yc����
    F.T=repmat(0,1,Np);%��ô�ͽ�F.T��Ϊ���������������Hͳһ����
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
%����ֵ�ĺ���
figure;
for i=1:Nc
    for j=1:Nc
        semilogx(F.freq,abs(squeeze(Data(i,j,:))),'b'); 
        hold on
        semilogx(F.freq,abs(squeeze(f(i,j,:))),'r.');
    end
end
hold off
xlabel('Ƶ�� Hz')
ylabel('��ֵ ')
legend('ԭ����','��Ͻ��')
xlim([F.freq(1) F.freq(end)])
%% �����ĺ���
figure;
for i=1:Nc
    for j=1:Nc
        residualError=abs(squeeze(Data(i,j,:)-f(i,j,:)));%ÿ����Ĳв�
%         residualError=abs(squeeze(abs(Data(i,j,:))-abs(f(i,j,:))));%ÿ����Ĳв��ģֵ�Ĳв�
        MFE(i,j)=max(residualError)/max(max(max(abs(squeeze(Data)),[],3)));%�ҳ��������������
%         MFE(i,j)=max(residualError)/max(abs(Data(i,j,:)));%�ҳ��������������
        SSR=sum(residualError.^2);%�в�ƽ����
        RMSE(i,j)=sqrt(SSR/length(residualError));%��������������
%         semilogx(F.freq,residualError./max(abs(squeeze(Data(i,j,:)))),'r'); %%������
        semilogx(F.freq,residualError./max(max(max(abs(squeeze(Data)),[],3))),'r'); %%������,��ĸ������H�����ģֵ
        hold on
    end
end
title('������')
xlabel('Ƶ�� Hz')
ylabel(' ')
xlim([F.freq(1) F.freq(end)])
hold off   
%% ����ǵĺ���
figure;
for i=1:Nc
    for j=1:Nc
        semilogx(F.freq,180*unwrap(angle(squeeze(Data(i,j,:))))/pi,'b'); 
        hold on
        semilogx(F.freq,180*unwrap(angle(squeeze(f(i,j,:))))/pi,'r.');
    end
end
hold off
xlabel('Ƶ�� Hz')
ylabel('��� ��')
legend('ԭ����','��Ͻ��')
xlim([F.freq(1) F.freq(end)])
end


function [ tauLeft,tauMPS, tauLarge ] = tauExtraction( omega,Hm,minusGammaL,threshold,l)
%[ tau ] = tauExtraction( Hm, gamma )
%   ������ȡģ���е���ʱ��
%   tauMPS��������С��λϵͳ��ȡ������ʱ�䣬�����������������ʱ���Ż�������
%   tauLarge�ǲ�������λ��������Ĵ���ʱ�䣬�����������������tauMPS������ʱ���Ż�������
%   omega�ǽ�Ƶ��,������
%   Hm����ά���󣬱�ʾģ������СNc��Nc
%   minusGammaL,��-��l����-�̣�YZ��l������ֵ����ά������Hmͬ��С
%   error ��ֵ����ʾȡ��ֵ�½���errorʱ��Ƶ���������
%   �ο�Gustavsen��2017��optimal time delay extraction for transmission line
%   modeling,���÷�ֱ�Ӽ��㷨
%   l����·����

[Nc Nc]=size(Hm(:,:,1));
tauMPS=zeros(1,Nc);
c=3e8;%����
for i=1:Nc %�������ÿ��ģ̬
    %%�ҳ�����Ƶ���_k
    k=find(abs(Hm(i,i,:))<=threshold);
    if length(k)==0
        k=length(Hm(i,i,:))-1;
        disp(['��',num2str(i),'ģ����ֵδ���½������辫�ȣ��������Ƶ�ʼ����'])
%         distance=log(omega(end)/omega(end-1));
%         A=[1,0,0;1,distance,distance;1,2*distance,2*distance];
%         b=zeros(3,1);
%         b(1)=log(abs(Hm(i,i,end-2))/abs(Hm(i,i,end-3)))/distance;
%         b(2)=log(abs(Hm(i,i,end-1))/abs(Hm(i,i,end-2)))/distance;
%         b(3)=log(abs(Hm(i,i,end))/abs(Hm(i,i,end-1)))/distance;
%         x=A\b;
%         f=@(wk)x(1)+x(2)*log(wk/omega(end-3))+x(3)*log(wk/omega(end-3))^2;
%         k=1;
%         while k>0
%             wk=exp(k*distance)*omega(end);
%             
    else
        k=k(1);%��ѡ�ļ���Ƶ������ֵ
        if k==length(Hm(i,i,:))
            k=length(Hm(i,i,:))-1;
        end
    end
  %% ����С���
    Ak1=log(abs(Hm(i,i,k+1))/abs(Hm(i,i,k)))/log(omega(k+1)/omega(k));
    sum1=pi/2*Ak1;%�ȱ����һ��Ak�����治������
    for j=1:length(omega)-1 %ÿ��Ƶ�ʶ�Ҫ�㣬������ۼӡ�����ĳ�4��Ƶ��
%         %���Ƕ�Hm��ĳƵ�㴦�����Ľ������о�����2019/3/4�޸���Hm���㷽ʽ��ע�͵�
%         if j==77
%             temp=Hm(1,1,j);
%             Hm(1,1,j)=Hm(2,2,j);
%             Hm(2,2,j)=temp;
%         end
        Aj1=log(abs(Hm(i,i,j+1))/abs(Hm(i,i,j)))/log(omega(j+1)/omega(j));
        Bj1=log(coth(1/2*abs(log((omega(j+1)+omega(j))/2/omega(k)))));
        sum1=sum1+1/pi*(Aj1-Ak1)*Bj1*log(omega(j+1)/omega(j));
    end 
    phik1=sum1;%�ۼӵõ��ؽ������
    tauMPS(i)=(phik1-imag(minusGammaL(i,i,k)))/omega(k);%
    tauLarge(i)=imag(minusGammaL(i,i,k))/-omega(k);
    %% ��ģ����ֵ���������޸ĺ�õ��Ħ�
%     HmAuxiliary=[squeeze(Hm(i,i,1:k)).',squeeze(Hm(i,i,k+1:end)).'.*(omega(k)./omega(k+1:end)).^Np];%����Ͻ������
%     Ak2=log(abs(HmAuxiliary(k+1))/abs(HmAuxiliary(k)))/log(omega(k+1)/omega(k));
%     sum2=pi/2*Ak2;%�ȱ����һ��Ak�����治������
%     for j=1:length(omega)-1 %ÿ��Ƶ�ʶ�Ҫ�㣬������ۼӡ�����ĳ�4��Ƶ��
%         Aj2=log(abs(HmAuxiliary(j+1))/abs(HmAuxiliary(j)))/log(omega(j+1)/omega(j));
%         Bj2=log(coth(1/2*abs(log((omega(j+1)+omega(j))/2/omega(k)))));
%         sum2=sum2+1/pi*(Aj2-Ak2)*Bj2*log(omega(j+1)/omega(j));
%     end
%     phik2=sum2;%�ۼӵõ��ؽ������
%     tauLeft(i)=(phik2-imag(minusGammaL(i,i,k)))/omega(k);
%     if tauLeft(i)<l/c
        tauLeft(i)=l/c;
%     end
end



end


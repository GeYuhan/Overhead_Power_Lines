function [ tauLeft,tauMPS, tauLarge ] = tauExtraction( omega,Hm,minusGammaL,threshold,l)
%[ tau ] = tauExtraction( Hm, gamma )
%   用来提取模量中的延时τ
%   tauMPS是依据最小相位系统提取出传播时间，行相量或标量。用于时间优化的下限
%   tauLarge是不经过相位补偿算出的传播时间，行向量或标量。大于tauMPS，用于时间优化的上限
%   omega是角频率,行向量
%   Hm是三维矩阵，表示模量，大小Nc乘Nc
%   minusGammaL,即-γl，是-√（YZ）l的特征值，三维矩阵，与Hm同大小
%   error 阈值，表示取幅值下降到error时的频点来计算τ
%   参考Gustavsen的2017年optimal time delay extraction for transmission line
%   modeling,采用非直接计算法
%   l是线路长度

[Nc Nc]=size(Hm(:,:,1));
tauMPS=zeros(1,Nc);
c=3e8;%光速
for i=1:Nc %逐个计算每个模态
    %%找出计算频点ω_k
    k=find(abs(Hm(i,i,:))<=threshold);
    if length(k)==0
        k=length(Hm(i,i,:))-1;
        disp(['第',num2str(i),'模量幅值未能下降到所设精度，将以最高频率计算τ'])
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
        k=k(1);%所选的计算频点索引值
        if k==length(Hm(i,i,:))
            k=length(Hm(i,i,:))-1;
        end
    end
  %% 求最小相角
    Ak1=log(abs(Hm(i,i,k+1))/abs(Hm(i,i,k)))/log(omega(k+1)/omega(k));
    sum1=pi/2*Ak1;%先保存第一项Ak，后面不断增加
    for j=1:length(omega)-1 %每个频率都要算，结果并累加。后面改成4个频程
%         %这是对Hm再某频点处发生的交换进行纠正，2019/3/4修改了Hm计算方式后，注释掉
%         if j==77
%             temp=Hm(1,1,j);
%             Hm(1,1,j)=Hm(2,2,j);
%             Hm(2,2,j)=temp;
%         end
        Aj1=log(abs(Hm(i,i,j+1))/abs(Hm(i,i,j)))/log(omega(j+1)/omega(j));
        Bj1=log(coth(1/2*abs(log((omega(j+1)+omega(j))/2/omega(k)))));
        sum1=sum1+1/pi*(Aj1-Ak1)*Bj1*log(omega(j+1)/omega(j));
    end 
    phik1=sum1;%累加得到重建的相角
    tauMPS(i)=(phik1-imag(minusGammaL(i,i,k)))/omega(k);%
    tauLarge(i)=imag(minusGammaL(i,i,k))/-omega(k);
    %% 对模量幅值函数进行修改后得到的τ
%     HmAuxiliary=[squeeze(Hm(i,i,1:k)).',squeeze(Hm(i,i,k+1:end)).'.*(omega(k)./omega(k+1:end)).^Np];%根拟合阶数相关
%     Ak2=log(abs(HmAuxiliary(k+1))/abs(HmAuxiliary(k)))/log(omega(k+1)/omega(k));
%     sum2=pi/2*Ak2;%先保存第一项Ak，后面不断增加
%     for j=1:length(omega)-1 %每个频率都要算，结果并累加。后面改成4个频程
%         Aj2=log(abs(HmAuxiliary(j+1))/abs(HmAuxiliary(j)))/log(omega(j+1)/omega(j));
%         Bj2=log(coth(1/2*abs(log((omega(j+1)+omega(j))/2/omega(k)))));
%         sum2=sum2+1/pi*(Aj2-Ak2)*Bj2*log(omega(j+1)/omega(j));
%     end
%     phik2=sum2;%累加得到重建的相角
%     tauLeft(i)=(phik2-imag(minusGammaL(i,i,k)))/omega(k);
%     if tauLeft(i)<l/c
        tauLeft(i)=l/c;
%     end
end



end


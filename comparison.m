%% 自编电缆程序与pscad结果对比
YMag=[reshape(YMag11,1,1,101),reshape(YMag12,1,1,101),reshape(YMag13,1,1,101);reshape(YMag21,1,1,101),reshape(YMag22,1,1,101),reshape(YMag23,1,1,101);reshape(YMag31,1,1,101),reshape(YMag32,1,1,101),reshape(YMag33,1,1,101)];
YAng=[reshape(YAng11,1,1,101),reshape(YAng12,1,1,101),reshape(YAng13,1,1,101);reshape(YAng21,1,1,101),reshape(YAng22,1,1,101),reshape(YAng23,1,1,101);reshape(YAng31,1,1,101),reshape(YAng32,1,1,101),reshape(YAng33,1,1,101)];
ZMag=[reshape(ZMag11,1,1,101),reshape(ZMag12,1,1,101),reshape(ZMag13,1,1,101);reshape(ZMag21,1,1,101),reshape(ZMag22,1,1,101),reshape(ZMag23,1,1,101);reshape(ZMag31,1,1,101),reshape(ZMag32,1,1,101),reshape(ZMag33,1,1,101)];
ZAng=[reshape(ZAng11,1,1,101),reshape(ZAng12,1,1,101),reshape(ZAng13,1,1,101);reshape(ZAng21,1,1,101),reshape(ZAng22,1,1,101),reshape(ZAng23,1,1,101);reshape(ZAng31,1,1,101),reshape(ZAng32,1,1,101),reshape(ZAng33,1,1,101)];
%%
HMag=[reshape(HMag11,1,1,101),reshape(HMag12,1,1,101),reshape(HMag13,1,1,101);reshape(HMag21,1,1,101),reshape(HMag22,1,1,101),reshape(HMag23,1,1,101);reshape(HMag31,1,1,101),reshape(HMag32,1,1,101),reshape(HMag33,1,1,101)];
HMagfit=[reshape(HMagfit11,1,1,101),reshape(HMagfit12,1,1,101),reshape(HMagfit13,1,1,101);reshape(HMagfit21,1,1,101),reshape(HMagfit22,1,1,101),reshape(HMagfit23,1,1,101);reshape(HMagfit31,1,1,101),reshape(HMagfit32,1,1,101),reshape(HMagfit33,1,1,101)];
HAng=[reshape(HAng11,1,1,101),reshape(HAng12,1,1,101),reshape(HAng13,1,1,101);reshape(HAng21,1,1,101),reshape(HAng22,1,1,101),reshape(HAng23,1,1,101);reshape(HAng31,1,1,101),reshape(HAng32,1,1,101),reshape(HAng33,1,1,101)];
HAngfit=[reshape(HAngfit11,1,1,101),reshape(HAngfit12,1,1,101),reshape(HAngfit13,1,1,101);reshape(HAngfit21,1,1,101),reshape(HAngfit22,1,1,101),reshape(HAngfit23,1,1,101);reshape(HAngfit31,1,1,101),reshape(HAngfit32,1,1,101),reshape(HAngfit33,1,1,101)];
HmMag=[reshape(HmMag1,1,1,101),zeros(1,1,101),zeros(1,1,101);zeros(1,1,101),reshape(HmMag2,1,1,101),zeros(1,1,101);zeros(1,1,101),zeros(1,1,101),reshape(HmMag3,1,1,101)];
HmAng=[reshape(HmAng1,1,1,101),zeros(1,1,101),zeros(1,1,101);zeros(1,1,101),reshape(HmAng2,1,1,101),zeros(1,1,101);zeros(1,1,101),zeros(1,1,101),reshape(HmAng3,1,1,101)];
HmMagfit=[reshape(HmMagfit1,1,1,101),zeros(1,1,101),zeros(1,1,101);zeros(1,1,101),reshape(HmMagfit2,1,1,101),zeros(1,1,101);zeros(1,1,101),zeros(1,1,101),reshape(HmMagfit3,1,1,101)];
HmAngfit=[reshape(HmAngfit1,1,1,101),zeros(1,1,101),zeros(1,1,101);zeros(1,1,101),reshape(HmAngfit2,1,1,101),zeros(1,1,101);zeros(1,1,101),zeros(1,1,101),reshape(HmAngfit3,1,1,101)];
YcMag=[reshape(YcMag11,1,1,101),reshape(YcMag12,1,1,101),reshape(YcMag13,1,1,101);reshape(YcMag21,1,1,101),reshape(YcMag22,1,1,101),reshape(YcMag23,1,1,101);reshape(YcMag31,1,1,101),reshape(YcMag32,1,1,101),reshape(YcMag33,1,1,101)];
YcMagfit=[reshape(YcMagfit11,1,1,101),reshape(YcMagfit12,1,1,101),reshape(YcMagfit13,1,1,101);reshape(YcMagfit21,1,1,101),reshape(YcMagfit22,1,1,101),reshape(YcMagfit23,1,1,101);reshape(YcMagfit31,1,1,101),reshape(YcMagfit32,1,1,101),reshape(YcMagfit33,1,1,101)];
YcAng=[reshape(YcAng11,1,1,101),reshape(YcAng12,1,1,101),reshape(YcAng13,1,1,101);reshape(YcAng21,1,1,101),reshape(YcAng22,1,1,101),reshape(YcAng23,1,1,101);reshape(YcAng31,1,1,101),reshape(YcAng32,1,1,101),reshape(YcAng33,1,1,101)];
YcAngfit=[reshape(YcAngfit11,1,1,101),reshape(YcAngfit12,1,1,101),reshape(YcAngfit13,1,1,101);reshape(YcAngfit21,1,1,101),reshape(YcAngfit22,1,1,101),reshape(YcAngfit23,1,1,101);reshape(YcAngfit31,1,1,101),reshape(YcAngfit32,1,1,101),reshape(YcAngfit33,1,1,101)];
for i=1:101
    HPscad(:,:,i)=HMag(:,:,i)*exp(1i*HAng(:,:,i)*pi/180);
    YcPscad(:,:,i)=YcMag(:,:,i)*exp(1i*YcAng(:,:,i)*pi/180);
end
% save dataComparing FN freq HMag HMagfit HAng HAngfit HmMag HmMagfit HmAng HmAngfit YcMag YcAng Y Z YMag YAng ZMag ZAng HPscad YcPscad
%%
load dataComparing
%% Y幅值比较
figure(11)
for i=1:3
    for j=1:3
        subplot(3,3,j+3*(i-1));
        yyaxis left
        semilogx(FN,shiftdim(YMag(i,j,:)),'b');%pscad,对应频率FN
        hold on
        semilogx(freq,shiftdim(abs(Y(i,j,:))),'g');%，自编程序，对应频率freq
        hold off
        ylabel('幅值');
        yyaxis right
        YMagError(i,j,:)=abs( YMag(i,j,:)-abs(Y(i,j,:)) )./abs(YMag(i,j,:));
        semilogx(freq,shiftdim(YMagError(i,j,:)),'r');
        ylabel('相对误差')
        xlim([freq(1),freq(end)]);
    end
end
suptitle('Y幅值');
%% Y相角比较
figure(14)
for i=1:3
    for j=1:3
        subplot(3,3,j+3*(i-1));
        yyaxis left
        semilogx(FN,shiftdim(YAng(i,j,:)),'b');%pscad
        hold on
        semilogx(freq,shiftdim(180/pi*angle(Y(i,j,:))),'g');
        hold off
        ylabel('相角/度');
        yyaxis right
        YAngError(i,j,:)=abs( YAng(i,j,:)-180/pi*angle(Y(i,j,:)) );
        semilogx(freq,shiftdim(YAngError(i,j,:)),'r');
        ylabel('绝对误差')
        xlim([freq(1),freq(end)]);
    end
end
suptitle('Y相角');
%% Z幅值比较
figure(12)
for i=1:3
    for j=1:3
        subplot(3,3,j+3*(i-1));
        yyaxis left
        semilogx(FN,shiftdim(ZMag(i,j,:)),'b');%pscad
        hold on
        semilogx(freq,shiftdim(abs(Z(i,j,:))),'g');
        hold off
        ylabel('幅值');
        yyaxis right
        ZMagError(i,j,:)=abs( ZMag(i,j,:)-abs(Z(i,j,:)) )./abs(ZMag(i,j,:));
        semilogx(freq,shiftdim(ZMagError(i,j,:)),'r');
        ylabel('相对误差')
        xlim([freq(1),freq(end)]);
    end
end
suptitle('Z幅值');
%% Z相角比较
figure(13)
for i=1:3
    for j=1:3
        subplot(3,3,j+3*(i-1));
        yyaxis left
        semilogx(FN,shiftdim(ZAng(i,j,:)),'b');%pscad
        hold on
        semilogx(freq,shiftdim(180/pi*angle(Z(i,j,:))),'g');
        hold off
        ylabel('相角/度');
        yyaxis right
        ZAngError(i,j,:)=abs( ZAng(i,j,:)-180/pi*angle(Z(i,j,:)) );
        semilogx(freq,shiftdim(ZAngError(i,j,:)),'r');
        ylabel('绝对误差')
        xlim([freq(1),freq(end)]);
    end
end
suptitle('Z相角');
%% 特征导纳Yc与传播矩阵(H与Hm)计算
l=100e3;
[H, Hm, Ti, gamma, minusGammaL,Yc]=deal(zeros(size(Y)));
omega=2*pi*freq;%计算ω
s=omega*1i;%s==jω
for i=1:length(s)
    sqrtYZ=sqrtm(Y(:,:,i)*Z(:,:,i));%计算√YZ
    Yc(:,:,i)=sqrtYZ/Z(:,:,i);%Yc=√YZ Z^-1
    exponent=-sqrtYZ*l;%exp（A）中的A
    H(:,:,i)=expm(exponent);%得到传播矩阵
    [Ti(:,:,i), gamma(:,:,i)]=eig(sqrtYZ);%√(YZ )的特征向量矩阵（即电流相模变换矩阵）和特征值γ
    minusGammaL(:,:,i)=-gamma(:,:,i)*l;%，即-γl，是exponet的指数，延时τ的提取要用
    Hm(:,:,i)=expm(minusGammaL(:,:,i));%得到模传播矩阵,Hm=expm(-γl),这是公式
end
%% H比较
figure(1)
for i=1:3
    for j=1:3
        subplot(3,3,j+3*(i-1));
        plot(shiftdim(abs(H(i,j,:))));
        hold  on
        plot(shiftdim(HMag(i,j,:)));
%         plot(shiftdim(abs(HPscad(i,j,:))));
        hold off
        legend('自计算','Pscad','组装')
    end
end
suptitle('H幅值')
figure(2)
for i=1:3
    for j=1:3
        subplot(3,3,j+3*(i-1));
        plot(shiftdim(unwrap(180/pi*angle(H(i,j,:)))));
        hold on
        plot(shiftdim(unwrap(HAng(i,j,:))));
%         plot(shiftdim(unwrap(180/pi*angle(HPscad(i,j,:)))));
        hold off
        legend('自计算','Pscad','组装')
    end
end
suptitle('H相角')
%% Yc比较
figure(3)
for i=1:3
    for j=1:3
        subplot(3,3,j+3*(i-1));
        plot(shiftdim(abs(Yc(i,j,:))));
        hold  on
        plot(shiftdim(YcMag(i,j,:)));
%         plot(shiftdim(abs(YcPscad(i,j,:))));
        hold off
        legend('自计算','Pscad')
    end
end
suptitle('Yc幅值')
figure(4)
for i=1:3
    for j=1:3
        subplot(3,3,j+3*(i-1));
        plot(shiftdim(unwrap(180/pi*angle(Yc(i,j,:)))));
        hold on
        plot(shiftdim(unwrap(YcAng(i,j,:))));
%         plot(shiftdim(unwrap(180/pi*angle(YcPscad(i,j,:)))));
        hold off
        legend('自计算','Pscad')
    end
end
suptitle('Yc相角')
%% Pscad拟合H比较
figure(5)
for i=1:3
    for j=1:3
        subplot(3,3,j+3*(i-1));
        yyaxis left
        plot(shiftdim(HMagfit(i,j,:)));
        hold  on
        plot(shiftdim(HMag(i,j,:)));
        hold off
        legend('Pscad拟合','Pscad')
        yyaxis right
        HMagfitError(i,j,:)=abs( HMag(i,j,:)-HMagfit(i,j,:) )./abs(max(HMag(i,j,:)));
        semilogx(freq,shiftdim(HMagfitError(i,j,:)),'r');
        ylabel('相对误差')
    end
end
suptitle('Pscad拟合H幅值')
figure(6)
for i=1:3
    for j=1:3
        subplot(3,3,j+3*(i-1));
        plot(shiftdim(unwrap(HAngfit(i,j,:))));
        hold on
        plot(shiftdim(unwrap(HAng(i,j,:))));
        hold off
        legend('Pscad拟合','Pscad')
    end
end
suptitle('Pscad拟合H相角')
%% Pscad拟合Yc比较
figure(7)
for i=1:3
    for j=1:3
        subplot(3,3,j+3*(i-1));
        yyaxis Left
        plot(shiftdim(YcMagfit(i,j,:)));
        hold  on
        plot(shiftdim(YcMag(i,j,:)));
        hold off
        legend('Pscad拟合','Pscad')
        yyaxis right
        YcMagfitError(i,j,:)=abs( YcMag(i,j,:)-YcMagfit(i,j,:) )./abs(max(YcMag(i,j,:)));
        semilogx(freq,shiftdim(YcMagfitError(i,j,:)),'r');
        ylabel('相对误差')
    end
end
suptitle('Pscad拟合Yc幅值')
figure(8)
for i=1:3
    for j=1:3
        subplot(3,3,j+3*(i-1));
        plot(shiftdim(unwrap(YcAngfit(i,j,:))));
        hold on
        plot(shiftdim(unwrap(YcAng(i,j,:))));
        hold off
        legend('Pscad拟合','Pscad')
    end
end
suptitle('Pscad拟合Yc相角')
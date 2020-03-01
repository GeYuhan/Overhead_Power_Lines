%%
tauPscad=[3.870208420195887E-004,  6.245788264517912E-004,  2.436715805576424E-003];
%%
p1=p1(2*(1:17)-1)+1i*p1(2*(1:17));
p2=p2(2*(1:23)-1)+1i*p2(2*(1:23));
p3=p3(2*(1:12)-1)+1i*p3(2*(1:12));
%%
PolesPscad{1}=p1.';
PolesPscad{2}=p2.';
PolesPscad{3}=p3.';
%%
b=b(1:2:end)+1i*b(2:2:end);
%%
clear b1 b2 b3
for i=1:9
    b1(i,:)=b(1+17*(i-1):1+17*(i-1)+16).';
    b2(i,:)=b(1+23*(i-1)+153:1+23*(i-1)+22+153).';
    b3(i,:)=b(1+12*(i-1)+207+153:1+12*(i-1)+11+207+153).';
end
%%
c=cat(3,reshape(d1,3,3,5),reshape(d2,3,3,6),reshape(d3,3,3,11));
%%
f=zeros(9,101);
% PolesPscad{1}=[-14226.6830821468,-132884.331506846,-666368.112382286,-2447692.76770130,-7493900.75818503];
% PolesPscad{2}=[-37.8520092115209,-716.795331034211,-46750.6674538162,-392972.361392121,-2585090.40120981,-51165748.4729594];
% PolesPscad{3}=[-2.89459240232613,-20.5074006560677,-112.018311415892,-527.184183274129,-2214.07346351439,-8434.91354148902,-28038.3774451594,-143817.030918459,-565773.658692029,-4074523.33962640-1279608.90154004i,-4074523.33962640+1279608.90154004i];
for i=1:9
    for j=1:length(FN)
        f(i,j)=f(i,j)+sum(b1(i,:)./(2i*pi*FN(j)-PolesPscad{1}))*exp(-2i*pi*FN(j)*tauPscad(1))+...
            sum(b2(i,:)./(2i*pi*FN(j)-PolesPscad{2}))*exp(-2i*pi*FN(j)*tauPscad(2))+...
            sum(b3(i,:)./(2i*pi*FN(j)-PolesPscad{3}))*exp(-2i*pi*FN(j)*tauPscad(3));
    end
end
%%
for i=9
    figure(7)
    semilogx(FN,abs(f(i,:)))
    xlim([FN(1),FN(end)])
    figure(8)
    semilogx(FN,180/pi*unwrap(angle(f(i,:))))
    xlim([FN(1),FN(end)])
end
%%
figure(9)
hold on
for i=1:3
    for j=1:3
%         figu
re(9)
        semilogx(FN,squeeze(abs(H(i,j,:))));
        xlim([FN(1),FN(end)])
%         figure(10)
%         semilogx(FN,180/pi*unwrap(angle(squeeze(H(i,j,:)))));
%         xlim([FN(1),FN(end)])
    end
end
%%
HFitPscad.freq=FN;
HFitPscad.P=[PolesPscad{1},PolesPscad{2},PolesPscad{3}];
HFitPscad.C=c;
HFitPscad.D=zeros(3);
HFitPscad.T=[repmat(tauPscad(1),1,5),repmat(tauPscad(2),1,6),repmat(tauPscad(3),1,11)];
HFitPscad.MRPR=1.22;
HFitPscad.MFE=0.001775;   
HFitPscad.RMSE=0.000424 ;
%%
HFit=Hresidue(FN,tauPscad,PolesPscad,H(1,3,:),0);
dataPlot(H(1,3,:),HFit);%绘图
%%
[SER,rmserr,bigHfit]=VFdriver(Hm(2,2,:),2i*pi*FN,PolesPscad{2},opts)
%%
HFit.C(1,1,:)=reshape([d1(1,:),d2(1,:),d3(1,:)],1,1,22);
HFit.C(1,2,:)=reshape([d1(4,:),d2(4,:),d3(4,:)],1,1,22);
HFit.C(1,3,:)=reshape([d1(7,:),d2(7,:),d3(7,:)],1,1,22);
HFit.C(2,1,:)=reshape([d1(2,:),d2(2,:),d3(2,:)],1,1,22);
HFit.C(2,2,:)=reshape([d1(5,:),d2(5,:),d3(5,:)],1,1,22);
HFit.C(2,3,:)=reshape([d1(8,:),d2(8,:),d3(8,:)],1,1,22);
HFit.C(3,1,:)=reshape([d1(3,:),d2(3,:),d3(3,:)],1,1,22);
HFit.C(3,2,:)=reshape([d1(6,:),d2(6,:),d3(6,:)],1,1,22);
HFit.C(3,3,:)=reshape([d1(9,:),d2(9,:),d3(9,:)],1,1,22);
dataPlot(H,HFit);%绘图
%%
m=4;
HFit.C=reshape([b1(m,:),b2(m,:),b3(m,:)],1,1,52);
dataPlot(H(1,2,:),HFit);%绘图
%% 拟合pscad拟合出的曲线
HFit=Hresidue(FN,tauPscad,PolesPscad,reshape(f,3,3,101),0);
dataPlot(reshape(f,3,3,101),HFit)%绘图
%% 整体自己拟合
HFit=Hresidue(FN,tauPscad,PolesPscad,H,0);
dataPlot(H,HFit);%绘图
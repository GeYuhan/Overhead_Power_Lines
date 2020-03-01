%% 导入pscad数据
%%优化后的延时
tauPscad=[  6.245747470990270E-004   2.201214748921703E-003   5.494120319644245E-003];
%%优化后的模量
Hm1=HmMag11.*exp(1i*HmAng11/180*pi);
Hm2=HmMag22.*exp(1i*HmAng11/180*pi);
Hm3=HmMag33.*exp(1i*HmAng11/180*pi);
Hm(1,1,:)=reshape(Hm1.*exp(-2i*pi*FN*tauPscad(1)),1,1,101);%恢复原相角
Hm(2,2,:)=reshape(Hm2.*exp(-2i*pi*FN*tauPscad(2)),1,1,101);
Hm(3,3,:)=reshape(Hm3.*exp(-2i*pi*FN*tauPscad(3)),1,1,101);
for i=1:101
    for j=1:3
        if Hm(j,j,i)==0
            Hm(j,j,i)=exp(-745);
        end
    end
end
%% 特征值
gamma1=LamdaMag11.*exp(1i*LamdaAng11/180*pi);
gamma2=LamdaMag22.*exp(1i*LamdaAng22/180*pi);
gamma3=LamdaMag33.*exp(1i*LamdaAng33/180*pi);
gamma(1,1,:)=reshape(sqrt(gamma1),1,1,101);
gamma(2,2,:)=reshape(sqrt(gamma2),1,1,101);
gamma(3,3,:)=reshape(sqrt(gamma3),1,1,101);
%% 相域传播函数
HMag=[reshape(HMag11,1,1,101),reshape(HMag12,1,1,101),reshape(HMag13,1,1,101);reshape(HMag21,1,1,101),reshape(HMag22,1,1,101),reshape(HMag23,1,1,101);reshape(HMag31,1,1,101),reshape(HMag32,1,1,101),reshape(HMag33,1,1,101)];
HAng=[reshape(HAng11,1,1,101),reshape(HAng12,1,1,101),reshape(HAng13,1,1,101);reshape(HAng21,1,1,101),reshape(HAng22,1,1,101),reshape(HAng23,1,1,101);reshape(HAng31,1,1,101),reshape(HAng32,1,1,101),reshape(HAng33,1,1,101)];
%% 特征导纳
YcMag=[reshape(YcMag11,1,1,101),reshape(YcMag12,1,1,101),reshape(YcMag13,1,1,101);reshape(YcMag21,1,1,101),reshape(YcMag22,1,1,101),reshape(YcMag23,1,1,101);reshape(YcMag31,1,1,101),reshape(YcMag32,1,1,101),reshape(YcMag33,1,1,101)];
YcAng=[reshape(YcAng11,1,1,101),reshape(YcAng12,1,1,101),reshape(YcAng13,1,1,101);reshape(YcAng21,1,1,101),reshape(YcAng22,1,1,101),reshape(YcAng23,1,1,101);reshape(YcAng31,1,1,101),reshape(YcAng32,1,1,101),reshape(YcAng33,1,1,101)];
%%
for j=1:101
    H(:,:,j)=HMag(:,:,j).*exp(1i*HAng(:,:,j)/180*pi);
    minusGammaL(:,:,j)=-gamma(:,:,j)*l;
    Hm(:,:,j)=expm(minusGammaL(:,:,j));
%     Yc(:,:,j)=YcMag(:,:,j).*exp(1i*YcAng(:,:,j)/180*pi);
end
%% pscad得到的留数极点模型
% PolesPscad={[-14226.6830821468,-132884.331506846,-666368.112382286,-2447692.76770130,-7493900.75818503]};
% poles=[-14226.6830821468,-132884.331506846,-666368.112382286,-2447692.76770130,-7493900.75818503];
% % for j=1:length(FN)%频率循环
% %     HmMPS(1,1,j)=Hm1(:,:,j)*exp(2i*pi*FN(j)*tauPscad(1));%MPS即 Minimum Phase System,最小相位系统,将tau替换成了tauTest(k)
% % end%频率循环结束
% opts.N=5;
% opts.Niter1=0;
% opts.Niter2=20;
% [SER,rmserr,bigHfit]=VFdriver(Hm1(1,1,:),2i*pi*FN,poles,opts);
% HmFitPscad(1).freq=freq;
% HmFitPscad(1).P=SER.poles;
% HmFitPscad(1).C=SER.R;
% HmFitPscad(1).D=0;
% HmFitPscad(1).T=0;
% HmFitPscad(1).MFE=[];
% HmFitPscad(1).RMSE=rmserr;
% % HmFitPscad(1)=Hresidue(FN,tauPscad(2),PolesPscad,Hm3,0);
%% 绘制一定极点数目下，最大RMSE、MFE随tau变化的图像
HmMFE=[];
for Np=5%极点数目递增
    tauTest=0.5*l/3e8:1e-8:3*l/3e8;
    for k=1:length(tauTest)%计算100个τ值       
        %% 模域
        %矩阵预先声明
        HmMPS=zeros(size(Hm));
        [~, Nc]=size(Hm(:,:,1));
        %逐个计算每个模量
        for i=1%应为Nc-----------------
%             fprintf('第%d个模量\n',i)
%             [tauMPS, tauLarge]=tauExtraction(omega,Hm(i,i,:),minusGammaL(i,i,:),0.001);%tau是行向量，选择当幅值下降到targetError时的频点
            %计算函数值，调用矢量拟合
            for j=1:length(s)%频率循环
                HmMPS(i,i,j)=Hm(i,i,j)*exp(s(j)*tauTest(k));%MPS即 Minimum Phase System,最小相位系统,将tau替换成了tauTest(k)
            end%频率循环结束
            %%2019/3/2
            opts.N=Np ;%           %Order of approximation. 
            opts.poletype='lincmplx'; %Mix of linearly spaced and logarithmically spaced poles
            opts.stable=1;      %不翻转不稳定极点
            opts.weightparam=1; %5 --> weighting with inverse magnitude norm
            opts.Niter1=7;    %Number of iterations for fitting sum of elements (fast!)
            opts.Niter2=4;    %Number of iterations for matrix fitting
            opts.asymp=1;      %Fitting includes D
            opts.plot=0;       %不绘图
            opts.screen=0;      %拟合过程中不输出提示
            poles=[];
            [SER,rmserr,bigHfit]=VFdriver(HmMPS(i,i,:),s,poles,opts);
            HmFit(i).freq=freq;
            HmFit(i).P=SER.poles;
            HmFit(i).C=SER.R;
            HmFit(i).D=0;
            HmFit(i).T=0;
            HmFit(i).MFE=[];
            HmFit(i).RMSE=rmserr;
            %%2019/3/2
%             HmFit(i)=matrixFit(freq,HmMPS(i,i,:),option);%HmFit是结构体数组,HmFit(i).P存着极点
            HmFit(i)=HmError(HmFit(i),tauTest(k),Hm(i,i,:));%最小相位的误差与模量的误差是一样的
%             %设置拟合选项（当前只有pattern、Np俩项）
%             option.pattern=0;
%             option.Np=Np;
%             option.relax=1;
%             HmFit(i)=matrixFit(freq,HmMPS(i,i,:),option);%HmFit是结构体数组,HmFit(i).P存着极点
            HmMFE(k,i,Np)=HmFit(i).MFE;%将均方误差赋给fu
            HmRMSE(k,i,Np)=HmFit(i).RMSE;
            %下面是验证Gustavsen的左括号实验，结果3000km没啥效果，与tauMPS完全重叠
%             [ tauMPS, tauLarge, tauRevised ] = tauExtraction( omega,Hm,minusGammaL,0.001,Np);
%             tauNew(:,Np)=tauRevised;
        end%逐个计算摸个模量结束
%         HPoles=cat(1,HmFit.P);%将极点拼接成矩阵，每一行对应一个τ，用于留数计算
%         %% 相域
%         HFit=Hresidue(freq,tau,HPoles,H);
%         HRMSE(k,:,Np)=HFit.RMSE;
%         HMFE(k,:,Np)=HFit.MFE;

    end
end
%%
HPoles{1}=[HmFit(1).P];%将极点拼接为元胞矩阵供相域拟合用
HPoles{2}=[HmFit(2).P];%将极点拼接为元胞矩阵供相域拟合用
HPoles{3}=[HmFit(3).P];%将极点拼接为元胞矩阵供相域拟合用
HFit=Hresidue(FN,tauPscad,HPoles,H,0);
% opts.N=22;
% opts.Niter1=0;
% opts.Niter2=0;
% HFit=VFdriver(H,2i*pi*FN,[HmFit(1).P;HmFit(2).P;HmFit(3).P],opts);
% %%2019/3/7
% clear a
% a.freq=FN;
% a.P=HFit.poles;
% a.C=HFit.R;
% a.D=0;
% a.T=0;
% a.MFE=[];
% a.RMSE=rmserr;
%% 绘图
figure
plot(tauTest(1:length(squeeze(HmMFE(:,i,Np)))),squeeze(HmMFE(:,i,Np)));
hold on
line([tauPscad(i),tauPscad(i)],[min(HmMFE(:,i,Np)),max(HmMFE(:,i,Np))]);
line([tau(i),tau(i)],[min(HmMFE(:,i,Np)),max(HmMFE(:,i,Np))]);
hold off
%% 绘制不同极点数目下模量随传播时间变化的图像
% figure
% for k=1:3%模量数目
%     subplot(3,1,k)
% %     semilogy(tauTest,HmRMSE(:,k,19).');
% %     hold on
% %     plot([tauNew(k,1),tauNew(k,1)],[1e-20,1],'-r');
% %     plot([tauMPS(k),tauMPS(k)],[1e-10,1],'-b');
%     for Npi=20
%     semilogy(tauTest(1:426),HmRMSE(1:426,k,Npi).');%426是第二模量的上界
%     hold on
% %     plot([tauNew(k,Npi),tauNew(k,Npi)],[1e-20,1],'-r');
% %     plot([tauMPS(k),tauMPS(k)],[1e-10,1],'-b');
%     end
%     hold off
%     string=sprintf('不同极点数目下拟合均方根误差与τ的关系',k);
%     title(string);
% %     legend('Np=5','Np=10','Np=15');
% end
% 
% % hold on
% % plot([tau tau]*1e6,[1e-7 1],'--r')
% % % xlim(1e6*[3.3e-4 tauTest(end)])
% % xlim(1e6*[tauTest(1) tauTest(end)])
% % xlabel('τ us')
% % ylabel('均方根误差')

%% ����pscad����
%%�Ż������ʱ
tauPscad=[  6.245747470990270E-004   2.201214748921703E-003   5.494120319644245E-003];
%%�Ż����ģ��
Hm1=HmMag11.*exp(1i*HmAng11/180*pi);
Hm2=HmMag22.*exp(1i*HmAng11/180*pi);
Hm3=HmMag33.*exp(1i*HmAng11/180*pi);
Hm(1,1,:)=reshape(Hm1.*exp(-2i*pi*FN*tauPscad(1)),1,1,101);%�ָ�ԭ���
Hm(2,2,:)=reshape(Hm2.*exp(-2i*pi*FN*tauPscad(2)),1,1,101);
Hm(3,3,:)=reshape(Hm3.*exp(-2i*pi*FN*tauPscad(3)),1,1,101);
for i=1:101
    for j=1:3
        if Hm(j,j,i)==0
            Hm(j,j,i)=exp(-745);
        end
    end
end
%% ����ֵ
gamma1=LamdaMag11.*exp(1i*LamdaAng11/180*pi);
gamma2=LamdaMag22.*exp(1i*LamdaAng22/180*pi);
gamma3=LamdaMag33.*exp(1i*LamdaAng33/180*pi);
gamma(1,1,:)=reshape(sqrt(gamma1),1,1,101);
gamma(2,2,:)=reshape(sqrt(gamma2),1,1,101);
gamma(3,3,:)=reshape(sqrt(gamma3),1,1,101);
%% ���򴫲�����
HMag=[reshape(HMag11,1,1,101),reshape(HMag12,1,1,101),reshape(HMag13,1,1,101);reshape(HMag21,1,1,101),reshape(HMag22,1,1,101),reshape(HMag23,1,1,101);reshape(HMag31,1,1,101),reshape(HMag32,1,1,101),reshape(HMag33,1,1,101)];
HAng=[reshape(HAng11,1,1,101),reshape(HAng12,1,1,101),reshape(HAng13,1,1,101);reshape(HAng21,1,1,101),reshape(HAng22,1,1,101),reshape(HAng23,1,1,101);reshape(HAng31,1,1,101),reshape(HAng32,1,1,101),reshape(HAng33,1,1,101)];
%% ��������
YcMag=[reshape(YcMag11,1,1,101),reshape(YcMag12,1,1,101),reshape(YcMag13,1,1,101);reshape(YcMag21,1,1,101),reshape(YcMag22,1,1,101),reshape(YcMag23,1,1,101);reshape(YcMag31,1,1,101),reshape(YcMag32,1,1,101),reshape(YcMag33,1,1,101)];
YcAng=[reshape(YcAng11,1,1,101),reshape(YcAng12,1,1,101),reshape(YcAng13,1,1,101);reshape(YcAng21,1,1,101),reshape(YcAng22,1,1,101),reshape(YcAng23,1,1,101);reshape(YcAng31,1,1,101),reshape(YcAng32,1,1,101),reshape(YcAng33,1,1,101)];
%%
for j=1:101
    H(:,:,j)=HMag(:,:,j).*exp(1i*HAng(:,:,j)/180*pi);
    minusGammaL(:,:,j)=-gamma(:,:,j)*l;
    Hm(:,:,j)=expm(minusGammaL(:,:,j));
%     Yc(:,:,j)=YcMag(:,:,j).*exp(1i*YcAng(:,:,j)/180*pi);
end
%% pscad�õ�����������ģ��
% PolesPscad={[-14226.6830821468,-132884.331506846,-666368.112382286,-2447692.76770130,-7493900.75818503]};
% poles=[-14226.6830821468,-132884.331506846,-666368.112382286,-2447692.76770130,-7493900.75818503];
% % for j=1:length(FN)%Ƶ��ѭ��
% %     HmMPS(1,1,j)=Hm1(:,:,j)*exp(2i*pi*FN(j)*tauPscad(1));%MPS�� Minimum Phase System,��С��λϵͳ,��tau�滻����tauTest(k)
% % end%Ƶ��ѭ������
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
%% ����һ��������Ŀ�£����RMSE��MFE��tau�仯��ͼ��
HmMFE=[];
for Np=5%������Ŀ����
    tauTest=0.5*l/3e8:1e-8:3*l/3e8;
    for k=1:length(tauTest)%����100����ֵ       
        %% ģ��
        %����Ԥ������
        HmMPS=zeros(size(Hm));
        [~, Nc]=size(Hm(:,:,1));
        %�������ÿ��ģ��
        for i=1%ӦΪNc-----------------
%             fprintf('��%d��ģ��\n',i)
%             [tauMPS, tauLarge]=tauExtraction(omega,Hm(i,i,:),minusGammaL(i,i,:),0.001);%tau����������ѡ�񵱷�ֵ�½���targetErrorʱ��Ƶ��
            %���㺯��ֵ������ʸ�����
            for j=1:length(s)%Ƶ��ѭ��
                HmMPS(i,i,j)=Hm(i,i,j)*exp(s(j)*tauTest(k));%MPS�� Minimum Phase System,��С��λϵͳ,��tau�滻����tauTest(k)
            end%Ƶ��ѭ������
            %%2019/3/2
            opts.N=Np ;%           %Order of approximation. 
            opts.poletype='lincmplx'; %Mix of linearly spaced and logarithmically spaced poles
            opts.stable=1;      %����ת���ȶ�����
            opts.weightparam=1; %5 --> weighting with inverse magnitude norm
            opts.Niter1=7;    %Number of iterations for fitting sum of elements (fast!)
            opts.Niter2=4;    %Number of iterations for matrix fitting
            opts.asymp=1;      %Fitting includes D
            opts.plot=0;       %����ͼ
            opts.screen=0;      %��Ϲ����в������ʾ
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
%             HmFit(i)=matrixFit(freq,HmMPS(i,i,:),option);%HmFit�ǽṹ������,HmFit(i).P���ż���
            HmFit(i)=HmError(HmFit(i),tauTest(k),Hm(i,i,:));%��С��λ�������ģ���������һ����
%             %�������ѡ���ǰֻ��pattern��Np���
%             option.pattern=0;
%             option.Np=Np;
%             option.relax=1;
%             HmFit(i)=matrixFit(freq,HmMPS(i,i,:),option);%HmFit�ǽṹ������,HmFit(i).P���ż���
            HmMFE(k,i,Np)=HmFit(i).MFE;%����������fu
            HmRMSE(k,i,Np)=HmFit(i).RMSE;
            %��������֤Gustavsen��������ʵ�飬���3000kmûɶЧ������tauMPS��ȫ�ص�
%             [ tauMPS, tauLarge, tauRevised ] = tauExtraction( omega,Hm,minusGammaL,0.001,Np);
%             tauNew(:,Np)=tauRevised;
        end%�����������ģ������
%         HPoles=cat(1,HmFit.P);%������ƴ�ӳɾ���ÿһ�ж�Ӧһ���ӣ�������������
%         %% ����
%         HFit=Hresidue(freq,tau,HPoles,H);
%         HRMSE(k,:,Np)=HFit.RMSE;
%         HMFE(k,:,Np)=HFit.MFE;

    end
end
%%
HPoles{1}=[HmFit(1).P];%������ƴ��ΪԪ���������������
HPoles{2}=[HmFit(2).P];%������ƴ��ΪԪ���������������
HPoles{3}=[HmFit(3).P];%������ƴ��ΪԪ���������������
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
%% ��ͼ
figure
plot(tauTest(1:length(squeeze(HmMFE(:,i,Np)))),squeeze(HmMFE(:,i,Np)));
hold on
line([tauPscad(i),tauPscad(i)],[min(HmMFE(:,i,Np)),max(HmMFE(:,i,Np))]);
line([tau(i),tau(i)],[min(HmMFE(:,i,Np)),max(HmMFE(:,i,Np))]);
hold off
%% ���Ʋ�ͬ������Ŀ��ģ���洫��ʱ��仯��ͼ��
% figure
% for k=1:3%ģ����Ŀ
%     subplot(3,1,k)
% %     semilogy(tauTest,HmRMSE(:,k,19).');
% %     hold on
% %     plot([tauNew(k,1),tauNew(k,1)],[1e-20,1],'-r');
% %     plot([tauMPS(k),tauMPS(k)],[1e-10,1],'-b');
%     for Npi=20
%     semilogy(tauTest(1:426),HmRMSE(1:426,k,Npi).');%426�ǵڶ�ģ�����Ͻ�
%     hold on
% %     plot([tauNew(k,Npi),tauNew(k,Npi)],[1e-20,1],'-r');
% %     plot([tauMPS(k),tauMPS(k)],[1e-10,1],'-b');
%     end
%     hold off
%     string=sprintf('��ͬ������Ŀ����Ͼ����������ӵĹ�ϵ',k);
%     title(string);
% %     legend('Np=5','Np=10','Np=15');
% end
% 
% % hold on
% % plot([tau tau]*1e6,[1e-7 1],'--r')
% % % xlim(1e6*[3.3e-4 tauTest(end)])
% % xlim(1e6*[tauTest(1) tauTest(end)])
% % xlabel('�� us')
% % ylabel('���������')

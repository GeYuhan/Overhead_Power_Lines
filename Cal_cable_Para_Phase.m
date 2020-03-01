%%%%%%%%%%%%%%%%%������²���%%%%%%%%
%% 
clc
clear all
Parameter_cable;
global s 
global h
global x_mul  h1 h2
% ff=[50 5000];
% ff=[0.01:1:99,100:10:1e3,(1e3+10):1000:1.01e6];
ff=[1:10:100];
%%%�����迹�Ĳ���
Z1=zeros(1,length(ff));
Z2=zeros(1,length(ff));
Z3=zeros(1,length(ff));
Z4=zeros(1,length(ff));
Z5=zeros(1,length(ff));
Z6=zeros(1,length(ff));
Z7=zeros(1,length(ff));
Z8=zeros(1,length(ff));
Z1m=zeros(1,length(ff));
Z2m=zeros(1,length(ff));
%%%%%����迹����
Zg1=zeros(1,length(ff));
Zg0=zeros(1,length(ff));
Zg2=zeros(1,length(ff));
Zg=zeros(1,length(ff));
%%%%��·�迹
Z11=zeros(1,length(ff));
Z22=zeros(1,length(ff));
Z12=zeros(1,length(ff));
Z21=zeros(1,length(ff));
Z23=zeros(1,length(ff));
Z32=zeros(1,length(ff));
Z33=zeros(1,length(ff));
%%%%�����迹����
ZRR=zeros(1,length(ff));
ZDD=zeros(1,length(ff));
ZMM=zeros(1,length(ff));
ZRD=zeros(1,length(ff));
ZDR=zeros(1,length(ff));
ZRM=zeros(1,length(ff));
ZMR=zeros(1,length(ff));
ZDM=zeros(1,length(ff));
ZMD=zeros(1,length(ff));
%%
for nnn=1:length(ff)
f=ff(nnn);%Hz������Ƶ��
Z2(nnn)=1j*2*pi*f*mu0/2/pi*log(r3/r2);%%���µ�о�����β�֮��
Z5(nnn)=1j*2*pi*f*mu0/2/pi*log(r5/r4);%%���β��뻥��֮��
Z8(nnn)=1j*2*pi*f*mu0/2/pi*log(r7/r6);%%��������֮��

%% ��·�迹����%%
%%%%%%%%%%%%%%%%%%��������迹%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1=sqrt(1j*2*pi*f*mu1/rou1);
if r1==0 
    Z1(nnn)=1/(2*pi*r2)*sqrt(1j*2*pi*f*mu1*rou1).*besseli(0,m1*r2)./besseli(1,m1*r2);%%ʵ�ĵ�о�����迹
else
    D1=besseli(1,m1*r2).*besselk(1,m1*r1)-besseli(1,m1*r1).*besselk(1,m1*r2);
    Z1(nnn)=rou1.*m1./(2*pi*r2*D1).*(besseli(0,m1*r2).*besselk(1,m1*r1)+besseli(1,m1*r1).*besselk(0,m1*r2));%%���ĵ�о������迹
end

m2=sqrt(1j*2*pi*f*mu2/rou2);
D2=besseli(1,m2*r4).*besselk(1,m2*r3)-besseli(1,m2*r3).*besselk(1,m2*r4);
Z3(nnn)=rou2.*m2./(2*pi*r3*D2).*(besseli(0,m2*r3).*besselk(1,m2*r4)+besseli(1,m2*r4).*besselk(0,m2*r3));%%���β��ڱ����迹
Z4(nnn)=rou2.*m2./(2*pi*r4*D2).*(besseli(0,m2*r4).*besselk(1,m2*r3)+besseli(1,m2*r3).*besselk(0,m2*r4));%%���β�������迹

m3=sqrt(1j*2*pi*f*mu3/rou3);
D3=besseli(1,m3*r6).*besselk(1,m3*r5)-besseli(1,m3*r5).*besselk(1,m3*r6);
Z6(nnn)=rou3.*m3./(2*pi*r5*D3).*(besseli(0,m3*r5).*besselk(1,m3*r6)+besseli(1,m3*r6).*besselk(0,m3*r5));%%�����ڱ����迹
Z7(nnn)=rou3.*m3./(2*pi*r6*D3).*(besseli(0,m3*r6).*besselk(1,m3*r5)+besseli(1,m3*r5).*besselk(0,m3*r6));%%����������迹

%%%%%��״�����������֮��ת���迹
Z1m(nnn)=rou2./(2*pi*r3*r4*D2);
Z2m(nnn)=rou3/(2*pi*r5*r6*D3);

%% ����Է����迹%%%%%%%
%%% ���Ƽ��㹫ʽ
me=sqrt(1j*2*pi*f*mu0/rou4);
%%%%% ���Ǵ��Ӱ���µĴ�ط����迹����%%%%%
h=ha;
s=1j*2*pi*f;%%%����Ƶ��
if Erec==1  %%%%����������������´��ڵ���ʱ,%% 
    if SOL==1 %�������������
        %%%%��ȷ���㷽��
        Zg1(nnn)=integral(@self_impedance,1e-6,1);
        %%%%���鹫ʽ
        Zg0(nnn)=1j*2*pi*f*mu0/2/pi.*(-log(0.57722*me*r7/2)+0.5-4*me*h/3);      
        %%%%%Pollaczek��ʽ%%%%%
        pp=1./sqrt(s*mu4/rou4);%%%%���������
        Zg2(nnn)=s*mu0/2/pi*(besselk(0,r7./pp)-besselk(0,2*h./pp)+2*integral(@self_impedance_Pollaczek,1e-6,1));
    elseif SOL==2   %%%˫�����������
        Zg1(nnn)=integral(@self_impedance_ground_double,1e-6,1);  
    end   
else   %%%%%%�����������������ڵ���ʱ%%%%%%%%%
    if SOL==1 %�������������
        %%%%%%��ȷ���㷽��
        Zg1(nnn)=integral(@self_impedance_air,1e-6,1);
        %%%%%%Carson��ʽ
        Zg2(nnn)=2*pi*f/pi*integral(@self_impedance_carson,1e-6,1);
    elseif SOL==2  %%%˫�����������
        Zg1(nnn)=integral(@self_impedance_air_double,1e-6,1);
    end
end

%% ��·�迹��%%%%%%%%
Zg(nnn)=Zg1(nnn);
Z11(nnn)=Z1(nnn)+Z2(nnn)+Z3(nnn);
Z22(nnn)=Z4(nnn)+Z5(nnn)+Z6(nnn);
Z33(nnn)=Z7(nnn)+Z8(nnn)+Zg(nnn);
Z12(nnn)=-Z1m(nnn);
Z21(nnn)=-Z1m(nnn);
Z23(nnn)=-Z2m(nnn);
Z32(nnn)=-Z2m(nnn);
% Zh=[Z11,Z12,0;Z21,Z22,Z23;0,Z32,Z33];
%%% �����迹��%%%%%%%%%%
ZRR(nnn)=Z11(nnn)+2*Z12(nnn)+Z22(nnn)+2*Z23(nnn)+Z33(nnn);
ZDD(nnn)=Z22(nnn)+2*Z23(nnn)+Z33(nnn);
ZMM(nnn)=Z33(nnn);
ZRD(nnn)=Z12(nnn)+Z22(nnn)+2*Z23(nnn)+Z33(nnn);
ZDR(nnn)=ZRD(nnn);
ZRM(nnn)=Z23(nnn)+Z33(nnn);
ZMR(nnn)=ZRM(nnn);
ZDM(nnn)=Z23(nnn)+Z33(nnn);
ZMD(nnn)=Z23(nnn)+Z33(nnn);
Zph(:,:,nnn)=[ZRR,ZRD,ZRM;ZDR,ZDD,ZDM;ZMR,ZMD,ZMM];%��������迹����


%%% ���ɾ���%%
CRD=2*pi*espo1/log(r3/r2);
CDM=2*pi*espo2/log(r5/r4);
CME=2*pi*espo3/log(r7/r6);
if Erec==2
    CME=2*pi*espo0/log(2*h/r7);
end
Yph(:,:,nnn)=1j*2*pi*f*[CRD,-CRD,0;-CRD,CRD+CDM,-CDM;0,-CDM,CDM+CME];%������µ��ɾ���

 %% �������ർ������%%%%%
 Zgabc=[];
 Ymul=[];
 Zg_p=zeros(1,3);
if Phase==3
    %%%%��ػ��迹�����%%%%%%%%    
    if Erec==1 %%%%%%%%������ص��µ����%%%%%%%%%%
        for ii=1:Phase
            %%%%�������Ĵ�ط����迹��
            h=habc(ii);
            if SOL==1 %�������������
                %%%%��ȷ���㷽��
                Zg_p(ii)=integral(@self_impedance,1e-6,1);    
            elseif SOL==2   %%%˫�����������
                Zg_p(ii)=integral(@self_impedance_ground_double,1e-6,1);  
            end   
            %%%%%%%��������Ļ��迹%%%%%%
            for jj=1:Phase
                if jj>ii
                    x_mul=abs(xabc(ii)-xabc(jj));
                    h1=habc(ii);
                    h2=habc(jj);
                    if SOL==1
                        Zgabc=[Zgabc,integral(@mutl_impedance,1e-4,1)];%% ���������ṹ�����
                    elseif SOL==2
                        Zgabc=[Zgabc,integral(@mutl_impedance_ground_double,1e-2,1)];%% ˫�������ṹ�����
                    end
                end                
            end
        end        
    else
        for ii=1:Phase
            %%%%�������Ĵ�ط����迹
            h=habc(ii);
             if SOL==1 %�������������
                 %%%%%%��ȷ���㷽��
                  Zg_p(ii)=integral(@self_impedance_air,1e-6,1);
             elseif SOL==2  %%%˫�����������
                  Zg_p(ii)=integral(@self_impedance_air_double,1e-6,1);
             end
             %%%%%%%��������Ļ��迹%%%%%%
            for jj=1:Phase
                if jj>ii
                    x_mul=abs(xabc(ii)-xabc(jj));
                    h1=habc(ii);
                    h2=habc(jj);
                    if SOL==1
                        Zgabc=[Zgabc,integral(@mutl_impedance_air,1e-4,1)];%% ���������ṹ�����
                    elseif SOL==2
                        Zgabc=[Zgabc,integral(@mutl_impedance_air_double,1e-4,1)];%% ˫�������ṹ�����
                    end
                end                
            end
        end        
    end
 %%%%%%%%����������迹��%%%%%%\ 
 %%%%%%A�������迹����%%%%%%%
    %%%% A�� ��·�迹��%%%%%%%%
    Zg_A(nnn)=Zg_p(1);
    Z33_A(nnn)=Z7(nnn)+Z8(nnn)+Zg_A(nnn);
    %%%A������迹��%%%%%%%%%%
    ZRR_A(nnn)=Z11(nnn)+2*Z12(nnn)+Z22(nnn)+2*Z23(nnn)+Z33_A(nnn);
    ZDD_A(nnn)=Z22(nnn)+2*Z23(nnn)+Z33_A(nnn);
    ZMM_A(nnn)=Z33_A(nnn);
    ZRD_A(nnn)=Z12(nnn)+Z22(nnn)+2*Z23(nnn)+Z33_A(nnn);
    ZDR_A(nnn)=ZRD_A(nnn);
    ZRM_A(nnn)=Z23(nnn)+Z33_A(nnn);
    ZMR_A(nnn)=ZRM_A(nnn);
    ZDM_A(nnn)=Z23(nnn)+Z33_A(nnn);
    ZMD_A(nnn)=Z23(nnn)+Z33_A(nnn);
    Zph_A=[ZRR_A(nnn) ZRD_A(nnn) ZRM_A(nnn);ZDR_A(nnn) ZDD_A(nnn) ZDM_A(nnn);ZMR_A(nnn) ZMD_A(nnn)  ZMM_A(nnn)];
  %%%%%%B�������迹����%%%%%%%
    %%%% B���·�迹��%%%%%%%%
    Zg_B(nnn)=Zg_p(2);
    Z33_B(nnn)=Z7(nnn)+Z8(nnn)+Zg_B(nnn);
    %%%B������迹��%%%%%%%%%%
    ZRR_B(nnn)=Z11(nnn)+2*Z12(nnn)+Z22(nnn)+2*Z23(nnn)+Z33_B(nnn);
    ZDD_B(nnn)=Z22(nnn)+2*Z23(nnn)+Z33_B(nnn);
    ZMM_B(nnn)=Z33_B(nnn);
    ZRD_B(nnn)=Z12(nnn)+Z22(nnn)+2*Z23(nnn)+Z33_B(nnn);
    ZDR_B(nnn)=ZRD_B(nnn);
    ZRM_B(nnn)=Z23(nnn)+Z33_B(nnn);
    ZMR_B(nnn)=ZRM_B(nnn);
    ZDM_B(nnn)=Z23(nnn)+Z33_B(nnn);
    ZMD_B(nnn)=Z23(nnn)+Z33_B(nnn);   
    Zph_B=[ZRR_B(nnn) ZRD_B(nnn) ZRM_B(nnn);ZDR_B(nnn) ZDD_B(nnn) ZDM_B(nnn);ZMR_B(nnn) ZMD_B(nnn)  ZMM_B(nnn)];
  %%%%%%C�������迹����%%%%%%%
    %%%% C���·�迹��%%%%%%%%
    Zg_C(nnn)=Zg_p(3);
    Z33_C(nnn)=Z7(nnn)+Z8(nnn)+Zg_C(nnn);
    %%%B������迹��%%%%%%%%%%
    ZRR_C(nnn)=Z11(nnn)+2*Z12(nnn)+Z22(nnn)+2*Z23(nnn)+Z33_C(nnn);
    ZDD_C(nnn)=Z22(nnn)+2*Z23(nnn)+Z33_C(nnn);
    ZMM_C(nnn)=Z33_C(nnn);
    ZRD_C(nnn)=Z12(nnn)+Z22(nnn)+2*Z23(nnn)+Z33_C(nnn);
    ZDR_C(nnn)=ZRD_C(nnn);
    ZRM_C(nnn)=Z23(nnn)+Z33_C(nnn);
    ZMR_C(nnn)=ZRM_C(nnn);
    ZDM_C(nnn)=Z23(nnn)+Z33_C(nnn);
    ZMD_C(nnn)=Z23(nnn)+Z33_C(nnn); 
    Zph_C=[ZRR_C(nnn) ZRD_C(nnn) ZRM_C(nnn);ZDR_C(nnn) ZDD_C(nnn) ZDM_C(nnn);ZMR_C(nnn) ZMD_C(nnn)  ZMM_C(nnn)];
    %%%%%%%%%���迹%%%%%%%
    ZAB(nnn)=Zgabc(1);
    ZAC(nnn)=Zgabc(2);
    ZBC(nnn)=Zgabc(3);
    
    Zab=[ZAB(nnn) ZAB(nnn) ZAB(nnn);ZAB(nnn) ZAB(nnn) ZAB(nnn);ZAB(nnn) ZAB(nnn) ZAB(nnn)];
    Zac=[ZAC(nnn) ZAC(nnn) ZAC(nnn);ZAC(nnn) ZAC(nnn) ZAC(nnn);ZAC(nnn) ZAC(nnn) ZAC(nnn)];
    Zbc=[ZBC(nnn) ZBC(nnn) ZBC(nnn);ZBC(nnn) ZBC(nnn) ZBC(nnn);ZBC(nnn) ZBC(nnn) ZBC(nnn)];
    Zba=Zab;
    Zca=Zac;
    Zcb=Zbc;    
    %%%%%%%%%%%������µ��迹��%%%%%    
    ZZ(:,:,nnn)=[Zph_A Zab Zac;Zba Zph_B Zbc;Zca Zcb Zph_C];%%��������迹��
    
    %%%%%%%%%%%%%%���㵼�ɵĽ��%%%%%%%
    if Erec==1 %%%%%%%%������ص��µ����%%%%%%%%%%
        Yph(:,:,nnn)=1j*2*pi*f*[CRD,-CRD,0;-CRD,CRD+CDM,-CDM;0,-CDM,CDM+CME];
        %%%%%����ĵ�����%%%%%%
        YY(:,:,nnn)=[Yph(:,:,nnn) zeros(3,3) zeros(3,3);zeros(3,3) Yph(:,:,nnn) zeros(3,3);zeros(3,3) zeros(3,3) Yph(:,:,nnn)];
    else
        for ii=1:Phase
            for jj=1:Phase
                if jj>ii
                    x_mul=abs(xabc(ii)-xabc(jj));
                    h1=habc(ii);
                    h2=habc(jj);
                    D_11=sqrt(x_mul.^2+(h1-h2).^2);
                    D_12=sqrt(x_mul.^2+(h1+h2).^2);
                    Ymul=[Ymul,1j*2*pi*f*2*pi*espo0*log(D_12/D_11)./(log(2*h1/r7).*log(2*h2/r7)-(log(D_12/D_11))^2)];
                end
            end
        end
        Ypha=1j*2*pi*f*[CRD,-CRD,0;-CRD,CRD+CDM,-CDM;0,-CDM,CDM+2*pi*espo0/log(2*habc(1)/r7)];
        Yphb=1j*2*pi*f*[CRD,-CRD,0;-CRD,CRD+CDM,-CDM;0,-CDM,CDM+2*pi*espo0/log(2*habc(2)/r7)];
        Yphc=1j*2*pi*f*[CRD,-CRD,0;-CRD,CRD+CDM,-CDM;0,-CDM,CDM+2*pi*espo0/log(2*habc(3)/r7)];
        Yab=[0 0 0;0 0 0;0 0 Ymul(1)];
        Yac=[0 0 0;0 0 0;0 0 Ymul(2)];
        Ybc=[0 0 0;0 0 0;0 0 Ymul(3)];
        YY(:,:,nnn)=[Ypha Yab Yac;Yab Yphb Ybc;Yac Ybc Yphc];%%������µ�����
    end
end

end
% ZZ_yuan=[Z1.' Z2.' Z3.' Z4.' Z5.' Z6.' Z7.' Z8.' Zg.'];
% ZZ_hl=[Z11.' Z12.' Z21.' Z22.' Z23.' Z33.'];
% ZZ_cab=[ZRR.' ZRD.' ZRM.' ZDR.' ZDD.' ZDM.' ZMR.' ZMD.' ZMM.'];
% [HH_y,LL_y]=size(ZZ_yuan);
% [HH_cab,LL_cab]=size(ZZ_cab);
% for jjj=1:LL_y
%     TEMM=Draw_resistance_figure(ff,ZZ_yuan(:,jjj),1);
% end
% 
% for jjj=1:LL_cab
%     TEMM=Draw_resistance_figure(ff,ZZ_cab(:,jjj),1);
% end



    
    
            
  






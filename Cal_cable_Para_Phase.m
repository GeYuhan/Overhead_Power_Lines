%%%%%%%%%%%%%%%%%计算电缆参数%%%%%%%%
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
%%%定义阻抗的参量
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
%%%%%大地阻抗参数
Zg1=zeros(1,length(ff));
Zg0=zeros(1,length(ff));
Zg2=zeros(1,length(ff));
Zg=zeros(1,length(ff));
%%%%回路阻抗
Z11=zeros(1,length(ff));
Z22=zeros(1,length(ff));
Z12=zeros(1,length(ff));
Z21=zeros(1,length(ff));
Z23=zeros(1,length(ff));
Z32=zeros(1,length(ff));
Z33=zeros(1,length(ff));
%%%%电缆阻抗矩阵
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
f=ff(nnn);%Hz，计算频率
Z2(nnn)=1j*2*pi*f*mu0/2/pi*log(r3/r2);%%电缆导芯与屏蔽层之间
Z5(nnn)=1j*2*pi*f*mu0/2/pi*log(r5/r4);%%屏蔽层与互层之间
Z8(nnn)=1j*2*pi*f*mu0/2/pi*log(r7/r6);%%护层与大地之间

%% 回路阻抗矩阵%%
%%%%%%%%%%%%%%%%%%内外表面阻抗%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1=sqrt(1j*2*pi*f*mu1/rou1);
if r1==0 
    Z1(nnn)=1/(2*pi*r2)*sqrt(1j*2*pi*f*mu1*rou1).*besseli(0,m1*r2)./besseli(1,m1*r2);%%实心导芯表面阻抗
else
    D1=besseli(1,m1*r2).*besselk(1,m1*r1)-besseli(1,m1*r1).*besselk(1,m1*r2);
    Z1(nnn)=rou1.*m1./(2*pi*r2*D1).*(besseli(0,m1*r2).*besselk(1,m1*r1)+besseli(1,m1*r1).*besselk(0,m1*r2));%%空心导芯外表面阻抗
end

m2=sqrt(1j*2*pi*f*mu2/rou2);
D2=besseli(1,m2*r4).*besselk(1,m2*r3)-besseli(1,m2*r3).*besselk(1,m2*r4);
Z3(nnn)=rou2.*m2./(2*pi*r3*D2).*(besseli(0,m2*r3).*besselk(1,m2*r4)+besseli(1,m2*r4).*besselk(0,m2*r3));%%屏蔽层内表面阻抗
Z4(nnn)=rou2.*m2./(2*pi*r4*D2).*(besseli(0,m2*r4).*besselk(1,m2*r3)+besseli(1,m2*r3).*besselk(0,m2*r4));%%屏蔽层外表面阻抗

m3=sqrt(1j*2*pi*f*mu3/rou3);
D3=besseli(1,m3*r6).*besselk(1,m3*r5)-besseli(1,m3*r5).*besselk(1,m3*r6);
Z6(nnn)=rou3.*m3./(2*pi*r5*D3).*(besseli(0,m3*r5).*besselk(1,m3*r6)+besseli(1,m3*r6).*besselk(0,m3*r5));%%护层内表面阻抗
Z7(nnn)=rou3.*m3./(2*pi*r6*D3).*(besseli(0,m3*r6).*besselk(1,m3*r5)+besseli(1,m3*r5).*besselk(0,m3*r6));%%护层外表面阻抗

%%%%%管状导体内外表面之间转移阻抗
Z1m(nnn)=rou2./(2*pi*r3*r4*D2);
Z2m(nnn)=rou3/(2*pi*r5*r6*D3);

%% 大地自返回阻抗%%%%%%%
%%% 近似计算公式
me=sqrt(1j*2*pi*f*mu0/rou4);
%%%%% 考虑大地影响下的大地返回阻抗计算%%%%%
h=ha;
s=1j*2*pi*f;%%%复角频率
if Erec==1  %%%%均匀土壤的情况电缆处于地下时,%% 
    if SOL==1 %单层土壤的情况
        %%%%精确计算方法
        Zg1(nnn)=integral(@self_impedance,1e-6,1);
        %%%%经验公式
        Zg0(nnn)=1j*2*pi*f*mu0/2/pi.*(-log(0.57722*me*r7/2)+0.5-4*me*h/3);      
        %%%%%Pollaczek公式%%%%%
        pp=1./sqrt(s*mu4/rou4);%%%%复镜像深度
        Zg2(nnn)=s*mu0/2/pi*(besselk(0,r7./pp)-besselk(0,2*h./pp)+2*integral(@self_impedance_Pollaczek,1e-6,1));
    elseif SOL==2   %%%双层土壤的情况
        Zg1(nnn)=integral(@self_impedance_ground_double,1e-6,1);  
    end   
else   %%%%%%均匀土壤，电缆置于地上时%%%%%%%%%
    if SOL==1 %单层土壤的情况
        %%%%%%精确计算方法
        Zg1(nnn)=integral(@self_impedance_air,1e-6,1);
        %%%%%%Carson公式
        Zg2(nnn)=2*pi*f/pi*integral(@self_impedance_carson,1e-6,1);
    elseif SOL==2  %%%双层土壤的情况
        Zg1(nnn)=integral(@self_impedance_air_double,1e-6,1);
    end
end

%% 回路阻抗阵%%%%%%%%
Zg(nnn)=Zg1(nnn);
Z11(nnn)=Z1(nnn)+Z2(nnn)+Z3(nnn);
Z22(nnn)=Z4(nnn)+Z5(nnn)+Z6(nnn);
Z33(nnn)=Z7(nnn)+Z8(nnn)+Zg(nnn);
Z12(nnn)=-Z1m(nnn);
Z21(nnn)=-Z1m(nnn);
Z23(nnn)=-Z2m(nnn);
Z32(nnn)=-Z2m(nnn);
% Zh=[Z11,Z12,0;Z21,Z22,Z23;0,Z32,Z33];
%%% 电缆阻抗阵%%%%%%%%%%
ZRR(nnn)=Z11(nnn)+2*Z12(nnn)+Z22(nnn)+2*Z23(nnn)+Z33(nnn);
ZDD(nnn)=Z22(nnn)+2*Z23(nnn)+Z33(nnn);
ZMM(nnn)=Z33(nnn);
ZRD(nnn)=Z12(nnn)+Z22(nnn)+2*Z23(nnn)+Z33(nnn);
ZDR(nnn)=ZRD(nnn);
ZRM(nnn)=Z23(nnn)+Z33(nnn);
ZMR(nnn)=ZRM(nnn);
ZDM(nnn)=Z23(nnn)+Z33(nnn);
ZMD(nnn)=Z23(nnn)+Z33(nnn);
Zph(:,:,nnn)=[ZRR,ZRD,ZRM;ZDR,ZDD,ZDM;ZMR,ZMD,ZMM];%单相电缆阻抗矩阵


%%% 导纳矩阵%%
CRD=2*pi*espo1/log(r3/r2);
CDM=2*pi*espo2/log(r5/r4);
CME=2*pi*espo3/log(r7/r6);
if Erec==2
    CME=2*pi*espo0/log(2*h/r7);
end
Yph(:,:,nnn)=1j*2*pi*f*[CRD,-CRD,0;-CRD,CRD+CDM,-CDM;0,-CDM,CDM+CME];%单相电缆导纳矩阵

 %% 考虑三相导体的情况%%%%%
 Zgabc=[];
 Ymul=[];
 Zg_p=zeros(1,3);
if Phase==3
    %%%%大地互阻抗的求解%%%%%%%%    
    if Erec==1 %%%%%%%%计算埋地电缆的情况%%%%%%%%%%
        for ii=1:Phase
            %%%%计算各相的大地返回阻抗；
            h=habc(ii);
            if SOL==1 %单层土壤的情况
                %%%%精确计算方法
                Zg_p(ii)=integral(@self_impedance,1e-6,1);    
            elseif SOL==2   %%%双层土壤的情况
                Zg_p(ii)=integral(@self_impedance_ground_double,1e-6,1);  
            end   
            %%%%%%%计算三相的互阻抗%%%%%%
            for jj=1:Phase
                if jj>ii
                    x_mul=abs(xabc(ii)-xabc(jj));
                    h1=habc(ii);
                    h2=habc(jj);
                    if SOL==1
                        Zgabc=[Zgabc,integral(@mutl_impedance,1e-4,1)];%% 单层土壤结构的情况
                    elseif SOL==2
                        Zgabc=[Zgabc,integral(@mutl_impedance_ground_double,1e-2,1)];%% 双层土壤结构的情况
                    end
                end                
            end
        end        
    else
        for ii=1:Phase
            %%%%计算各相的大地返回阻抗
            h=habc(ii);
             if SOL==1 %单层土壤的情况
                 %%%%%%精确计算方法
                  Zg_p(ii)=integral(@self_impedance_air,1e-6,1);
             elseif SOL==2  %%%双层土壤的情况
                  Zg_p(ii)=integral(@self_impedance_air_double,1e-6,1);
             end
             %%%%%%%计算三相的互阻抗%%%%%%
            for jj=1:Phase
                if jj>ii
                    x_mul=abs(xabc(ii)-xabc(jj));
                    h1=habc(ii);
                    h2=habc(jj);
                    if SOL==1
                        Zgabc=[Zgabc,integral(@mutl_impedance_air,1e-4,1)];%% 单层土壤结构的情况
                    elseif SOL==2
                        Zgabc=[Zgabc,integral(@mutl_impedance_air_double,1e-4,1)];%% 双层土壤结构的情况
                    end
                end                
            end
        end        
    end
 %%%%%%%%构造三相的阻抗阵%%%%%%\ 
 %%%%%%A相的相关阻抗数据%%%%%%%
    %%%% A相 回路阻抗阵%%%%%%%%
    Zg_A(nnn)=Zg_p(1);
    Z33_A(nnn)=Z7(nnn)+Z8(nnn)+Zg_A(nnn);
    %%%A相电缆阻抗阵%%%%%%%%%%
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
  %%%%%%B相的相关阻抗数据%%%%%%%
    %%%% B相回路阻抗阵%%%%%%%%
    Zg_B(nnn)=Zg_p(2);
    Z33_B(nnn)=Z7(nnn)+Z8(nnn)+Zg_B(nnn);
    %%%B相电缆阻抗阵%%%%%%%%%%
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
  %%%%%%C相的相关阻抗数据%%%%%%%
    %%%% C相回路阻抗阵%%%%%%%%
    Zg_C(nnn)=Zg_p(3);
    Z33_C(nnn)=Z7(nnn)+Z8(nnn)+Zg_C(nnn);
    %%%B相电缆阻抗阵%%%%%%%%%%
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
    %%%%%%%%%互阻抗%%%%%%%
    ZAB(nnn)=Zgabc(1);
    ZAC(nnn)=Zgabc(2);
    ZBC(nnn)=Zgabc(3);
    
    Zab=[ZAB(nnn) ZAB(nnn) ZAB(nnn);ZAB(nnn) ZAB(nnn) ZAB(nnn);ZAB(nnn) ZAB(nnn) ZAB(nnn)];
    Zac=[ZAC(nnn) ZAC(nnn) ZAC(nnn);ZAC(nnn) ZAC(nnn) ZAC(nnn);ZAC(nnn) ZAC(nnn) ZAC(nnn)];
    Zbc=[ZBC(nnn) ZBC(nnn) ZBC(nnn);ZBC(nnn) ZBC(nnn) ZBC(nnn);ZBC(nnn) ZBC(nnn) ZBC(nnn)];
    Zba=Zab;
    Zca=Zac;
    Zcb=Zbc;    
    %%%%%%%%%%%三相电缆的阻抗阵%%%%%    
    ZZ(:,:,nnn)=[Zph_A Zab Zac;Zba Zph_B Zbc;Zca Zcb Zph_C];%%三相电缆阻抗阵
    
    %%%%%%%%%%%%%%计算导纳的结果%%%%%%%
    if Erec==1 %%%%%%%%计算埋地电缆的情况%%%%%%%%%%
        Yph(:,:,nnn)=1j*2*pi*f*[CRD,-CRD,0;-CRD,CRD+CDM,-CDM;0,-CDM,CDM+CME];
        %%%%%三相的导纳阵%%%%%%
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
        YY(:,:,nnn)=[Ypha Yab Yac;Yab Yphb Ybc;Yac Ybc Yphc];%%三相电缆导纳阵
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



    
    
            
  






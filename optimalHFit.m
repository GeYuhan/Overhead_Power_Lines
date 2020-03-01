function [HFit,HmFit,tauOptimal]=optimalHFit(freq,Hm,minusGammaL,H,tolerance,l)
%��ȡ����ʱ�Ż����ԣ�BrentѰ�ţ��õ����������
%   ���������
%       freqƵ��������
%       Hmģ����ά����
%       minusGammaL��-��l������ת�ƺ�����ָ��
%       H������ά����
%       tolerance���������������
%       l��·����
%   ���������
%       HFit,������Ͻṹ��
%       tauOptimal�Ż����������ʱ

%����maxFittingerror�Ѽ��뵽��Ͻ���ṹ����
%����matrixFit�����������Np��������option�ṹ����
%F = matrixFit(frequency,Data,option)
%��ǰoption���������ʽ�����裺
%option.pattern
%option.Np
%option.passivity
%option.weight
%option.realFit
%option.tolerance
%
%F�ǽṹ��
%   F.freq==frequency��Ƶ�������� 1 Ns
%   F.P�Ǽ�����������1 by Np
%   F.C������������άNc Nc Np��ÿһҳ��Ӧһ������
%   F.D�ǳ��������Nc Nc��HΪ��
%   F.T����ʱ������Nk*Np��YcΪ�㣬HΪ���㡣���ڱ������Ƕ�Yc����С��λHmMPS��ϣ��ʾ�Ϊ�㡣���������и�ֵ
%   F.MFE����������maxFittingError������������ֵ����ǲ�����
%   F.RMSE�� ��� �ľ�������������
%
%�Ѷ�tauExtraction����������
%[ tauMPS��tauLarge ] = tauExtraction( omega,Hm,minusGammaL,threshold)
%tauMPS��������С��λϵͳ��ȡ������ʱ�䣬����ʱ���Ż�������
%tauLarge�ǲ�������λ��������Ĵ���ʱ�䣬����tauMPS������ʱ���Ż�������

%Ƶ�ʲ�������
omega=2*pi*freq;%�����
s=omega*1i;%s==j��
%% ��ʱ����ȡ
% tau=tauExtraction(omega,Hm,minusGammaL,0.01);%tau����������0.01�Ǽ���Ƶ��ѡ�����ֵ 
%�ж��Ƿ�Ҫ�ϲ�


%% ģ�����
[tauLeft, tau, tauRight]=tauExtraction(omega,Hm,minusGammaL,0.01,l);%tau����������ѡ�񵱷�ֵ�½���targetErrorʱ��Ƶ��;
%tauLeft=tauMPS;%�Ż���������
%tauRight=tauLarge;%�Ż���������
[tauMin,tauIndex]=min(tau);
theta=25;%�ϲ��ĽǶ���ֵ����λ��
for i=1:length(tau)
    if omega(end)*abs(tau(i)-tauMin)<theta/180*pi & i~=tauIndex
        tau(i)=inf;
%         tauRight(i)=max(tauRight(i),tauRight(tauIndex));
        Hm(tauIndex,tauIndex,:)=Hm(i,i,:)+Hm(tauIndex,tauIndex,:);
    end
end
index=find(tau~=inf);
tau=tau(index);
tauRight=tauRight(index);
tauLeft=tauLeft(index);

tauMiddle=tau;%�Ż������е㡣��Brent's Method��Ҫ������������������ڲ�
for i=1:length(tauMiddle)%��ֹ������ȡ������ʱ���Ż�������
    if tauMiddle(i)<=tauLeft(i)|tauMiddle(i)>=tauRight(i)
        tauMiddle(i)=tauLeft(i)+0.618*(tauRight(i)-tauLeft(i));
    end
end   
Hm=Hm(index,index,:);%����Ӧ������ɾȥ�ˣ������������������С
HmMPS=zeros(size(Hm));

[~, Nc]=size(Hm(:,:,1));
HPoles=cell(Nc,1);%�����������Ԫ���������ڴ���ģ���ļ��㣬�����������
global NpMaxH
NpMax=NpMaxH;%50;%��󼫵���Ŀ%NpMaxH��fittingSet��ȡ
%%�ǵûָ�--------------------------------------------------------@@@@@@@@@@@@@
% NpOptimal=ones(1,Nc)*12;
NpOptimal=zeros(1,Nc);%������Ŀ�ݴ��������ڴ���ÿ��ģ���ļ�����Ŀ�����㿪ʼ�������һ�����б��Զ���һ��ʹ�ü�����Ŀ��1��ʼ���ϵ���
targetError=0.25;%ģ���������������������pscadĬ��0.25��
MFEOptimal=ones(1,Nc)+targetError;%��ʼֵ��Զ����������ʹ�����������
%�������
% tolerance=0.2e-2;%�����������������������0.2%�ο���pscad��
%Brent's Method��������
iterationMax=100;%ԭ��100��
xmin=[];%����Ǽ�Сֵ�����꣬Ҳ����tau
tol=1e-6;%Brent's Method�ľ���
GOLD=0.381966;%�ƽ�ָ��
Nk=50;%�����Ż�������10�����ӣ���0.5��^10>9.7e-04,�ܰ�ģ�򾫶ȷֵ��㹻С
for k=1:Nk
    fprintf('��ʼ��%d���Ż�\n',k)
    %% ģ��
    %����Ԥ������
    %%brentѰ��
    %�������ÿ��ģ��
    for i=1: Nc%ӦΪNc-----------------
        fprintf('��%d��ģ��\n',i)
        if MFEOptimal(i)<=targetError
            fprintf('��%d��ģ�������㾫�ȣ�����\n',i)
            continue   
        else
            if NpOptimal(i)<NpMax
                Np=NpOptimal(i)+1;%����һ���Ż�ʹ�õļ�����Ŀ�����Ͻ���
            else
                fprintf('��%dģ��������Ŀ�Ѵ����ֵ%d\n',i,NpMax);
                continue
            end
        end
        while Np<=NpMax
            %�Ȳ��ı��ʼ���䣬����������к��ٸ�
            %��ʼ���������ֵӦ���Ǿ����ĳ��Ԫ��-----------------------
            %����fa��fb��fc
            %����[a,b,c]��[fa,fb,fc]----------------------------
            %Ҫ���м�ֵ������ֵ��С
%             tau(1)=tauLeft(i);
%             tau(3)=tauMiddle(i);
%             tau(2)=tauRight(i);
%             fprintf('���Ҳ�ֵ%f��targetError%f\n',tau(1)-tau(3),targetError);
            %���㺯��ֵ������ʸ�����
            
            %ȷ������
%             [tauLeft(i), ~, tauRight(i)]=tauExtraction(omega,Hm(i,i,:),minusGammaL(i,i,:),0.001,Np,l);
            for j=1:length(s)%Ƶ��ѭ��
                HmMPS(i,i,j)=Hm(i,i,j)*exp(s(j)*tauMiddle(i));%MPS�� Minimum Phase System,��С��λϵͳ,��tau�滻����u
            end%Ƶ��ѭ������
            %�������ѡ���ǰֻ��pattern��Np���
%             option.pattern=0;
%             option.Np=Np;
%             option.relax=1;
%             option.weight=1;
            %%2019/3/2
            opts.N=Np ;%           %Order of approximation. 
            opts.poletype='linlogcmplx'; %Mix of linearly spaced and logarithmically spaced poles
            opts.stable=1;      %��ת���ȶ�����
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
            HmFit(i)=Error(HmFit(i),HmMPS(i,i,:));%��С��λ�������ģ���������һ����
            MFEMiddle=HmFit(i).MFE;%����������fu
            %Brent's Method Code
            a=min(tauLeft(i),tauRight(i));
            b=max(tauLeft(i),tauRight(i));
            [x,w,v]=deal(tauMiddle(i));
            [fx,fw,fv]=deal(MFEMiddle);
            e=0;
           %% Ѱ�ŵ���
            for iteration=1:iterationMax%Ѱ�ŵ���
                xm=0.5*(a+b);
                tol1=tol*abs(x)+eps;
                tol2=2*tol1;%��������Сֵ��������ĳ���
                %�������㹻С�������Ѱ�ң����ؼ�Сֵ����xmin
                if abs(x-xm)<=(tol2-0.5*(b-a))
                    xmin(i)=x;
                    break
                end
                if abs(e)>tol1%������䳤���ܹ��ƶ���һ�β�����һ�룬�������������ڲ�
                    r=(x-w)*(fx-fv);
                    q=(x-v)*(fx-fw);
                    p=(x-v)*q-(x-w)*r;
                    q=2*(q-r);
                    if q>0
                        p=-p;
                    end
                    q=abs(q);
                    eTemp=e;
                    e=d;
                    %�ϳ�����ʹ�ûƽ�ָ
                    if abs(p)>=abs(0.5*q*eTemp)||p<=q*(a-x)||p>=q*(b-x)
                        if x>=xm
                            e=a-x;%
                            d=GOLD*e;
                        else
                            e=b-x;
                            d=GOLD*e;
                        end
                    else
                        d=p/q;
                        u=x+d;
                        %u must not be evaluated too close to a or b
                        if u-a<tol2||b-u<tol2
                            d=(sign(xm-x)+((xm-x)==0))*tol1;
                        end
                    end%�ϳ�����ʹ�ûƽ�ָ ����
                else%Ѱ���ٶ�̫����ʹ�ûƽ�ָ(���Լ��µ�)
                    if x>=xm
                        e=a-x;
                    else
                        e=b-x;
                    end
                    d=GOLD*e;
                end
                %The u must not be evaluated too close to x
                u=x+(sign(d)+(d==0))*max(abs(d), tol1);
                %���㺯��ֵ������ʸ�����
                for j=1:length(s)%Ƶ��ѭ��
                    HmMPS(i,i,j)=Hm(i,i,j)*exp(s(j)*u);%MPS�� Minimum Phase System,��С��λϵͳ,��tau�滻����u
                end%Ƶ��ѭ������
                %%2019/3/2
                [SER,rmserr,bigHfit]=VFdriver(HmMPS(i,i,:),s,poles,opts);
                HmFitTemp.freq=freq;
                HmFitTemp.P=SER.poles;
                HmFitTemp.C=SER.R;
                HmFitTemp.D=0;
                HmFitTemp.T=0;
                HmFitTemp.MFE=[];
                HmFitTemp.RMSE=rmserr;
                %%2019/3/2
%                 HmFitTemp=matrixFit(freq,HmMPS(i,i,:),option);%HmFit�ǽṹ������,HmFit(i).P���ż���
                HmFitTemp=Error(HmFitTemp,HmMPS(i,i,:));%����2018��12��28�պ󲹵ģ������������ģ��������ǰһֱ����С��λ
                fu=HmFitTemp.MFE;%����������fu
               %% ��������
                if fu <= fx
                    if u >= x
                        a = x; 
                    else
                        b = x; 
                    end
                    v = w; 
                    fv = fw;
                    w = x; 
                    fw = fx;
                    x = u; 
                    fx = fu;
                    HmFit(i)=HmFitTemp;
                else % fu > fx
                    if u < x
                        a = u; 
                    else
                        b = u;
                    end
                    if ( (fu <= fw) || (w == x) )
                        v = w;
                        fv = fw;
                        w = u; 
                        fw = fu;
                    elseif ( (fu <= fv) || (v == x) || (v == w) )
                        v = u; 
                        fv = fu;
                    end
                end
            end%Ѱ�ŵ�������
            if iteration==iterationMax
                fprintf('Brent�ﵽ������\n');
                xmin(i)=x;
            end
           %% �ж�ģ�򾫶��Ƿ�����Ҫ�󣬲�����Ҫ���Ӽ���
            if HmFit(i).MFE>targetError
                if Np<NpMax
                    Np=Np+1;
                    fprintf('������Ŀ��һ��Ϊ%d\n',Np)
                    continue%���Ӽ��㣬������һ��Ѱ��---------------------------------------------
                else
                    fprintf('��%dģ��������Ŀ�Ѵ����ֵ%d\n',i,NpMax);
                    break
                end
            else
                xmin(i)=x;%�����Сֵ�ĺ����꣬���Ż���Ĵ���ʱ���
                fprintf('��%dģ���ڼ�����Ŀ%dѰ�ҳɹ�\n',i,Np)
%                 dataPlot(HmMPS(i,i,:),HmFit(i));%��ͼ
                break
            end
        end%����ѭ������
        tauOptimal(i)=xmin(i);%���Ż���Ĵ���ʱ�������
%         if Np>=NpMax
%             fprintf('��%dģ��������Ŀ�Ѵ����ֵ%d\n',i,NpMax);
%         end
        NpOptimal(i)=min(Np,NpMax);%����ÿ��ģ���ļ�����Ŀ����Ϊ��һ���Ż��Ļ���;��������Ŀ������󼫵���Ŀ����������
        MFEOptimal(i)=HmFit(i).MFE;
        HPoles{i}=[HmFit(i).P];%������ƴ��ΪԪ���������������
    end%�����������ģ������
    %% ����
    HFit=Hresidue(freq,tauOptimal,HPoles,H,0);  
    if HFit.MFE>tolerance
        if sum(NpOptimal)<Nc*NpMax%�������ģ����û�дﵽ������Ŀ����
            targetError=targetError/2;
            continue%targetError���룬������һ�����
        else
            fprintf('��%d���Ż�H������Ŀ�Ѵ����ֵ%d,�ﵽ����%f��Ŀ�꾫��%f\n',k,NpMax*Nc,HFit.MFE,tolerance);
            break;
        end
    else
        fprintf('��%d���Ż�H�ﵽĿ�꾫��%f\n',k,tolerance);
%         dataPlot(Hm(1,1,:),HmFit(1));
        break%�ﵽ���辫�ȣ��Ż���ֹ
    end
end%k=100,�Ż�ѭ������
if k==Nk
    fprintf('�Ż������Ѵ�����%d��\n',k);
end
%% ��Ͻ����������ʾ
global MRPRtolerance
if HFit.MRPR>MRPRtolerance
    HFit.MRPR
    HFit=Hresidue(freq,tauOptimal,HPoles,H,1);
end
fprintf('��Ͻ���\n')
end%��������


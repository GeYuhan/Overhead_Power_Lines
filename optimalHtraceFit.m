function [HFit,HtraceFit,tauOptimal]=optimalHtraceFit(freq,Hm,traceH,minusGammaL,H,tolerance,l)
%% Ƶ�ʲ�������
omega=2*pi*freq;%�����
s=omega*1i;%s==j��
%% ģ�����
[tauLeft, tau, tauRight]=tauExtraction(omega,Hm,minusGammaL,0.001,l);%tau����������ѡ�񵱷�ֵ�½���targetErrorʱ��Ƶ��;
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
tauOptimal=tauMiddle;%��������ֵ
Hm=Hm(index,index,:);%����Ӧ������ɾȥ�ˣ������������������С
HmMPS=zeros(size(Hm));

[~, Nc]=size(Hm(:,:,1));
HPoles=cell(Nc,1);%�����������Ԫ���������ڴ���ģ���ļ��㣬�����������
global NpMaxH
NpMax=NpMaxH;%50;%��󼫵���Ŀ%NpMaxH��fittingSet��ȡ
%%�ǵûָ�--------------------------------------------------------@@@@@@@@@@@@@
% NpOptimal=ones(1,Nc)*12;
NpOptimal=0;%������Ŀ�ݴ��������ڴ���ÿ��ģ���ļ�����Ŀ�����㿪ʼ�������һ�����б��Զ���һ��ʹ�ü�����Ŀ��1��ʼ���ϵ���
targetError=0.25;%ģ�����������������������pscadĬ��0.25��һ�룩
MFEOptimal=1+targetError;%��ʼֵ��Զ����������ʹ�����������
%�������
% tolerance=0.2e-2;%�����������������������0.2%�ο���pscad��
%Brent's Method��������
iterationMax=100;%ԭ��100��
xmin=[];%����Ǽ�Сֵ�����꣬Ҳ����tau
tol=1e-3;%Brent's Method�ľ���
GOLD=0.381966;%�ƽ�ָ��
Nk=40;%�����Ż�������10�����ӣ���0.5��^10>9.7e-04,�ܰ�ģ�򾫶ȷֵ��㹻С
for k=1:Nk
    fprintf('��ʼ��%d���Ż�\n',k)
    %% ģ��
    %����Ԥ������
    %%brentѰ��
    %�������ÿ��ģ��
    
    if MFEOptimal<=targetError
        continue
    else
        Np=NpOptimal+1;%����һ���Ż�ʹ�õļ�����Ŀ�����Ͻ���
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
        %%2019/3/1����ϲ���Ҫ����tauֵ�仯����ԭ����
        %             for j=1:length(s)%Ƶ��ѭ��
        %                 HmMPS(i,i,j)=Hm(i,i,j)*exp(s(j)*tauMiddle(i));%MPS�� Minimum Phase System,��С��λϵͳ,��tau�滻����u
        %             end%Ƶ��ѭ������
        %�������ѡ���ǰֻ��pattern��Np���
        for i=1:Nc%ӦΪNc-----------------
            fprintf('��%d��ģ��\n',i)
            option.pattern=0;
            option.Np=Np;
            option.relax=1;
            option.weight=1;
            HtraceFit=traceFit(freq,tauMiddle,traceH,option);%HmFit�ǽṹ������,HmFit(i).P���ż���
            MFEMiddle=HtraceFit.MFE;%����������fu
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
                tauOptimal(i)=u;
                HtraceFitTemp=traceFit(freq,tauOptimal,traceH,option);%HmFit�ǽṹ������,HmFit(i).P���ż���
                fu=HtraceFitTemp.MFE;%����������fu
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
                    HtraceFit=HtraceFitTemp;
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
        end%�����������ģ������
        %% �ж�ģ�򾫶��Ƿ�����Ҫ�󣬲�����Ҫ���Ӽ���
        if HtraceFit.MFE>targetError
            Np=Np+1;
            fprintf('������Ŀ��һ��Ϊ%d\n',Np)
            continue%���Ӽ��㣬������һ��Ѱ��---------------------------------------------
        else
            xmin(i)=x;%�����Сֵ�ĺ����꣬���Ż���Ĵ���ʱ���
            fprintf('H�ļ��ڼ�����Ŀ%dѰ�ҳɹ�\n',Np)
            %                 dataPlot(HmMPS(i,i,:),HmFit(i));%��ͼ
            break
        end
    end%����ѭ������
    tauOptimal(i)=xmin(i);%���Ż���Ĵ���ʱ�������
    if Np>=NpMax
        fprintf('H�ļ�������Ŀ�Ѵ����ֵ%d\n',NpMax);
    end
    NpOptimal=min(Np,NpMax);%����ÿ��ģ���ļ�����Ŀ����Ϊ��һ���Ż��Ļ���
    MFEOptimal=HtraceFit.MFE;
    for i=1:Nc
        HPoles{i,1}=[HtraceFit.P(HtraceFit.T==tauOptimal(i))];%������ƴ��ΪԪ���������������
    end
    %% ����
    HFit=Hresidue(freq,tauOptimal,HPoles,H,0);
    if HFit.MFE>tolerance
        if NpOptimal<Nc*NpMax%�������ģ����û�дﵽ������Ŀ����
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
if HFit.MRPR>100
    HFit=Hresidue(freq,tauOptimal,HPoles,H,1);
end
fprintf('��Ͻ���\n')
end%��������
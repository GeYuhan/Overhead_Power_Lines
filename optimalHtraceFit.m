function [HFit,HtraceFit,tauOptimal]=optimalHtraceFit(freq,Hm,traceH,minusGammaL,H,tolerance,l)
%% 频率参数生成
omega=2*pi*freq;%计算ω
s=omega*1i;%s==jω
%% 模域参数
[tauLeft, tau, tauRight]=tauExtraction(omega,Hm,minusGammaL,0.001,l);%tau是行向量，选择当幅值下降到targetError时的频点;
%tauLeft=tauMPS;%优化区间下限
%tauRight=tauLarge;%优化区间上限
[tauMin,tauIndex]=min(tau);
theta=25;%合并的角度阈值，单位度
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

tauMiddle=tau;%优化区间中点。因Brent's Method需要三个点进行逆抛物线内插
for i=1:length(tauMiddle)%防止出现提取出的延时在优化区间外
    if tauMiddle(i)<=tauLeft(i)|tauMiddle(i)>=tauRight(i)
        tauMiddle(i)=tauLeft(i)+0.618*(tauRight(i)-tauLeft(i));
    end
end
tauOptimal=tauMiddle;%储存最优值
Hm=Hm(index,index,:);%将对应的行列删去了，矩阵行数和列数会减小
HmMPS=zeros(size(Hm));

[~, Nc]=size(Hm(:,:,1));
HPoles=cell(Nc,1);%声明极点矩阵（元胞），用于储存模量的极点，供相域拟合用
global NpMaxH
NpMax=NpMaxH;%50;%最大极点数目%NpMaxH从fittingSet读取
%%记得恢复--------------------------------------------------------@@@@@@@@@@@@@
% NpOptimal=ones(1,Nc)*12;
NpOptimal=0;%极点数目暂存区，用于储存每个模量的极点数目。从零开始，程序第一次运行便自动加一，使得极点数目从1开始往上递增
targetError=0.25;%模域拟合允许的最大相对拟合误差（是pscad默认0.25的一半）
MFEOptimal=1+targetError;%初始值永远大于允许误差，使程序得以运行
%相域参数
% tolerance=0.2e-2;%相域拟合允许的最大相对拟合误差（0.2%参考自pscad）
%Brent's Method参数设置
iterationMax=100;%原来100次
xmin=[];%这就是极小值的坐标，也就是tau
tol=1e-3;%Brent's Method的精度
GOLD=0.381966;%黄金分割比
Nk=40;%相域优化次数。10次足矣，（0.5）^10>9.7e-04,能把模域精度分的足够小
for k=1:Nk
    fprintf('开始第%d次优化\n',k)
    %% 模域
    %矩阵预先声明
    %%brent寻优
    %逐个计算每个模量
    
    if MFEOptimal<=targetError
        continue
    else
        Np=NpOptimal+1;%在上一次优化使用的极点数目基础上进行
    end
    while Np<=NpMax
        %先不改变初始区间，程序可以运行后再改
        %初始化，后面的值应该是矩阵的某个元素-----------------------
        %计算fa，fb，fc
        %计算[a,b,c]与[fa,fb,fc]----------------------------
        %要求中间值比两边值都小
        %             tau(1)=tauLeft(i);
        %             tau(3)=tauMiddle(i);
        %             tau(2)=tauRight(i);
        %             fprintf('左右差值%f与targetError%f\n',tau(1)-tau(3),targetError);
        %计算函数值，调用矢量拟合
        
        %确定区间
        %             [tauLeft(i), ~, tauRight(i)]=tauExtraction(omega,Hm(i,i,:),minusGammaL(i,i,:),0.001,Np,l);
        %%2019/3/1迹拟合不需要依据tau值变化计算原数据
        %             for j=1:length(s)%频率循环
        %                 HmMPS(i,i,j)=Hm(i,i,j)*exp(s(j)*tauMiddle(i));%MPS即 Minimum Phase System,最小相位系统,将tau替换成了u
        %             end%频率循环结束
        %设置拟合选项（当前只有pattern、Np俩项）
        for i=1:Nc%应为Nc-----------------
            fprintf('第%d个模量\n',i)
            option.pattern=0;
            option.Np=Np;
            option.relax=1;
            option.weight=1;
            HtraceFit=traceFit(freq,tauMiddle,traceH,option);%HmFit是结构体数组,HmFit(i).P存着极点
            MFEMiddle=HtraceFit.MFE;%将均方误差赋给fu
            %Brent's Method Code
            a=min(tauLeft(i),tauRight(i));
            b=max(tauLeft(i),tauRight(i));
            [x,w,v]=deal(tauMiddle(i));
            [fx,fw,fv]=deal(MFEMiddle);
            e=0;
            %% 寻优迭代
            for iteration=1:iterationMax%寻优迭代
                xm=0.5*(a+b);
                tol1=tol*abs(x)+eps;
                tol2=2*tol1;%这就是最后极小值所在区间的长度
                %若区间足够小，则结束寻找，返回极小值坐标xmin
                if abs(x-xm)<=(tol2-0.5*(b-a))
                    xmin(i)=x;
                    break
                end
                if abs(e)>tol1%如果区间长度能够移动上一次步长的一半，则尝试逆抛物线内插
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
                    %较长区间使用黄金分割法
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
                    end%较长区间使用黄金分割法 结束
                else%寻找速度太慢，使用黄金分割法(我自己猜的)
                    if x>=xm
                        e=a-x;
                    else
                        e=b-x;
                    end
                    d=GOLD*e;
                end
                %The u must not be evaluated too close to x
                u=x+(sign(d)+(d==0))*max(abs(d), tol1);
                %计算函数值，调用矢量拟合
                tauOptimal(i)=u;
                HtraceFitTemp=traceFit(freq,tauOptimal,traceH,option);%HmFit是结构体数组,HmFit(i).P存着极点
                fu=HtraceFitTemp.MFE;%将均方误差赋给fu
                %% 更新区间
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
            end%寻优迭代结束
            if iteration==iterationMax
                fprintf('Brent达到最大次数\n');
                xmin(i)=x;
            end
        end%逐个计算摸个模量结束
        %% 判断模域精度是否满足要求，不满足要增加极点
        if HtraceFit.MFE>targetError
            Np=Np+1;
            fprintf('极点数目加一，为%d\n',Np)
            continue%增加极点，继续新一轮寻优---------------------------------------------
        else
            xmin(i)=x;%输出最小值的横坐标，迹优化后的传播时间τ
            fprintf('H的迹在极点数目%d寻找成功\n',Np)
            %                 dataPlot(HmMPS(i,i,:),HmFit(i));%绘图
            break
        end
    end%极点循环结束
    tauOptimal(i)=xmin(i);%将优化后的传播时间存起来
    if Np>=NpMax
        fprintf('H的迹极点数目已达最大值%d\n',NpMax);
    end
    NpOptimal=min(Np,NpMax);%储存每个模量的极点数目，作为下一次优化的基础
    MFEOptimal=HtraceFit.MFE;
    for i=1:Nc
        HPoles{i,1}=[HtraceFit.P(HtraceFit.T==tauOptimal(i))];%将极点拼接为元胞矩阵供相域拟合用
    end
    %% 相域
    HFit=Hresidue(freq,tauOptimal,HPoles,H,0);
    if HFit.MFE>tolerance
        if NpOptimal<Nc*NpMax%如果所有模量都没有达到极点数目极限
            targetError=targetError/2;
            continue%targetError减半，继续新一轮拟合
        else
            fprintf('第%d次优化H极点数目已达最大值%d,达到精度%f，目标精度%f\n',k,NpMax*Nc,HFit.MFE,tolerance);
            break;
        end
    else
        fprintf('第%d次优化H达到目标精度%f\n',k,tolerance);
        %         dataPlot(Hm(1,1,:),HmFit(1));
        break%达到所需精度，优化终止
    end
end%k=100,优化循环结束
if k==Nk
    fprintf('优化次数已达上限%d次\n',k);
end
%% 拟合结束，输出提示
if HFit.MRPR>100
    HFit=Hresidue(freq,tauOptimal,HPoles,H,1);
end
fprintf('拟合结束\n')
end%函数结束
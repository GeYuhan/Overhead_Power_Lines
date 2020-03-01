classdef transmissionLine
    %UNTITLED11 Summary of this class goes here
    %   Detailed explanation goes here
    %   给定T，deltaT，tau,HFit，YcFit，生成传输线类。同时外部也会用这些数据生成控制参数
    %   含有成员 当前时刻k，与外部控制时间同步
    %   给定时刻k，自我更新，可访问iHist0与GGTilde
    %   能够访问vTemp与iTemp
    
    properties
        k;%当前离散时刻
        tauInteger;%tau是时间步长的整数倍
        NpYc;
        NpH;
%         T;%仿真总时间
%         deltaT;%仿真时间步长
%         Niteration;%迭代次数fix(T/deltaT)
%         tauMax;%最大的τ
%         T= 0.5e-3;% Observation time
%         deltaT=0.5e-6;%单位秒,通常是周期的1/20
       %% 状态变量，会变化
        YA;
        YB;
        XA;
        XB;
        X1A;
        X1B;
        X2A;
        X2B;
        X3A;
        X3B;
        b2A;
        b2B;
        %% 插值计算参数,会变化
        vTemp;%waibu
        vATemp;
        vBTemp;
        iTemp;%waibu
        iATemp;
        iBTemp;
        ifarA;
        ifarB;
        %对（t-τ）进行线性插值的系数
        epsilon;
        alpha2;
        alpha1;         
        %% 状态变量参数生成
        ai;
        au;
        G;
        GiTilde
        aki;
        aku;
        lambda1;
        judge1
        judge0;
        RkiTilde;
        lambda0;
        %%等效电路
        iHist0;%外部
        GGTilde%外部
    end
    
    methods
        function obj=transmissionLine(T,deltaT,tau,HFit,YcFit)%初始化构造函数
            %%离散时间
            obj.tauInteger=fix(HFit.T./deltaT);%向量
            Niteration=fix(T/deltaT);%迭代次数
            %         tauInteger= fix(tau./deltaT);    % 得到了tau对Δt的倍数。tau是一个行向量，存储着三个延时，%下面直接算出了此值，构造插值系数
            tauMax = fix(max(tau)./deltaT);    %最大的τ
            %%矩阵大小参数
            [Nc,~]=size(HFit.C(:,:,1));
            obj.NpYc = length(YcFit.P); % Number of poles of Yc
            obj.NpH= length(HFit.P); % H的极点数目
            NpYc=obj.NpYc;
            NpH=obj.NpH;
            %%状态变量初始化
            %Yc
            obj.YA = zeros(Nc,NpYc); % State variables
            obj.YB  = zeros(Nc,NpYc);  % State variables
            %H
            obj.XA  = zeros(Nc,NpH);  %ihhA1状态变量
            obj.XB  = zeros(Nc,NpH);  %ihhB1状态变量
            obj.X2A  = zeros(Nc,NpH);  %ihhA2状态变量
            obj.X2B  = zeros(Nc,NpH);  %ihhB2状态变量
            obj.X3A  = zeros(Nc,NpH);  %ihhA3状态变量
            obj.X3B  = zeros(Nc,NpH);  %ihhB3状态变量
            obj.b2A  = zeros(Nc,NpH);
            obj.b2B  = zeros(Nc,NpH);
            %%插值用的参数
            %建立电压、电流记忆区，存储整数倍时刻前的值，那么长度为tauMax+2
            obj.vTemp = zeros(2*Nc,tauMax+2);  % 储存τ时刻之前的电压值
            obj.vATemp=zeros(Nc,tauMax+2);
            obj.vBTemp=zeros(Nc,tauMax+2);
            obj.iTemp = zeros(2*Nc,tauMax+2);  % 储存τ时刻之前的电压值
            obj.iATemp=zeros(Nc,tauMax+2);
            obj.iBTemp=zeros(Nc,tauMax+2);
            obj.ifarA=zeros(Nc,tauMax+3);
            obj.ifarB=zeros(Nc,tauMax+3);
            %%线性插值系数
            % obj.epsilon =tau - tauInteger*deltaT;  % Time for the interpolatI0n，ε，除以Δt就是论文中的ξ，向量，三个值
            obj.epsilon =HFit.T - fix(HFit.T/deltaT)*deltaT;%NpH长度向量
            obj.alpha2 = obj.epsilon./deltaT;%
            obj.alpha1 = 1-obj.alpha2; %
            %%状态变量所用参数
            %note：YC —— Ai、GiTilde(即Gi~)，obj.G
            obj.ai = (1+(deltaT/2)*YcFit.P)./(1-(deltaT/2)*YcFit.P);%就是ppt中的ai
            obj.au = ((deltaT/2)./(1-(deltaT/2)*YcFit.P)); %乘上留数就是ppt中的Gi~
            
            obj.G = zeros(Nc); %等效电路中的导纳
            for nm = 1:NpYc
                obj.GiTilde(:,:,nm) = YcFit.C(:,:,nm)*obj.au(nm);   %就是ppt中的Gi^
                obj.G  = obj.G + obj.GiTilde(:,:,nm); %
            end
            obj.G  = obj.G + YcFit.D;   % Admitance of the Ish，就是最终等效电路里的  obj.G
            % Constants for the states obj.YA and obj.YB ，就是H相关
            
            %note：H —— Aki、obj.RkiTilde(即Rki~)
            obj.aki = (1+(deltaT/2)*HFit.P)./(1-(deltaT/2)*HFit.P);   %就是ppt中的a_k,i
            obj.aku = (((deltaT/2))./(1-(deltaT/2)*HFit.P));    %乘上留数就是ppt中的R_k,i~
            obj.lambda1=zeros(Nc);%先用来记录RkiTilde的和
            obj.judge1=fix(HFit.T./deltaT)==0;%1xNpH行向量，Δt>τ，为1；Δt<τ,为0
            obj.judge0=abs(obj.judge1-1);%1xNpH行向量，Δt>τ，为0；Δt<τ,为1
            for k = 1: NpH
                obj.RkiTilde(:,:,k)=  HFit.C(:,:,k).*obj.aku(k);  %就是ppt中的R_k,i~
                obj.lambda1=obj.lambda1+ obj.judge1(k)*obj.RkiTilde(:,:,k)*obj.alpha1(k);%若Δt>τ，为1
            end
            %插值与迭代卷积后得到的参数
            obj.lambda0=obj.G*obj.lambda1;
        end
        function obj=iteration(obj)%状态迭代函数
            %Yc状态变量，一共3个
            %Y1A对应
            vn=obj.vTemp(:,1);
            vAn=vn(1:end/2);
            vBn=vn(end/2+1:end);
            %%这里增加in
            in=obj.iTemp(:,1);
            iAn=in(1:end/2);
            iBn=in(end/2+1:end);
            %%2019/3/23
%             correctionA=0;
%             correctionB=0;
            for m = 1:obj.NpYc
                obj.YA(:,m) = obj.ai(m)*obj.YA(:,m) + obj.GiTilde(:,:,m)*(obj.ai(m)+1)*vAn;%就是Gi^,vn即上一时刻的电压
                obj.YB(:,m) = obj.ai(m)*obj.YB(:,m) + obj.GiTilde(:,:,m)*(obj.ai(m)+1)*vBn;
                %%2019/3/23
%                 correctionA=correctionA+obj.GiTilde(:,:,m)*obj.ai(m)*vAn;
%                 correctionB=correctionB+obj.GiTilde(:,:,m)*obj.ai(m)*vBn;
            end
            %这里ihyTemp改为Ifar
            %先更新一下，第一列永远保存最新值
            obj.ifarB(:,2:end)=obj.ifarB(:,1:end-1);
            obj.ifarA(:,2:end)=obj.ifarA(:,1:end-1);
            iyA=sum(obj.YA,2);
            iyB=sum(obj.YB,2);
            obj.ifarB(:,1)=iAn+obj.G*vAn+sum(obj.YA,2);%用于B端计算
            obj.ifarA(:,1)=iBn+obj.G*vBn+sum(obj.YB,2);%用于A端计算
%             ihyA=sum(obj.YA,2)-correctionA;
%             ihyB=sum(obj.YB,2)-correctionB;           
            %H状态变量，一共4个
            obj.b2A=0;
            obj.b2B=0;
            obj.vATemp=obj.vTemp(1:end/2,:);
            obj.vBTemp=obj.vTemp(end/2+1:end,:);
            obj.iATemp=obj.iTemp(1:end/2,:);
            obj.iBTemp=obj.iTemp(end/2+1:end,:);
            for m = 1:obj.NpH
                %K1就是a_k,i, K2就是R~_k,i
                %这里先留个坑，就是将H拟合改为区分τ，这里就使用结构体的数据
                %但这样做的话，就不能做到HFit与YcFit的统一。后面没空就不改了
                %ihhA1,ihhB1
                
                %修改之后
                obj.XA(:,m)=obj.aki(m)*obj.XA(:,m) + obj.RkiTilde(:,:,m)*( obj.alpha1(m)*obj.ifarA(:,obj.tauInteger(m))+obj.ifarA(:,obj.tauInteger(m)+1)+obj.alpha2(m)*obj.ifarA(:,obj.tauInteger(m)+2) );
                obj.XB(:,m)=obj.aki(m)*obj.XB(:,m) + obj.RkiTilde(:,:,m)*( obj.alpha1(m)*obj.ifarB(:,obj.tauInteger(m))+obj.ifarB(:,obj.tauInteger(m)+1)+obj.alpha2(m)*obj.ifarB(:,obj.tauInteger(m)+2) );
%                 
%                 obj.X1A(:,m)=obj.aki(m)*obj.X1A(:,m) + obj.RkiTilde(:,:,m)*(obj.aki(m)*obj.alpha1(m)+1)*obj.vBTemp(:,obj.tauInteger(m)+1)+obj.RkiTilde(:,:,m)*obj.alpha2(m)*obj.vBTemp(:,obj.tauInteger(m)+2);%fix(HFit.T(m)./deltaT)即相应τ的整数部分
%                 obj.X1B(:,m)=obj.aki(m)*obj.X1B(:,m) + obj.RkiTilde(:,:,m)*(obj.aki(m)*obj.alpha1(m)+1)*obj.vATemp(:,obj.tauInteger(m)+1)+obj.RkiTilde(:,:,m)*obj.alpha2(m)*obj.vATemp(:,obj.tauInteger(m)+2);%fix(HFit.T(m)./deltaT)即相应τ的整数部分
%                 %ihhA2,ihhB2
%                 obj.X2A(:,m)=obj.aki(m)*obj.X2A(:,m) + obj.RkiTilde(:,:,m)*(obj.aki(m)*obj.alpha1(m)+1)*obj.iBTemp(:,obj.tauInteger(m)+1)+obj.RkiTilde(:,:,m)*obj.alpha2(m)*obj.iBTemp(:,obj.tauInteger(m)+2);%fix(HFit.T(m)./deltaT)即相应τ的整数部分
%                 obj.X2B(:,m)=obj.aki(m)*obj.X2B(:,m) + obj.RkiTilde(:,:,m)*(obj.aki(m)*obj.alpha1(m)+1)*obj.iATemp(:,obj.tauInteger(m)+1)+obj.RkiTilde(:,:,m)*obj.alpha2(m)*obj.iATemp(:,obj.tauInteger(m)+2);%fix(HFit.T(m)./deltaT)即相应τ的整数部分
%                 %ihhA3,ihhB3
%                 obj.X3A(:,m)=obj.aki(m)*obj.X3A(:,m) + obj.RkiTilde(:,:,m)*obj.alpha1(m)*obj.ihyBTemp(:,obj.tauInteger(m)+1)+ obj.RkiTilde(:,:,m)*obj.ihyBTemp(:,obj.tauInteger(m)+2)+ obj.RkiTilde(:,:,m)*obj.alpha2(m)*obj.ihyBTemp(:,obj.tauInteger(m)+3);
%                 obj.X3B(:,m)=obj.aki(m)*obj.X3B(:,m) + obj.RkiTilde(:,:,m)*obj.alpha1(m)*obj.ihyATemp(:,obj.tauInteger(m)+1)+ obj.RkiTilde(:,:,m)*obj.ihyATemp(:,obj.tauInteger(m)+2)+ obj.RkiTilde(:,:,m)*obj.alpha2(m)*obj.ihyATemp(:,obj.tauInteger(m)+3);
%                 %ihhA3,ihhB3(这是因为考虑到ihy现在的值（m）是可以读取的，只不过要在语法上适应τ<Δt)
% %                 obj.X3A(:,m)=obj.aki(m)*obj.X3A(:,m) + obj.RkiTilde(:,:,m)*obj.alpha1(m)*obj.ihyBTemp(:,obj.tauInteger(m))+ obj.RkiTilde(:,:,m)*obj.ihyBTemp(:,obj.tauInteger(m)+1)+ obj.RkiTilde(:,:,m)*obj.alpha2(m)*obj.ihyBTemp(:,obj.tauInteger(m)+2);
% %                 obj.X3B(:,m)=obj.aki(m)*obj.X3B(:,m) + obj.RkiTilde(:,:,m)*obj.alpha1(m)*obj.ihyATemp(:,obj.tauInteger(m))+ obj.RkiTilde(:,:,m)*obj.ihyATemp(:,obj.tauInteger(m)+1)+ obj.RkiTilde(:,:,m)*obj.alpha2(m)*obj.ihyATemp(:,obj.tauInteger(m)+2);
%                 %λλ02,λλ12,
%                 obj.b2A=obj.b2A+ obj.judge0(m)*obj.RkiTilde(:,:,m)*obj.alpha1(m)*( obj.iBTemp(:,obj.tauInteger(m)+obj.judge1(m))+obj.G*obj.vBTemp(:,obj.tauInteger(m)+obj.judge1(m)) );%若Δt>τ，为0,v(n-k)
%                 obj.b2B=obj.b2B+ obj.judge0(m)*obj.RkiTilde(:,:,m)*obj.alpha1(m)*( obj.iATemp(:,obj.tauInteger(m)+obj.judge1(m))+obj.G*obj.vATemp(:,obj.tauInteger(m)+obj.judge1(m)) );
%                 %%2019/3/23
%                 %%A端的两个修正项
% %                 correctionA= obj.RkiTilde(:,:,m)*obj.alpha1(m)*obj.aki(m)*( obj.iBTemp(:,obj.tauInteger(m)+1)+obj.G*obj.vBTemp(:,obj.tauInteger(m)+1) );
% %                 correctionB= obj.RkiTilde(:,:,m)*obj.alpha1(m)*obj.aki(m)*( obj.iATemp(:,obj.tauInteger(m)+1)+obj.G*obj.vATemp(:,obj.tauInteger(m)+1) );
%                 %%我觉的judge1可以去掉，因为一旦judge0等于0，后面的无意义；不可以，去掉运行会报错，除非改成判断，judge0等于1才执行
%                 obj.b2A=obj.b2A;
%                 obj.b2B=obj.b2B;
% %                 obj.b2A=obj.b2A-correctionA;
% %                 obj.b2B=obj.b2B-correctionB;
            end
            %更新记忆区
            obj.vTemp(:,2:end)=obj.vTemp(:,1:end-1);
            obj.iTemp(:,2:end)=obj.iTemp(:,1:end-1);
%             obj.ihyATemp(:,2:end)=obj.ihyATemp(:,1:end-1);
%             obj.ihyBTemp(:,2:end)=obj.ihyBTemp(:,1:end-1);
%             %状态空间输出
%             ihhA1=sum(obj.X1A,2);
%             ihhB1=sum(obj.X1B,2);
%             ihhA2=sum(obj.X2A,2);
%             ihhB2=sum(obj.X2B,2);
%             ihhA3=sum(obj.X3A,2);
%             ihhB3=sum(obj.X3B,2);
%             ihhA=obj.G*ihhA1+ihhA2+ihhA3;
%             ihhB=obj.G*ihhB1+ihhB2+ihhB3;
            %历史电流源第一次计算
            iHistA=sum(obj.XA,2)-iyA;%初步得到历史电流源
            iHistB=sum(obj.XB,2)-iyB;
            %% 矩阵组装
            obj.GGTilde=[obj.G,zeros(size(obj.G));zeros(size(obj.G)),obj.G];
%             doubleLambda0=[zeros(size(obj.lambda0)),obj.lambda0;obj.lambda0,zeros(size(obj.lambda0))];
%             doubleLambda1=[obj.lambda1,zeros(size(obj.lambda1));zeros(size(obj.lambda1)),obj.lambda1];
%             A=eye(size(doubleLambda1))+doubleLambda1;%电流的系数
%             obj.GGTilde=A\(GG-doubleLambda0);%等效电路的电导
%             b2=[obj.b2A;obj.b2B];
%             obj.iHist0=A\([iHistA;iHistB]+b2);
            obj.iHist0=[iHistA;iHistB];
        end
        function obj=setvTemp(obj,vn)
            obj.vTemp(:,1)=vn;
        end
        function obj=setiTemp(obj,in)
            obj.iTemp(:,1)=in;
        end
        function out=getiHist0(obj)
            out=obj.iHist0;
        end
        function out=getGGTilde(obj)
            out=obj.GGTilde;
        end
    end
    
end


classdef transmissionLine
    %UNTITLED11 Summary of this class goes here
    %   Detailed explanation goes here
    %   ����T��deltaT��tau,HFit��YcFit�����ɴ������ࡣͬʱ�ⲿҲ������Щ�������ɿ��Ʋ���
    %   ���г�Ա ��ǰʱ��k�����ⲿ����ʱ��ͬ��
    %   ����ʱ��k�����Ҹ��£��ɷ���iHist0��GGTilde
    %   �ܹ�����vTemp��iTemp
    
    properties
        k;%��ǰ��ɢʱ��
        tauInteger;%tau��ʱ�䲽����������
        NpYc;
        NpH;
%         T;%������ʱ��
%         deltaT;%����ʱ�䲽��
%         Niteration;%��������fix(T/deltaT)
%         tauMax;%���Ħ�
%         T= 0.5e-3;% Observation time
%         deltaT=0.5e-6;%��λ��,ͨ�������ڵ�1/20
       %% ״̬��������仯
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
        %% ��ֵ�������,��仯
        vTemp;%waibu
        vATemp;
        vBTemp;
        iTemp;%waibu
        iATemp;
        iBTemp;
        ifarA;
        ifarB;
        %�ԣ�t-�ӣ��������Բ�ֵ��ϵ��
        epsilon;
        alpha2;
        alpha1;         
        %% ״̬������������
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
        %%��Ч��·
        iHist0;%�ⲿ
        GGTilde%�ⲿ
    end
    
    methods
        function obj=transmissionLine(T,deltaT,tau,HFit,YcFit)%��ʼ�����캯��
            %%��ɢʱ��
            obj.tauInteger=fix(HFit.T./deltaT);%����
            Niteration=fix(T/deltaT);%��������
            %         tauInteger= fix(tau./deltaT);    % �õ���tau�Ԧ�t�ı�����tau��һ�����������洢��������ʱ��%����ֱ������˴�ֵ�������ֵϵ��
            tauMax = fix(max(tau)./deltaT);    %���Ħ�
            %%�����С����
            [Nc,~]=size(HFit.C(:,:,1));
            obj.NpYc = length(YcFit.P); % Number of poles of Yc
            obj.NpH= length(HFit.P); % H�ļ�����Ŀ
            NpYc=obj.NpYc;
            NpH=obj.NpH;
            %%״̬������ʼ��
            %Yc
            obj.YA = zeros(Nc,NpYc); % State variables
            obj.YB  = zeros(Nc,NpYc);  % State variables
            %H
            obj.XA  = zeros(Nc,NpH);  %ihhA1״̬����
            obj.XB  = zeros(Nc,NpH);  %ihhB1״̬����
            obj.X2A  = zeros(Nc,NpH);  %ihhA2״̬����
            obj.X2B  = zeros(Nc,NpH);  %ihhB2״̬����
            obj.X3A  = zeros(Nc,NpH);  %ihhA3״̬����
            obj.X3B  = zeros(Nc,NpH);  %ihhB3״̬����
            obj.b2A  = zeros(Nc,NpH);
            obj.b2B  = zeros(Nc,NpH);
            %%��ֵ�õĲ���
            %������ѹ���������������洢������ʱ��ǰ��ֵ����ô����ΪtauMax+2
            obj.vTemp = zeros(2*Nc,tauMax+2);  % �����ʱ��֮ǰ�ĵ�ѹֵ
            obj.vATemp=zeros(Nc,tauMax+2);
            obj.vBTemp=zeros(Nc,tauMax+2);
            obj.iTemp = zeros(2*Nc,tauMax+2);  % �����ʱ��֮ǰ�ĵ�ѹֵ
            obj.iATemp=zeros(Nc,tauMax+2);
            obj.iBTemp=zeros(Nc,tauMax+2);
            obj.ifarA=zeros(Nc,tauMax+3);
            obj.ifarB=zeros(Nc,tauMax+3);
            %%���Բ�ֵϵ��
            % obj.epsilon =tau - tauInteger*deltaT;  % Time for the interpolatI0n���ţ����Ԧ�t���������еĦΣ�����������ֵ
            obj.epsilon =HFit.T - fix(HFit.T/deltaT)*deltaT;%NpH��������
            obj.alpha2 = obj.epsilon./deltaT;%
            obj.alpha1 = 1-obj.alpha2; %
            %%״̬�������ò���
            %note��YC ���� Ai��GiTilde(��Gi~)��obj.G
            obj.ai = (1+(deltaT/2)*YcFit.P)./(1-(deltaT/2)*YcFit.P);%����ppt�е�ai
            obj.au = ((deltaT/2)./(1-(deltaT/2)*YcFit.P)); %������������ppt�е�Gi~
            
            obj.G = zeros(Nc); %��Ч��·�еĵ���
            for nm = 1:NpYc
                obj.GiTilde(:,:,nm) = YcFit.C(:,:,nm)*obj.au(nm);   %����ppt�е�Gi^
                obj.G  = obj.G + obj.GiTilde(:,:,nm); %
            end
            obj.G  = obj.G + YcFit.D;   % Admitance of the Ish���������յ�Ч��·���  obj.G
            % Constants for the states obj.YA and obj.YB ������H���
            
            %note��H ���� Aki��obj.RkiTilde(��Rki~)
            obj.aki = (1+(deltaT/2)*HFit.P)./(1-(deltaT/2)*HFit.P);   %����ppt�е�a_k,i
            obj.aku = (((deltaT/2))./(1-(deltaT/2)*HFit.P));    %������������ppt�е�R_k,i~
            obj.lambda1=zeros(Nc);%��������¼RkiTilde�ĺ�
            obj.judge1=fix(HFit.T./deltaT)==0;%1xNpH����������t>�ӣ�Ϊ1����t<��,Ϊ0
            obj.judge0=abs(obj.judge1-1);%1xNpH����������t>�ӣ�Ϊ0����t<��,Ϊ1
            for k = 1: NpH
                obj.RkiTilde(:,:,k)=  HFit.C(:,:,k).*obj.aku(k);  %����ppt�е�R_k,i~
                obj.lambda1=obj.lambda1+ obj.judge1(k)*obj.RkiTilde(:,:,k)*obj.alpha1(k);%����t>�ӣ�Ϊ1
            end
            %��ֵ����������õ��Ĳ���
            obj.lambda0=obj.G*obj.lambda1;
        end
        function obj=iteration(obj)%״̬��������
            %Yc״̬������һ��3��
            %Y1A��Ӧ
            vn=obj.vTemp(:,1);
            vAn=vn(1:end/2);
            vBn=vn(end/2+1:end);
            %%��������in
            in=obj.iTemp(:,1);
            iAn=in(1:end/2);
            iBn=in(end/2+1:end);
            %%2019/3/23
%             correctionA=0;
%             correctionB=0;
            for m = 1:obj.NpYc
                obj.YA(:,m) = obj.ai(m)*obj.YA(:,m) + obj.GiTilde(:,:,m)*(obj.ai(m)+1)*vAn;%����Gi^,vn����һʱ�̵ĵ�ѹ
                obj.YB(:,m) = obj.ai(m)*obj.YB(:,m) + obj.GiTilde(:,:,m)*(obj.ai(m)+1)*vBn;
                %%2019/3/23
%                 correctionA=correctionA+obj.GiTilde(:,:,m)*obj.ai(m)*vAn;
%                 correctionB=correctionB+obj.GiTilde(:,:,m)*obj.ai(m)*vBn;
            end
            %����ihyTemp��ΪIfar
            %�ȸ���һ�£���һ����Զ��������ֵ
            obj.ifarB(:,2:end)=obj.ifarB(:,1:end-1);
            obj.ifarA(:,2:end)=obj.ifarA(:,1:end-1);
            iyA=sum(obj.YA,2);
            iyB=sum(obj.YB,2);
            obj.ifarB(:,1)=iAn+obj.G*vAn+sum(obj.YA,2);%����B�˼���
            obj.ifarA(:,1)=iBn+obj.G*vBn+sum(obj.YB,2);%����A�˼���
%             ihyA=sum(obj.YA,2)-correctionA;
%             ihyB=sum(obj.YB,2)-correctionB;           
            %H״̬������һ��4��
            obj.b2A=0;
            obj.b2B=0;
            obj.vATemp=obj.vTemp(1:end/2,:);
            obj.vBTemp=obj.vTemp(end/2+1:end,:);
            obj.iATemp=obj.iTemp(1:end/2,:);
            obj.iBTemp=obj.iTemp(end/2+1:end,:);
            for m = 1:obj.NpH
                %K1����a_k,i, K2����R~_k,i
                %�����������ӣ����ǽ�H��ϸ�Ϊ���֦ӣ������ʹ�ýṹ�������
                %���������Ļ����Ͳ�������HFit��YcFit��ͳһ������û�վͲ�����
                %ihhA1,ihhB1
                
                %�޸�֮��
                obj.XA(:,m)=obj.aki(m)*obj.XA(:,m) + obj.RkiTilde(:,:,m)*( obj.alpha1(m)*obj.ifarA(:,obj.tauInteger(m))+obj.ifarA(:,obj.tauInteger(m)+1)+obj.alpha2(m)*obj.ifarA(:,obj.tauInteger(m)+2) );
                obj.XB(:,m)=obj.aki(m)*obj.XB(:,m) + obj.RkiTilde(:,:,m)*( obj.alpha1(m)*obj.ifarB(:,obj.tauInteger(m))+obj.ifarB(:,obj.tauInteger(m)+1)+obj.alpha2(m)*obj.ifarB(:,obj.tauInteger(m)+2) );
%                 
%                 obj.X1A(:,m)=obj.aki(m)*obj.X1A(:,m) + obj.RkiTilde(:,:,m)*(obj.aki(m)*obj.alpha1(m)+1)*obj.vBTemp(:,obj.tauInteger(m)+1)+obj.RkiTilde(:,:,m)*obj.alpha2(m)*obj.vBTemp(:,obj.tauInteger(m)+2);%fix(HFit.T(m)./deltaT)����Ӧ�ӵ���������
%                 obj.X1B(:,m)=obj.aki(m)*obj.X1B(:,m) + obj.RkiTilde(:,:,m)*(obj.aki(m)*obj.alpha1(m)+1)*obj.vATemp(:,obj.tauInteger(m)+1)+obj.RkiTilde(:,:,m)*obj.alpha2(m)*obj.vATemp(:,obj.tauInteger(m)+2);%fix(HFit.T(m)./deltaT)����Ӧ�ӵ���������
%                 %ihhA2,ihhB2
%                 obj.X2A(:,m)=obj.aki(m)*obj.X2A(:,m) + obj.RkiTilde(:,:,m)*(obj.aki(m)*obj.alpha1(m)+1)*obj.iBTemp(:,obj.tauInteger(m)+1)+obj.RkiTilde(:,:,m)*obj.alpha2(m)*obj.iBTemp(:,obj.tauInteger(m)+2);%fix(HFit.T(m)./deltaT)����Ӧ�ӵ���������
%                 obj.X2B(:,m)=obj.aki(m)*obj.X2B(:,m) + obj.RkiTilde(:,:,m)*(obj.aki(m)*obj.alpha1(m)+1)*obj.iATemp(:,obj.tauInteger(m)+1)+obj.RkiTilde(:,:,m)*obj.alpha2(m)*obj.iATemp(:,obj.tauInteger(m)+2);%fix(HFit.T(m)./deltaT)����Ӧ�ӵ���������
%                 %ihhA3,ihhB3
%                 obj.X3A(:,m)=obj.aki(m)*obj.X3A(:,m) + obj.RkiTilde(:,:,m)*obj.alpha1(m)*obj.ihyBTemp(:,obj.tauInteger(m)+1)+ obj.RkiTilde(:,:,m)*obj.ihyBTemp(:,obj.tauInteger(m)+2)+ obj.RkiTilde(:,:,m)*obj.alpha2(m)*obj.ihyBTemp(:,obj.tauInteger(m)+3);
%                 obj.X3B(:,m)=obj.aki(m)*obj.X3B(:,m) + obj.RkiTilde(:,:,m)*obj.alpha1(m)*obj.ihyATemp(:,obj.tauInteger(m)+1)+ obj.RkiTilde(:,:,m)*obj.ihyATemp(:,obj.tauInteger(m)+2)+ obj.RkiTilde(:,:,m)*obj.alpha2(m)*obj.ihyATemp(:,obj.tauInteger(m)+3);
%                 %ihhA3,ihhB3(������Ϊ���ǵ�ihy���ڵ�ֵ��m���ǿ��Զ�ȡ�ģ�ֻ����Ҫ���﷨����Ӧ��<��t)
% %                 obj.X3A(:,m)=obj.aki(m)*obj.X3A(:,m) + obj.RkiTilde(:,:,m)*obj.alpha1(m)*obj.ihyBTemp(:,obj.tauInteger(m))+ obj.RkiTilde(:,:,m)*obj.ihyBTemp(:,obj.tauInteger(m)+1)+ obj.RkiTilde(:,:,m)*obj.alpha2(m)*obj.ihyBTemp(:,obj.tauInteger(m)+2);
% %                 obj.X3B(:,m)=obj.aki(m)*obj.X3B(:,m) + obj.RkiTilde(:,:,m)*obj.alpha1(m)*obj.ihyATemp(:,obj.tauInteger(m))+ obj.RkiTilde(:,:,m)*obj.ihyATemp(:,obj.tauInteger(m)+1)+ obj.RkiTilde(:,:,m)*obj.alpha2(m)*obj.ihyATemp(:,obj.tauInteger(m)+2);
%                 %�˦�02,�˦�12,
%                 obj.b2A=obj.b2A+ obj.judge0(m)*obj.RkiTilde(:,:,m)*obj.alpha1(m)*( obj.iBTemp(:,obj.tauInteger(m)+obj.judge1(m))+obj.G*obj.vBTemp(:,obj.tauInteger(m)+obj.judge1(m)) );%����t>�ӣ�Ϊ0,v(n-k)
%                 obj.b2B=obj.b2B+ obj.judge0(m)*obj.RkiTilde(:,:,m)*obj.alpha1(m)*( obj.iATemp(:,obj.tauInteger(m)+obj.judge1(m))+obj.G*obj.vATemp(:,obj.tauInteger(m)+obj.judge1(m)) );
%                 %%2019/3/23
%                 %%A�˵�����������
% %                 correctionA= obj.RkiTilde(:,:,m)*obj.alpha1(m)*obj.aki(m)*( obj.iBTemp(:,obj.tauInteger(m)+1)+obj.G*obj.vBTemp(:,obj.tauInteger(m)+1) );
% %                 correctionB= obj.RkiTilde(:,:,m)*obj.alpha1(m)*obj.aki(m)*( obj.iATemp(:,obj.tauInteger(m)+1)+obj.G*obj.vATemp(:,obj.tauInteger(m)+1) );
%                 %%�Ҿ���judge1����ȥ������Ϊһ��judge0����0������������壻�����ԣ�ȥ�����лᱨ�����Ǹĳ��жϣ�judge0����1��ִ��
%                 obj.b2A=obj.b2A;
%                 obj.b2B=obj.b2B;
% %                 obj.b2A=obj.b2A-correctionA;
% %                 obj.b2B=obj.b2B-correctionB;
            end
            %���¼�����
            obj.vTemp(:,2:end)=obj.vTemp(:,1:end-1);
            obj.iTemp(:,2:end)=obj.iTemp(:,1:end-1);
%             obj.ihyATemp(:,2:end)=obj.ihyATemp(:,1:end-1);
%             obj.ihyBTemp(:,2:end)=obj.ihyBTemp(:,1:end-1);
%             %״̬�ռ����
%             ihhA1=sum(obj.X1A,2);
%             ihhB1=sum(obj.X1B,2);
%             ihhA2=sum(obj.X2A,2);
%             ihhB2=sum(obj.X2B,2);
%             ihhA3=sum(obj.X3A,2);
%             ihhB3=sum(obj.X3B,2);
%             ihhA=obj.G*ihhA1+ihhA2+ihhA3;
%             ihhB=obj.G*ihhB1+ihhB2+ihhB3;
            %��ʷ����Դ��һ�μ���
            iHistA=sum(obj.XA,2)-iyA;%�����õ���ʷ����Դ
            iHistB=sum(obj.XB,2)-iyB;
            %% ������װ
            obj.GGTilde=[obj.G,zeros(size(obj.G));zeros(size(obj.G)),obj.G];
%             doubleLambda0=[zeros(size(obj.lambda0)),obj.lambda0;obj.lambda0,zeros(size(obj.lambda0))];
%             doubleLambda1=[obj.lambda1,zeros(size(obj.lambda1));zeros(size(obj.lambda1)),obj.lambda1];
%             A=eye(size(doubleLambda1))+doubleLambda1;%������ϵ��
%             obj.GGTilde=A\(GG-doubleLambda0);%��Ч��·�ĵ絼
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


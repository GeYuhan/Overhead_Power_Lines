classdef node
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    %   �������ڵ絼�����Դ��ÿһ�μ��㶼���������
    %   ������������Ч����Դ�뵼�ɣ�����ά�������룬����ά�ȱ�ʾ�м�����Ч����
    properties
        Nc;%������Ŀ
        v;%����������ѹ
        i;%������������
        k;%ʱ��
    end
    methods
        function obj=node(Nc)%���캯��
            obj.Nc=Nc;
            obj.v=zeros(Nc,1);
            obj.i=zeros(Nc,1);
            obj.k=1;%��ʼֵ
        end
        function [obj,vn,in]=caculation(obj,admittance,source,GGTilde,iHist0)
            [~,~,page]=size(GGTilde);%�����˼�����Ч��·
            k=obj.k;
            %��ѹ����
            Gsum=sum(GGTilde,3)+admittance;
            vn=Gsum\(sum(iHist0,3)+source);
            obj.v(:,k)=vn;
            %��������
            for p=1:page
                in(:,p)=GGTilde(:,:,p)*vn-iHist0(:,:,p);
                obj.i(:,k,p)=in(:,p);
            end
            obj.k=k+1;%����ʱ��
        end
        function [v,i]=getUI(obj)
            v=obj.v;
            i=obj.i;
        end     
    end
    
end


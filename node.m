classdef node
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    %   不储存内电导与电流源，每一次计算都靠外界输入
    %   可以输入外界等效电流源与导纳，以三维矩阵输入，第三维度表示有几个等效输入
    properties
        Nc;%导体数目
        v;%列向量，电压
        i;%列向量，电流
        k;%时刻
    end
    methods
        function obj=node(Nc)%构造函数
            obj.Nc=Nc;
            obj.v=zeros(Nc,1);
            obj.i=zeros(Nc,1);
            obj.k=1;%初始值
        end
        function [obj,vn,in]=caculation(obj,admittance,source,GGTilde,iHist0)
            [~,~,page]=size(GGTilde);%检查接了几个等效电路
            k=obj.k;
            %电压计算
            Gsum=sum(GGTilde,3)+admittance;
            vn=Gsum\(sum(iHist0,3)+source);
            obj.v(:,k)=vn;
            %电流计算
            for p=1:page
                in(:,p)=GGTilde(:,:,p)*vn-iHist0(:,:,p);
                obj.i(:,k,p)=in(:,p);
            end
            obj.k=k+1;%更新时刻
        end
        function [v,i]=getUI(obj)
            v=obj.v;
            i=obj.i;
        end     
    end
    
end


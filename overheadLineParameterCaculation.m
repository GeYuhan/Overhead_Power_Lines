function [Z,Y]=overheadLineParameterCaculation
%% 加载架空线数据
overheadLineParameter
global freq
%% 参数计算
%%求出子导体数目与分裂圆半径
numOfLines=N*numOfPhase;%输电线数目
numOfConductors=numOfLines+numOfGroundLines;%导体数目（包括输电线与地线）
radiusOfCircle=d/2/sin(pi/N);%分裂圆的半径，单位m
%%求出分裂导线子导线相对中心的位置
XrelativeToPhase=cell2mat(arrayfun(@(x) radiusOfCircle*cos(pi/4*(x-1)),1:N,'UniformOutput', false));%行向量
YrelativeToPhase=cell2mat(arrayfun(@(x) radiusOfCircle*sin(pi/4*(x-1)) ,1:N,'UniformOutput',false));
%%求出所有导体（输电线与地线）的位置坐标
XcoordinateOfConductor=ones(N,1)*XcoordinateOfPhase+repmat(XrelativeToPhase.',1,numOfPhase);%矩阵,numOfphase代表3相
XcoordinateOfConductor=[reshape(XcoordinateOfConductor,1,numOfLines),XcoordinateOfGround];%行向量,地线在末尾
YcoordinateOfConductor=ones(N,1)*YcoordinateOfPhase+repmat(YrelativeToPhase.',1,numOfPhase);
YcoordinateOfConductor=[reshape(YcoordinateOfConductor,1,numOfLines),YcoordinateOfGround];%行向量,地线在末尾
%%若考察地线，可将地线位置加在向量末尾，向量长度 numOfLines+numOfGroundLines
%%半径向量与电导率向量生成，在Z矩阵计算时用
radiusOfConductor=[repmat(radiusOfLine,1,numOfLines),repmat(radiusOfGroundLine,1,numOfGroundLines)];%存放所有导体的半径，行向量
conductivityOfConductor=[repmat(conductivityOfLine,1,numOfLines),repmat(conductivityOfGroundLine,1,numOfGroundLines)];%存放所有导体的电导率，行向量
%%土壤磁导率计算
Epsilon0 =8.854187817e-12;% 8.854187817e-12;%真空介电常数，单位F/m
miu0=4*pi*1e-7;%真空磁导率，单位H/m
miu=miur*miu0;%土壤磁导率
%%复频率与透射深度计算
omega=2*pi*freq;%计算ω
s=omega*1i;%s==jω
p=1./sqrt(s*miu*conductivityOfGround);%复透射深度,用到了土壤电导率
%% Y矩阵生成
P=zeros(numOfConductors);%电位系数矩阵
for i=1:numOfConductors
    for j=1:numOfConductors
        Xdifference=(XcoordinateOfConductor(i)-XcoordinateOfConductor(j));
        Ydifference=(YcoordinateOfConductor(i)-YcoordinateOfConductor(j));
        Ysum=(YcoordinateOfConductor(i)+YcoordinateOfConductor(j));
        if i==j
            P(i,j)=1/(2*pi*Epsilon0)*log(Ysum/radiusOfConductor(i));
        else
            P(i,j)=1/(2*pi*Epsilon0)*0.5*log((Xdifference^2+Ysum^2)/(Xdifference^2+Ydifference^2));
        end
    end
end
B=inv(P);%静电感应系数矩阵
Y=arrayfun(@(x) s(x)*B ,1:length(s),'UniformOutput',false);
Y=cat(3,Y{:});
%%给Y补个G
% Y=Y+1e-10;
%% Z矩阵生成
%%carson公式假设：（已改为复镜像）
%%（1）忽略土壤位移电流
%%（2）土壤为半无限大均匀媒质
%%（3）导体是相互平行和无限长的
Z=zeros(numOfConductors,numOfConductors,length(s));%阻抗矩阵,三维

for i=1:numOfConductors
    for j=1:numOfConductors
        Xdifference=(XcoordinateOfConductor(i)-XcoordinateOfConductor(j));
        Ydifference=(YcoordinateOfConductor(i)-YcoordinateOfConductor(j));
        Ysum=(YcoordinateOfConductor(i)+YcoordinateOfConductor(j));
        if i==j
            %%自阻抗
%             %%积分函数定义,复镜像去掉
%             f33=@(u,y,k) exp(-2*y*u)./(u+sqrt(u.^2+1/(p(k)^2)));%被积函数
%             f33(u,YcoordinateOfConductor(i),k)
%             f3=arrayfun(@(k)s(k)*miu0/pi*integral(@(u) f33(u,YcoordinateOfConductor(i),k),0,inf),1:length(s),'UniformOutput',false);%公式第三项，元胞数组
%             selfImpedance3=cat(3,f3{:});%将元胞数组转化为三维矩阵
%             figure(i);
%             plot(squeeze(selfImpedance3));
%             if i==2
%                 E1=selfImpedance3;
%             end
            %%贝塞尔函数项
            f22=@(conductivity,radius,k) sqrt(s(k)*miu/conductivity)/(2*pi*radius)*besseli(0,radius*sqrt(s(k)*miu*conductivity))/besseli(1,radius*sqrt(s(k)*miu*conductivity));
            f2=arrayfun(@(k) f22(conductivityOfConductor(i),radiusOfConductor(i),k),1:length(s),'UniformOutput',false);%这里输入了线路的的电导率与半径
            selfImpedance2(i,j,:)=cat(3,f2{:});%将元胞数组转化为三维矩阵
            %%对数项
            f1=arrayfun(@(k) s(k)*miu0/2/pi*log(2*(YcoordinateOfConductor(i)+p(k))/radiusOfConductor(i)),1:length(s),'UniformOutput',false);
            selfImpedance1(i,j,:)=cat(3,f1{:});%将元胞数组转化为三维矩阵
            %%三项求和
            Z(i,j,:)=selfImpedance1(i,j,:)+selfImpedance2(i,j,:);%三维矩阵
            %将selfImpedance1与selfImpedance2后面加上（i，j，：）是因为缺省，为了查看中间情况
        else
            %%互阻抗
%             %%积分项,由自阻抗的的第三项修改而来
%             f33=@(u,ysum,xdifference,k) cos(u*xdifference).*exp(-ysum*u)./(u+sqrt(u.^2+1/(p(k)^2)));%被积函数
%             f3=arrayfun(@(k)s(k)*miu0/pi*integral(@(u) f33(u,Ysum,Xdifference,k),0,inf),1:length(s),'UniformOutput',false);
%             mutualImpedance3=cat(3,f3{:});%将元胞数组转化为三维矩阵
            %%对数项，由自阻抗的第一项修改而来
            f1=arrayfun(@(k) s(k)*miu0/2/pi*log(sqrt((Xdifference^2+(Ysum+2*p(k))^2)/(Xdifference^2+Ydifference^2))),1:length(s),'UniformOutput',false);
            mutualImpedance1=cat(3,f1{:});%将元胞数组转化为三维矩阵
            %%两项求和
            Z(i,j,:)=mutualImpedance1;%三维矩阵 
            
        end
    end
end
%% 地线消除与分裂导线合并
%%采取老师的办法，可以自适应不同相
%%构造合并矩阵
[P1, P2]=deal(zeros(numOfConductors));%P1左乘矩阵，P2右乘矩阵，numOfphase代表三相
PP=diag(ones(1,N));%构造N单相N分裂导线化简用矩阵
PP(1,:)=ones(1,N);
for n=1:numOfPhase
    index1=(n-1)*N;
    index2=n*N;
    P1(index1+1:index2,index1+1:index2)=PP;
    P2(index1+1:index2,index1+1:index2)=PP.';
end
%%Y合并，不考虑地线
Y=arrayfun(@(x) P1*Y(:,:,x)*P2 ,1:length(s),'UniformOutput',false);%得到的是元胞数组
Y=cat(3,Y{:});%转化为三维矩阵
%%化简后的矩阵（取合并后的元素）
Yphase=zeros(numOfPhase,numOfPhase,length(s));
for i=1:numOfPhase
    for j=1:numOfPhase
        Yphase(i,j,:)=Y((i-1)*N+1,(j-1)*N+1,:);
    end
end

%%Z合并，不考虑地线
inverseZ=arrayfun(@(x) P1*Z(:,:,x)^-1*P2 ,1:length(s),'UniformOutput',false);%得到的是元胞数组
inverseZ=cat(3,inverseZ{:});%转化为三维矩阵
%%%化简后的矩阵（取合并后的元素）
simplifyInverseZ=zeros(numOfPhase,numOfPhase,length(s));
for i=1:numOfPhase
    for j=1:numOfPhase
        simplifyInverseZ(i,j,:)=inverseZ((i-1)*N+1,(j-1)*N+1,:);
    end
end
Zphase=arrayfun(@(x) simplifyInverseZ(:,:,x)^-1 ,1:length(s),'UniformOutput',false);
Zphase=cat(3,Zphase{:});%转化为三维矩阵
%% 序分量
%不考虑换位
if isTransposition
    %考虑换位
    Zs=1/3*(Zphase(1,1,:)+Zphase(2,2,:)+Zphase(3,3,:));%对角线元素
    Zm=1/3*(Zphase(1,2,:)+Zphase(2,3,:)+Zphase(3,1,:));%非对角线元素
    Ys=1/3*(Yphase(1,1,:)+Yphase(2,2,:)+Yphase(3,3,:));
    Ym=1/3*(Yphase(1,2,:)+Yphase(2,3,:)+Yphase(3,1,:));
    Zpositive=Zs-Zm;
    Zzero=Zs+2*Zm;
    Ypositive=Ys-Ym;
    Yzero=Ys+2*Ym;
    for i=1:length(s)
    Zphase(:,:,i)=diag([1 1 1]).*Zs(:,:,i)+(ones(numOfPhase)-diag([1,1,1])).*Zm(:,:,i);
    Yphase(:,:,i)=diag([1,1,1]).*Ys(:,:,i)+(ones(numOfPhase)-diag([1,1,1])).*Ym(:,:,i);
    end
else
    %未换位
    a=exp(2j*pi/3);
    A=[1 1 1;1 a^2 a;1 a a^2];
    for i=1:length(s)        
        Zsym(:,:,i)=A\Zphase(:,:,i)*A;%对称分量参数阻抗矩阵
        Ysym(:,:,i)=A\Yphase(:,:,i)*A;%对称分量参数导纳矩阵
    end
end
%% 输出
Z=Zphase;
Y=Yphase;
end
%% 线路几何参数
%%输电线
%导线型号LGJ-500/35 外径30mm
N=8; % 分裂数
d=0.4;%分裂间距，单位m
numOfPhase=3;%相数
radiusOfLine=15e-3;%导线半径（外径的一半），单位m
%%水平排列三相，从左至右ABC，对地高度30m，两两水平间距26.25m
XcoordinateOfPhase=[-26.56 0 26.56];%排列顺序为A B C
YcoordinateOfPhase=[30 30 30];%排列顺序为A B C

%%地线
numOfGroundLines=2;%地线的数目
radiusOfGroundLine=10e-3;%地线半径（外径的一半），单位m
XcoordinateOfGround=[-28.56 28.56];%地线的横坐标，从左至右
YcoordinateOfGround=[50 50];%地线的纵坐标，从左至右
% %% 线路长度
%  l=100e3;%线路长度，单位m
%% 线路材料参数
resistanceOfLine=0.05812e-3;%导线单位长度电阻，单位Ω/m
sectionAreaOfLine=pi*radiusOfLine^2;%导线截面积，单位m^2
conductivityOfLine=1/resistanceOfLine/sectionAreaOfLine;%导线的电导率,单位S/m

resistanceOfGroundLine=0.3601e-3;%地线单位长度电阻，单位Ω/m
sectionAreaOfGroundLine=pi*radiusOfGroundLine^2;%地线截面积，单位m^2
conductivityOfGroundLine=1/resistanceOfGroundLine/sectionAreaOfGroundLine;%地线的电导率,单位S/m
%% 土壤参数
miur=1;%土壤相对磁导率,magnetic permeability 
conductivityOfGround=0.01;%土壤电导率,单位S/m
%% 是否换位
isTransposition=0;%换位为1，未换位为0
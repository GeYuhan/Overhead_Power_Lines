%% ��·���β���
%%�����
%�����ͺ�LGJ-500/35 �⾶30mm
N=8; % ������
d=0.4;%���Ѽ�࣬��λm
numOfPhase=3;%����
radiusOfLine=15e-3;%���߰뾶���⾶��һ�룩����λm
%%ˮƽ�������࣬��������ABC���Եظ߶�30m������ˮƽ���26.25m
XcoordinateOfPhase=[-26.56 0 26.56];%����˳��ΪA B C
YcoordinateOfPhase=[30 30 30];%����˳��ΪA B C

%%����
numOfGroundLines=2;%���ߵ���Ŀ
radiusOfGroundLine=10e-3;%���߰뾶���⾶��һ�룩����λm
XcoordinateOfGround=[-28.56 28.56];%���ߵĺ����꣬��������
YcoordinateOfGround=[50 50];%���ߵ������꣬��������
% %% ��·����
%  l=100e3;%��·���ȣ���λm
%% ��·���ϲ���
resistanceOfLine=0.05812e-3;%���ߵ�λ���ȵ��裬��λ��/m
sectionAreaOfLine=pi*radiusOfLine^2;%���߽��������λm^2
conductivityOfLine=1/resistanceOfLine/sectionAreaOfLine;%���ߵĵ絼��,��λS/m

resistanceOfGroundLine=0.3601e-3;%���ߵ�λ���ȵ��裬��λ��/m
sectionAreaOfGroundLine=pi*radiusOfGroundLine^2;%���߽��������λm^2
conductivityOfGroundLine=1/resistanceOfGroundLine/sectionAreaOfGroundLine;%���ߵĵ絼��,��λS/m
%% ��������
miur=1;%������Դŵ���,magnetic permeability 
conductivityOfGround=0.01;%�����絼��,��λS/m
%% �Ƿ�λ
isTransposition=0;%��λΪ1��δ��λΪ0
function [ poles ] = initialPoles( frequency,number,option)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
upLimit=log10(frequency(end))+3;%���㷶Χ��������Ƶ��
if option==1
    beta=2*pi*logspace(-1,upLimit,fix(number/2));%�����˼�������
%     beta=2*pi*linspace(frequency(1),frequency(end),fix(number/2));%%�������������
    alpha=beta./-100;
    poles=[alpha+i.*beta,alpha-i.*beta,-1.*ones(number-2*length(beta))];
else
    %ʵ����
    beta=2*pi*logspace(-1,upLimit,number);
%     beta=2*pi*linspace(frequency(2),frequency(end),number);%%�������������
    poles=-beta;
    %�뽫pscad��������ģ����ϳ����й�
%     if number==8
%         poles=[ -5.36051826087555, -1256.36776235495,-25803.7679465045,-123450.955675195,...
%     -428849.477279687,-1378201.70454145,-4452764.55597626,-19690295.1279037];
%     elseif number==12
%         poles=[-23.0342716070336, -121.160092578432, -1133.71630362679, -4626.72047784827,...
%     -15405.7800151817, -44833.0746364702, -118639.475856658, -299713.867826582,...
%     -733016.343706784, -1746723.42496681, -9006409.78198609+3661008.53361868i, -9006409.78198609-3661008.53361868i];
%     end
    
end

end


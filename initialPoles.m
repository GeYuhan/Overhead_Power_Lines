function [ poles ] = initialPoles( frequency,number,option)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
upLimit=log10(frequency(end))+3;%极点范围扩大三个频程
if option==1
    beta=2*pi*logspace(-1,upLimit,fix(number/2));%扩大了极点区间
%     beta=2*pi*linspace(frequency(1),frequency(end),fix(number/2));%%解决单复数问题
    alpha=beta./-100;
    poles=[alpha+i.*beta,alpha-i.*beta,-1.*ones(number-2*length(beta))];
else
    %实极点
    beta=2*pi*logspace(-1,upLimit,number);
%     beta=2*pi*linspace(frequency(2),frequency(end),number);%%解决单复数问题
    poles=-beta;
    %与将pscad极点输入模域拟合程序有关
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


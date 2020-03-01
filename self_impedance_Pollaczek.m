function yy=self_impedance_Pollaczek(x)
global s h
Parameter_cable;
p=1./sqrt(s*mu4/rou4);%%%%∏¥æµœÒ…Ó∂»
zeta=2*h/abs(p);
elta=r7/(2*h);
f1=exp(-zeta.*sqrt(x.^2+1j))./(sqrt(x.^2+1j)+x).*cos(zeta.*elta.*x);
F1=exp(-zeta.*sqrt(1./x.^2+1j))./(sqrt(1./x.^2+1j)+1./x)./x.^2.*cos(zeta.*elta./x);
yy=f1+F1;
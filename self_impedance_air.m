function y=self_impedance_air(u)
global s h
Parameter_cable;

% s=1j*2*pi*f;%%%¸´½ÇÆµÂÊ
alfa0=sqrt(u.^2+s*mu0*(gama0+s*espo0)+(s).^2*mu0*espo0);
alfa1=sqrt(u.^2+s*mu4*(1./rou4+s*espo4)+(s).^2*mu4*espo4);
f1=cos(u*r7)./alfa0.*((mu4*alfa0+mu0*alfa1)+(mu4*alfa0-mu0*alfa1).*exp(-2*alfa0*h))./(mu4*alfa0+mu0*alfa1);

alfa00=sqrt(1./u.^2+s*mu0*(gama0+s*espo0)+(s).^2*mu0*espo0);
alfa11=sqrt(1./u.^2+s*mu4*(1./rou4+s*espo4)+(s).^2*mu4*espo4);
F1=cos(1./u*r7)./alfa00.*((mu4*alfa00+mu0*alfa11)+(mu4*alfa00-mu0*alfa11).*exp(-2*alfa00*h))./(mu4*alfa00+mu0*alfa11)./u.^2;

y=s*mu0/2/pi*(f1+F1);
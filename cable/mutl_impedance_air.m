function y=mutl_impedance_air(u)
global s  x_mul  h1 h2
Parameter_cable;

% s=1j*2*pi*f;%%%¸´½ÇÆµÂÊ
alfa0=sqrt(u.^2+s*mu0*(gama0+s*espo0)+(s).^2*mu0*espo0);
alfa1=sqrt(u.^2+s*mu4*(1./rou4+s*espo4)+(s).^2*mu4*espo4);
f1=cos(u*x_mul)./alfa0.*((mu4*alfa0+mu0*alfa1).*exp(-alfa0.*abs(h1-h2))+(mu4*alfa0-mu0*alfa1).*exp(-alfa0.*(h1+h2)))./(mu4*alfa0+mu0*alfa1);

alfa00=sqrt(1./u.^2+s*mu0*(gama0+s*espo0)+(s).^2*mu0*espo0);
alfa11=sqrt(1./u.^2+s*mu4*(1./rou4+s*espo4)+(s).^2*mu4*espo4);
F1=cos(1./u*x_mul)./alfa00.*((mu4*alfa00+mu0*alfa11).*exp(-alfa00.*abs(h1-h2))+(mu4*alfa00-mu0*alfa11).*exp(-alfa00.*(h1+h2)))./(mu4*alfa00+mu0*alfa11)./u.^2;

y=s*mu0/2/pi*(f1+F1);

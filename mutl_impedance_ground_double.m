function y=mutl_impedance_ground_double(u)
global s  x_mul  h1 h2
Parameter_cable;
alfa0=sqrt(u.^2+s*mu0*(gama0+s*espo0)+(s).^2*mu0*espo0);
alfa1=sqrt(u.^2+s*mu4*(1./rou4+s*espo4)+(s).^2*mu4*espo4);
alfa2=sqrt(u.^2+s*mu5*(1./rou5+s*espo5)+(s).^2*mu5*espo5);
S10=mu4*alfa0+mu0*alfa1;
D10=mu4*alfa0-mu0*alfa1;
S21=mu5*alfa1+mu4*alfa2;
D21=mu5*alfa1-mu4*alfa2;
f1=cos(u*x_mul)./alfa1.*(S10.*S21.*exp(-alfa1.*abs(h1-h2))-S21.*D10.*exp(-alfa1.*(h1+h2))+D21.*S10.*exp(-alfa1.*(2*d1-h1-h2))-D10.*D21.*exp(-alfa1.*(2*d1-abs(h1-h2))))./(D10.*D21.*exp(-2*alfa1*d1)+S10.*S21);

alfa00=sqrt(1./u.^2+s*mu0*(gama0+s*espo0)+(s).^2*mu0*espo0);
alfa11=sqrt(1./u.^2+s*mu4*(1./rou4+s*espo4)+(s).^2*mu4*espo4);
alfa22=sqrt(1./u.^2+s*mu5*(1./rou5+s*espo5)+(s).^2*mu5*espo5);
S10_2=mu4*alfa00+mu0*alfa11;
D10_2=mu4*alfa00-mu0*alfa11;
S21_2=mu5*alfa11+mu4*alfa22;
D21_2=mu5*alfa11-mu4*alfa22;
F1=cos(1./u*x_mul)./alfa11.*(S10_2.*S21_2.*exp(-alfa11.*abs(h1-h2))-S21_2.*D10_2.*exp(-alfa11.*(h1+h2))+D21_2.*S10_2.*exp(-alfa11.*(2*d1-h1-h2))-D10_2.*D21_2.*exp(-alfa11.*(2*d1-abs(h1-h2))))./(D10_2.*D21_2.*exp(-2*alfa11*d1)+S10_2.*S21_2)./u.^2;

y=s*mu4/2/pi*(f1+F1);
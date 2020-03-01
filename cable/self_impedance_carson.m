function y=self_impedance_carson(u)
global s h
Parameter_cable;
pp=1./sqrt(s*mu4/rou4);
f1=exp(-2*h.*u)./(u+sqrt(u.^2+1./pp.^2));
F1=exp(-2*h./u)./(1./u+sqrt(1./u.^2+1./pp.^2))./u.^2;
y=f1+F1;
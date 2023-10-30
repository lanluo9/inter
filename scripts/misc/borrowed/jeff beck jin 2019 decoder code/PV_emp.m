function [prefs] = PV_emp( Y,X )

%empirically estimates the ideal 'preferred stimuli' of a population vector
%decoder by empirically minimizing error.

c=cos(Y')*X;
s=sin(Y')*X;


a=sqrt(c.^2./(c.^2+s.^2)).*sign(c);
b=sqrt(s.^2./(c.^2+s.^2)).*sign(s);

prefs=angle(a'+sqrt(-1)*b');

%PVs=angle(X*cos(a)+sqrt(-1)*X*sin(b))


end


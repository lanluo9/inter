
function [out] = logbesseli(z)
   idx = find(z<200);
   
   out = z-1/2*log(2*pi*z)+log(1+1/8./z+9/2./(8*z).^2+9*25/6./(8*z).^3);
   out(idx)=log(besseli(0,z(idx)));
   
end
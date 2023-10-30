function [ out ] = fcn( a,X )

    n=size(a,1);

    Y=X(:,1);
    X=X(:,2:end);
    
    b=a(n/2+1:end);
    a=a(1:n/2);

    kappa = sqrt((X*a).^2+(X*b).^2);

    out = cos(Y)'*(X*a) + sin(Y)'*(X*b) - sum(logbesseli(kappa));
    out=-out;
    
end

function [out] = logbesseli(z)
   idx = find(z<80);
   
   out = z-1/2*log(2*pi*z)+log(1+1/8./z+9/2./(8*z).^2+9*25/6./(8*z).^3);
   out(idx)=log(besseli(0,z(idx)));
   
end
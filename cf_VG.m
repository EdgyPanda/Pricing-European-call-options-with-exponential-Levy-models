function y = cf_VG(u,T,sigma,nu,theta)

%Lewis formulation. 

%omega =  (1-theta.* nu - sigma.^2 .* nu.*0.5).^(-T./nu);

%temp = (1-1i.*u.*theta.*nu+sigma.^2.*0.5.*nu.*u.^2).^(-T./nu);

%y = exp((1i.*u.*omega.*T)) .* temp;

w  =  1/nu .*log(1 - theta.*nu - 0.5.*nu.*sigma.^2);
y  = exp(1i.*w.*u.*T) .* (1 - 1i.*theta.*nu.*u + .5*sigma^2.*nu.*(u.^2) ) .^ (-T/nu);

end


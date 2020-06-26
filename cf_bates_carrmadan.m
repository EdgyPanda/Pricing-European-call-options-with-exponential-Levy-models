function y = cf_bates_carrmadan(u,S0,T,r,d,V0,theta,kappa,sigma,rho,a,b,lambda)
% Bates
    x = cf_heston(u,S0,T,r,d,V0,theta,kappa,sigma,rho)...
    + cf_jumplognormal(u,a,b,lambda,T);

    y = exp(x);
end
%Schoutens_cf(u,S0,T,r,q,v0,kappa,theta,sigma,rho)
function y = cf_heston(u,S0,T,r,d,V0,theta,kappa,sigma,rho)
    
alfa = -.5*(u.*u + u.*1i);
beta = kappa - rho.*sigma.*u.*1i;
sigma2 = sigma .* sigma;
gamma = .5 .* sigma2;

D = sqrt(beta .* beta - 4.0 .* alfa .* gamma);

bD = beta - D;
eDt = exp(- D .* T);

G = bD ./ (beta + D);
B = (bD ./ sigma2) .* ((1.0 - eDt) ./ (1.0 - G .* eDt));
psi = (G .* eDt - 1.0) ./(G - 1.0);
A = ((kappa .* theta) / (sigma2)) .* (bD .* T - 2.0 .* log(psi));

y = A + B.*V0 + 1i.*u.*(log(S0)+(r-d).*T);

end

function yJump = cf_jumplognormal(u,a,b,lambda,T)
% LogNormalJump for Merton and Bates
    a=log(1+a)-0.5.*b.^2;

    yJump = lambda.*T.*(-a.*u.*1i + (exp(u.*1i.*log(1.0+a)+0.5.*b.*b.*u.*1i.*(u.*1i-1.0))-1.0));
end

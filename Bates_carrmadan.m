function [call] = Bates_carrmadan(S0,K,T,r,q,V0,kappa,theta,sigma,rho,alpha,b,lambda)
%Carr madan approach:

a = 0.75;


 %cf_bates_carrmadan(u,S0,T,r,d,V0,theta,kappa,sigma,rho,a,b,lambda)
density = @(u) (exp(-r.*T).*cf_bates_carrmadan(u-(a+1).*1i,S0,T,r,q,V0, theta, kappa, sigma, rho, alpha,b,lambda))...
    ./(a.^2+a-u.^2+1i.*(2.*a+1).*u);


integrand = @(u) real(exp(-1i.*u.*log(K)).*density(u));

call = (exp(-a.*log(K))./pi) .* integral(integrand,0,inf,'ArrayValued',true); 


end
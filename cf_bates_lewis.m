function [x] = cf_bates_lewis(u,T,v0,theta,kappa,sigma,rho,a,b,lambda)

    x = cf_jumplognormal(u,a,b,lambda,T) .* Schoutens_cf_lewis(u,T,v0,kappa,theta,sigma,rho);
    

end




function yJump = cf_jumplognormal(u,a,b,lambda,T)
% LogNormalJump for Merton and Bates
    %yJump = exp(lambda*T*(-1i.*u.*a.*(exp(u*1i*log(1.0+a)+0.5*b*b*u*1i.*(u*1i-1.0))-1.0)));
    w = -lambda.*(exp(a+b.^2 .* 0.5) - 1.0);
    yJump = exp(1i.*u.*T.*w + lambda.*T.*(exp(1i.*u.*a - u.^2.*b.^2.*0.5) - 1.0));
   %Brix %yJump = exp(1i.*u.*T-T.*lambda.* (exp(a+0.5.*b.^2)-1.0)) .* exp(-1i*u.*T.*lambda.*(exp(a+0.5.*b.^2)-1)+T.*lambda.*(exp(1i.*u.*a+0.5.*b.^2.*u.^2)-1));

   
end

function [cf] = Schoutens_cf_lewis(u,T,v0,kappa,theta,sigma,rho)
%---------------------------correction term:-------------------------------
epsiloncor = kappa - sigma.*rho;

dcor = sqrt((-epsiloncor).^2);

g2cor = (epsiloncor-dcor)./(epsiloncor+dcor);

innerterm2cor = (1-g2cor.*exp(-dcor.*T))./(1-g2cor);

term2cor = ((kappa.*theta)./(sigma.^2)) .* ((epsiloncor-dcor).*T - 2.* log(innerterm2cor));

term3cor = (v0./sigma.^2) .* (epsiloncor-dcor) .* ((1-exp(-dcor.*T))./(1-g2cor.*exp(-dcor.*T)));

w = term2cor+term3cor;

%--------------------------------------------------------------------------
epsilon = kappa - sigma.*rho.*1i.*u;

d = sqrt((-epsilon).^2 - sigma.^2.*(-u.^2-1i.*u));

g2 = (epsilon-d)./(epsilon+d);

term1 = 1i.*u.*T.*w;

innerterm2 = (1-g2.*exp(-d.*T))./(1-g2);

term2 = ((kappa.*theta)./(sigma.^2)) .* ((epsilon-d).*T - 2.* log(innerterm2));

term3 = (v0./sigma.^2) .* (epsilon-d) .* ((1-exp(-d.*T))./(1-g2.*exp(-d.*T)));

cf = exp(term1+term2+term3);

end
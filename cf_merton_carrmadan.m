function y = cf_merton_carrmadan(u,S0,T,r,d,sigma,a,b,lambda)
% Merton Jump Diffusion
    x = cf_bs(u,S0,T,r,d,sigma) ...
        + cf_jumplognormal(u,a,b,lambda,T);
    y = exp(x);
    
end


function y = cf_bs(u,S0,T,r,d,sigma)
    %y = exp(i*u*(lnS+(r-d-0.5*sigma*sigma)*T) - 0.5*sigma*sigma*u.*u*T);
    y = 1i.*u.*(log(S0)+(r-d-0.5.*sigma.*sigma).*T) - 0.5.*sigma.*sigma.*u.*u.*T;
end

function yJump = cf_jumplognormal(u,a,b,lambda,T)
% LogNormalJump for Merton and Bates
    yJump = lambda.*T.*(-a.*u.*1i + (exp(u.*1i.*log(1.0+a)+0.5.*b.*b.*u.*1i.*(u.*1i-1.0))-1.0));
end
function y = cf_merton(u,T,sigma,a,b,lambda)
% Merton Jump Diffusion
    x = cf_bs(u,T,sigma) +  ...
        cf_jumplognormal(u,a,b,lambda,T);
    
    y = exp(x);%exp(temp(u,T,sigma,a,b,lambda));
    %y = exp(temp(u,T,sigma,a,b,lambda));
end

%Setting lambda =0 we get the Black-Scholes equation which works. 
%Correct
function y = cf_bs(u,T,sigma)
    %y = exp(i*u*(lnS+(r-d-0.5*sigma*sigma)*T) - 0.5*sigma*sigma*u.*u*T);
    y = 1i.*u.*T.*(-0.5.*sigma.^2) - 0.5.*sigma.^2.*u.^2.*T;
end


function yJump = cf_jumplognormal(u,a,b,lambda,T)
% LogNormalJump for Merton and Bates
    %yJump = lambda*T*(-1i.*u.*a.*(exp(u*1i*log(1.0+a)+0.5*b*b*u*1i.*(u*1i-1.0))-1.0));
    w = -lambda.*(exp(a+b.^2 .* 0.5) - 1.0);
    yJump = 1i.*u.*T.*w + lambda.*T.*(exp(1i.*u.*a - u.^2.*b.^2.*0.5) - 1.0);
end




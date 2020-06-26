function cpl = skew_lewis(T, impvol,func, varargin)
%setting z0=0.5.
%setting integral limits from 0 to 1000. 

%S = repmat(S, size(K,1),1);


int = @(u) tempintegrand(u,func,T,varargin{:});
%int = @(u) u.^0.00001.*k; %just a tester
lewisintegral = integral(int,0,145000,'ArrayValued',true);

cpl = -exp(impvol.^2.*0.1250.*T).*sqrt(2./pi).*1./(sqrt(T)).* lewisintegral;

end

function [ret] = tempintegrand(u,func,varargin)

ret = u.*imag((feval(func, ...
    u-0.5*1i, varargin{:})))./(u.*u+.25);

end


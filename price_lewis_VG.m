function cpl = price_lewis_VG(S, K, T, r, d, varargin)
%setting z0=0.5.
%setting integral limits from 0 to 1000. 
%Contains "func" to specify any CF. 

%S = repmat(S, size(K,1),1);

k = log(S./K) + (r-d) .* T;

int = @(u) tempintegrand(u,k,T,varargin{:});
%int = @(u) u.^0.00001.*k; %just a tester
lewisintegral = integral(int,0,5000,'ArrayValued',true);

cpl = S.*exp(-d.*T)-(sqrt(S.*K)./pi).*exp(-(r+d).*T.*0.5) ...
    .* lewisintegral;

end

function [ret] = tempintegrand(u,k, varargin)

ret = real(exp(1i.*u.*k) ...
    .*(feval(@cf_VG, ...
    u-0.5*1i, varargin{:})))./(u.*u+.25);

end


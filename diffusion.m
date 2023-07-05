function out = diffusion(t,D,A,C,nmax,rres,rlim,sigma)
j0 = besselzero(0,nmax,1); % get nmax zeros of 0th order bessel function of first kind

c0 = zeros(1,nmax); % constants defined by boundary condition
for ii = 1:nmax
    c0(ii) = (C*2)/(j0(ii)*besselj(1,j0(ii)));
end

r = 0:rres:rlim;


if isscalar(t)
    u = zeros(1,numel(r));
    for ii = 1:nmax
        u = u + c0(ii).*besselj(0,j0(ii)/A.*r).*exp(-(j0(ii)/A)^2*D.*t);
    end
    u = (C-u).*exp(-r.^2/(2*sigma^2)).*r*...
        (-1/(sigma^2*(exp(-rlim^2/(2*sigma^2)) - 1)));
    out = trapz(r,u);
else
    [R,T] = meshgrid(r,t);
    u = zeros(numel(t),numel(r));
    for ii = 1:nmax
        u = u + c0(ii).*besselj(0,j0(ii)/A.*R).*exp(-(j0(ii)/A)^2*D.*T);
    end
    u = (C-u).*exp(-R.^2/(2*sigma^2)).*r*...
        (-1/(sigma^2*(exp(-rlim^2/(2*sigma^2)) - 1)));
    out = trapz(r,u,2);
end
end
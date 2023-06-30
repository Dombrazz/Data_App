function [u,t,r] = diffusionModel(A,C,nmax,D,r,t)

j0 = besselzero(0,nmax,1); % get nmax zeros of 0th order bessel function of first kind

c0 = zeros(1,nmax); % constants defined by boundary condition
for ii = 1:nmax
    c0(ii) = (C*2)/(j0(ii)*besselj(1,j0(ii)));
end

if nargin == 4
    nt = nmax;
    nr = nmax;
    Dt = 10;
    r = (0:nr-1)*A/nr;
    t = (0:nt-1)*Dt/nr;
    [R,T] = meshgrid(r,t);
    u = zeros(nt,nr);
    for ii = 1:nmax
        u = u + c0(ii).*besselj(0,j0(ii)/A.*R).*exp(-(j0(ii)/A)^2*D.*T);
    end
    
    u = -u+C;
elseif nargin == 6
    u = 0;
    for ii = 1:nmax
        u = u + c0(ii).*besselj(0,j0(ii)/A*r).*exp(-(j0(ii)/A)^2*D*t);
    end
    u = -u+C;
end
end
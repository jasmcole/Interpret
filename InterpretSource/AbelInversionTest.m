clear
%y = 0:0.001:2;
y = linspace(0,2,1000);
N = length(y);
dy = y(2) - y(1);
r = y;
f = 1 - y;
f(f<0) = 0;

%f = f + 0.5*(rand(size(f)) - 0.5);

Af = -(2/pi)*log(r./(1 + sqrt(1 - r.^2)));
Af(~isreal(Af)) = 0;
Af = real(Af);

dfdy = gradient(f,1);
dfdy = dfdy/dy;

disp(['Numerical integration:'])
tic

for r_index = 1:length(r)
    for y_index = r_index:length(y)
        y_int = (double(y_index) - 0.9)*dy;
        r_int = (double(r_index) - 1.0)*dy;
        integrand(y_index - r_index + 1) = dfdy(y_index)*(y_int^2 - r_int^2)^-0.5;
    end
    Af_num(r_index) = -(2/pi)*dy*trapz(integrand);
end

toc

disp(['Bessel decomposition:'])
tic

k = 0:0.1:1000;
dk = k(2) - k(1);

for k_index = 1:length(k)
    k_int = k(k_index);
    Ff(k_index) = (2/pi)*trapz(cos(k_int*y).*f*dy);
end

% fac = 4;
% frackeep = .1;
% Ff = dct(f,fac*length(f));
% k = linspace(0,pi/dy/fac,length(y));
% dk = k(2) - k(1);
% k = k(1:round(frackeep*N));
% Ff = Ff(1:round(frackeep*N));
% Ff(1) = Ff(1)/sqrt(2);
% Ff = Ff*sqrt(pi/N/(4/fac));

% Ff = real(fft(f));
% Ff = Ff(1:round(end/2));
% k = linspace(0,2*pi/dy,length(Ff));
% dk = k(2) - k(1);

for r_index = 1:length(r)
    r_int = (double(r_index) - 1)*dy;
    Af_bessel(r_index) = trapz(besselj(0, k*r_int).*k.*Ff*dk);
end

%Af_bessel = Af_bessel + pi;
%Af_bessel = Af_bessel/sqrt(length(y));

toc


%plot(r,Af,r,Af_num,r,Af_bessel); ylim([0 5])
plot(r,Af_num - Af,r,Af_bessel - Af); ylim([-0.2 0.2])

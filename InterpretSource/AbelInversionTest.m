clear
y = 0:0.01:2;
dy = y(2) - y(1);
r = y;
f = 1 - y;
f(f<0) = 0;

Af = -(2/pi)*log(r./(1 + sqrt(1 - r.^2)));
Af(~isreal(Af)) = 0;

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

k = 0:0.01:100;
dk = k(2) - k(1);

for k_index = 1:length(k)
    k_int = k(k_index);
    Ff(k_index) = (2/pi)*trapz(cos(k_int*y).*f*dy);
end

for r_index = 1:length(r)
    r_int = (double(r_index) - 1)*dy;
    Af_bessel(r_index) = trapz(besselj(0, k*r_int).*k.*Ff*dk);
end

toc


plot(r,Af,r,Af_num,r,Af_bessel)

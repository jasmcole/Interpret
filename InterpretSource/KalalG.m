%rho from 0 to 1
function g = KalalG(k,rhovec)

for n = 1:length(rhovec)
    rho = rhovec(n);
    if(rho == 1)
        g(n) = 0;
    elseif(rho == 0)
        g(n) = (2/pi)*sinint(pi*k);
    else
        t = linspace(0,sqrt(1-rho^2),1000);
        integrand = (2/pi)*(1./sqrt(t.^2 + rho^2)).*sin(k*pi*sqrt(t.^2 + rho.^2));
        g(n) = trapz(t,integrand);
    end
end

end
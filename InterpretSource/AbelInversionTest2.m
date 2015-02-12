clear
load('g.mat')

tnum = 0;
tbessel = 0;

%for line = 1:1000

y = linspace(0,2,100);
R = max(y);
dy = y(2) - y(1);
r = y;
f = 1 - y;
f(f<0) = 0;
f = f + 0.5*(rand(size(f)) - 0.5);
%im(:,line) = f;

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
tnum = tnum + toc;
%imnum(:,line) = Af_num;
disp(['Bessel decomposition:'])

%k = 0:0.1:1000;
%dk = k(2) - k(1);
% for k_index = 1:length(k)
%     k_int = k(k_index);
%     Ff(k_index) = (2/pi)*trapz(cos(k_int*y).*f*dy);
% end

x = y/R;
Ff = dct(f);
k = 1:10;
dctinds = k+1;
ak = 2*Ff(dctinds)/sqrt(2*length(x));

Af_bessel = zeros(size(x));

tic
for n = 1:length(k)
   Af_bessel  = Af_bessel + (pi/2)*k(n)*ak(n)*g(n,:);
end

tbessel = tbessel + toc;
%imbessel(:,line) = Af_bessel;

% subplot(121)
% imagesc(imnum)
% subplot(122)
% imagesc(imbessel)
% drawnow
subplot(121)
fs14
plot(y,f)
ylim([-0.1 1.1])
subplot(122)
fs14
plot(r,Af,r,Af_num,r,Af_bessel); ylim([-0.5 5])
%plot(r,Af_num - Af,r,Af_bessel - Af); ylim([-0.2 0.2])
%end

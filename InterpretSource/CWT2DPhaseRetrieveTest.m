clear
x = linspace(-0.5,0.5,480);
y = linspace(-0.5,0.5,640);
kx = linspace(-pi/(x(2) - x(1)),pi/(x(2) - x(1)), length(x));
ky = linspace(-pi/(y(2) - y(1)),pi/(y(2) - y(1)), length(y));
[xg yg] = meshgrid(x,y);
Nfringes = 50;
lambda = (max(x) - min(x))/Nfringes;
alo = lambda/5;
ahi = 3*lambda;
phireal = 2*pi*xg/lambda + 30*exp(-yg.^2/0.2^2).*exp(-(xg - 0.5).^2/0.6^2);
phireal = phireal + 10*exp(-(yg - 0.3).^2/0.01^2).*exp(-(xg + 0.2).^2/0.1^2);
phireal = phireal - min(min(phireal));
I = 1 + cos(phireal) + 2*perlin_noise(size(phireal)) + 0.3*randn(size(phireal));
phireal = phireal - 2*pi*xg/lambda;
m = 2;
n = 2;

subplot(3,2,5)
imagesc(phireal); fs14; title('\phi_{real}')
subplot(3,2,6)
imagesc(I); fs14; title('I_{real}')
drawnow
%%

theta = linspace(-pi/2, 0.99*pi/2, 11);
a = linspace(alo, ahi ,15);
Wold = zeros(size(xg));
thetapeak = Wold;
apeak = Wold;

for tind = 1:length(theta)
    for aind = 1:length(a)
        zeta = exp(-(m/(a(aind)^2))*(xg.^2 + yg.^2) + (1i*2*pi/a(aind))*(xg*cos(theta(tind)) + yg*sin(theta(tind))));
        Wnew = (1/(a(aind)^2))*ifft2(fft2(I).*fft2(zeta));
        newinds = (abs(Wnew) > abs(Wold));
        thetapeak(newinds) = theta(tind);
        apeak(newinds) = a(aind);
        Wold(newinds) = Wnew(newinds); 
    end
end

apeak = fftshift(apeak);
thetapeak = fftshift(thetapeak);
W = fftshift(Wold);

phi = angle(W);
phi = unwrap(phi')';

phi = phi - 2*pi*xg/lambda;

subplot(3,2,1)
imagesc(apeak); fs14; title('\lambda_{calc}')
subplot(3,2,2)
imagesc(thetapeak); fs14; title('\theta_{calc}')
subplot(3,2,3)
imagesc(phi); fs14; title('\phi_{calc}')
subplot(3,2,4)
imagesc(real(W)); fs14; title('I_{calc}')
subplot(3,2,5)
imagesc(phireal); fs14; title('\phi_{real}')
subplot(3,2,6)
imagesc(I); fs14; title('I_{real}')

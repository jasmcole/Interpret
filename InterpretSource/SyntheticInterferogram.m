clear
ncrit = 1750e24;
lambda = 800e-9;
ymax = 2000e-6;
Ny = 300;
y = linspace(-ymax, ymax, Ny);
dy = y(2) - y(1);
x = 0:dy:1e-2;
Nx = length(x);
nmin = 1e23;
nmax = 5e24;
Rmin = 0;
Rmax = 1000e-6;
k = 2*pi/(40*dy);

nmaxarray = linspace(nmin, nmax, Nx);
Rarray = linspace(Rmin, Rmax, Nx);

for n = 1:Nx
    nmax = nmaxarray(n);
    R = Rarray(n);
    
    phi(:,n) = (2*pi/lambda)*(nmax/ncrit)*( 0.5*sqrt(R^2 - y.^2) - (0.5*y.^2/R).*log(R./abs(y) + sqrt(R^2./y.^2 - 1)) );
    phi(:,n) = real(phi(:,n));
    fringes(:,n) = k*x(n)*ones(Ny,1);
    ne(:,n) = nmax*(1 - abs(y)/R);
end

ne(ne < 0) = 0;

interferogram = 0.5*(1 + cos(fringes + phi));
interferogram = interferogram + 0*rand(size(interferogram));
interferogram = interferogram/max(max(abs(interferogram)));
reference = 0.5*(1 + cos(fringes));
reference = reference + 0*rand(size(reference));
reference = reference/max(max(abs(reference)));

subplot(1,2,1)
imagesc(x,y,phi)
subplot(1,2,2)
imagesc(x,y,interferogram)

imwrite(interferogram, 'SyntheticInterferogram.tiff', 'tiff')
imwrite(reference, 'SyntheticReference.tiff', 'tiff')

[rho1 peakdensity] = AnalyseInterferogram('SyntheticInterferogram.tiff', 'SyntheticReference.tiff', 'CWT');

ydiff = 0.5*(length(ne(:,1)) - length(rho1(:,1)));
xdiff = 0.5*(length(ne(1,:)) - length(rho1(1,:)));
ne1 = ne(1+ydiff:end-ydiff,1+xdiff:end-xdiff);

[rho2 peakdensity] = AnalyseInterferogram('SyntheticInterferogram.tiff', 'SyntheticReference.tiff', 'FFT');

ydiff = 0.5*(length(ne(:,1)) - length(rho2(:,1)));
xdiff = 0.5*(length(ne(1,:)) - length(rho2(1,:)));
ne2 = ne(1+ydiff:end-ydiff,1+xdiff:end-xdiff);

[rho3 peakdensity] = AnalyseInterferogram('SyntheticInterferogram.tiff', 'SyntheticReference.tiff', 'Hilbert');

ydiff = 0.5*(length(ne(:,1)) - length(rho3(:,1)));
xdiff = 0.5*(length(ne(1,:)) - length(rho3(1,:)));
ne3 = ne(1+ydiff:end-ydiff,1+xdiff:end-xdiff);

figure
subplot(1,3,1)
contour(ne1, linspace(nmin, nmax, 10))
title('CWT')
hold on
contour(rho1, linspace(nmin, nmax, 10))
subplot(1,3,2)
contour(ne2, linspace(nmin, nmax, 10))
title('FFT')
hold on
contour(rho2, linspace(nmin, nmax, 10))
subplot(1,3,3)
contour(ne3, linspace(nmin, nmax, 10))
title('Hilbert')
hold on
contour(rho3, linspace(nmin, nmax, 10))



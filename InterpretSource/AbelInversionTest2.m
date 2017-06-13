clear
% Contains pre-computed basis functions
load('/Users/jmc208/Dropbox/Experiments/2016TA2/Server Files/Interpret/InterpretSource/g.mat')

% Co-ordinate in measurement space
y = linspace(0,2,100);
R = max(y);
dy = y(2) - y(1);
% Co-ordinate in Abel-inverted space
r = y;
% Test input with analytic solution
f = 1 - y;
f(f < 0) = 0;
freal = f;
% Add some noise
f = f + 0.1*(rand(size(f)) - 0.5);

% Analytic Abel inversion
Af = -(2/pi)*log(r./(1 + sqrt(1 - r.^2)));
Af(~isreal(Af)) = 0;
Af = real(Af);

dfdy = gradient(f,1);
dfdy = dfdy/dy;

% Direct numerical calculation of Abel inversion
for r_index = 1:length(r)
    for y_index = r_index:length(y)
        y_int = (double(y_index) - 0.9)*dy;
        r_int = (double(r_index) - 1.0)*dy;
        integrand(y_index - r_index + 1) = dfdy(y_index)*(y_int^2 - r_int^2)^-0.5;
    end
    Af_num(r_index) = -(2/pi)*dy*trapz(integrand);
end

% Normalised co-ordinate
x = y/R;
% Discrete cosine transform of input f
Ff = dct(f);
% Number of modes to use
k = 1:10;
% Shift index by 1 for some reason
dctinds = k+1;
% Calculate coefficients
ak = 2*Ff(dctinds)/sqrt(2*length(x));

% Fourier Abel inversion
Af_bessel = zeros(size(x));

for n = 1:length(k)
   Af_bessel  = Af_bessel + (pi/2)*k(n)*ak(n)*g(n,:);
end

subplot(121)
plot(y,freal,y,f)
ylim([-0.1 1.1])
xlabel('$y$'), ylabel('$f(y)$')
legend('Actual input', 'Noisy input')
FormatPlot();
subplot(122)
plot(r,Af,r,Af_num,r,Af_bessel); ylim([-0.5 5])
xlabel('$r$'), ylabel('$\mathcal{A}^{-1}f(r)$')
legend('Actual solution', 'Direct numerical inversion', 'Fourier inversion')
FormatPlot();


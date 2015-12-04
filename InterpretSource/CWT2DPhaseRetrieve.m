function phase = CWT2DPhaseRetrieve(handles)

% From Jason 20/01/2014
% Make sure fringes are vertical in image and you've put the right number
% of fringes in (should be automated later...)

datafile = double(handles.dataimage);
intref = double(handles.refimage);
calibdata = handles.calibdata;

Nfringes = calibdata.nfringes;

[ysize xsize] = size(datafile);

x = linspace(-0.5,0.5,xsize);
y = linspace(-0.5,0.5,ysize);

[xg yg] = meshgrid(x,y);

m = 3;
n = 2;

lambda = (max(x) - min(x))/Nfringes;
alo = lambda*calibdata.CWT_fringelo;
ahi = lambda*calibdata.CWT_fringehi;

theta = linspace(-pi/2, 0.99*pi/2, calibdata.CWT_ntheta);
a = linspace(alo, ahi,calibdata.CWT_nlambda);

Wold = zeros(size(xg));
Woldref = zeros(size(xg));
thetapeak = Wold;
thetapeakref = Woldref;
apeak = Wold;
apeakref = Woldref;

for tind = 1:length(theta)
    for aind = 1:length(a)
        zeta = exp(-(m/(a(aind)^2))*(xg.^2 + yg.^2) + (1i*2*pi/a(aind))*(xg*cos(theta(tind)) + yg*sin(theta(tind))));
        Wnew = (1/(a(aind)^2))*ifft2(fft2(datafile).*fft2(zeta));
        Wnewref = (1/(a(aind)^2))*ifft2(fft2(intref).*fft2(zeta));
        newinds = (abs(Wnew) > abs(Wold));
        newindsref = (abs(Wnewref) > abs(Woldref));
        thetapeak(newinds) = theta(tind);
        thetapeakref(newindsref) = theta(tind);
        apeak(newinds) = a(aind);
        apeakref(newindsref) = a(aind);
        Wold(newinds) = Wnew(newinds);
        Woldref(newindsref) = Wnewref(newindsref);
    end    
    set(handles.StatusBox, 'String', ['2D CWT ' num2str(round(100*tind/length(theta))) '% complete']); drawnow
end

apeak = fftshift(apeak);
apeakref = fftshift(apeakref);
thetapeak = fftshift(thetapeak);
thetapeakref = fftshift(thetapeakref);
W = fftshift(Wold);
Wref = fftshift(Woldref);

phase = angle(W);
phaseref = angle(Wref);

%Take ends off to avoid phase errors
phase(:,1:10) = [];
phase(:,end-9:end) = [];
phaseref(:,1:10) = [];
phaseref(:,end-9:end) = [];

% subplot(3,2,1)
% imagesc(apeak); fs14; title('\lambda_{calc}')
% subplot(3,2,2)
% imagesc(thetapeak); fs14; title('\theta_{calc}')
% subplot(3,2,3)
% imagesc(phi); fs14; title('\phi_{calc}')
% subplot(3,2,4)
% imagesc(real(W)); fs14; title('I_{calc}')
% subplot(3,2,6)
% imagesc(I); fs14; title('I_{real}')

phase = phase - phaseref;

axes(handles.PhasediagAxes)
imagesc(real(W))
axis image xy
axes(handles.PhaseAxes)
imagesc(phase)
axis image xy
drawnow

end
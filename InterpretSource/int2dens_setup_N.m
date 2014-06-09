%function [Phi, xim, yim, mag, lineoutstruct] = int2dens_setup_N(file,reffile, shotdate)
function [Phi, xim, yim, mag,P] = int2dens_setup_N(file,reffile, shotdate,useoldarea,useoldfilter)
%Nick's version of function getting rid of unnecessary bits etc

% call, for example
%   quickintdens('datafile.tiff','reffile.tiff','20090527',0)
% and follow instructions
%
% FOR TAP: I use flippicture = 1. That flips up to down, so in probe1if, the laser still comes from the left.
%

expensiveunwrapflag = 0;
%expensiveunwrapflag = 'GoldStein';    % or 'Quality' or 'GoldStein', 0 for the usual
%expensiveunwrapflag = 'Quality';    % or 'Quality' or 'GoldStein', 0 for the usual
rotateflag = 0;
plasmaflag = 0;     % use 0 for neutral gas shots
flipphase = 0;      % 0 if don't flip, 1 if do flip
%useoldarea = 0;    % 0 for redefine, 1 for use the one saved in cutarea, -1 for use the full image
%useoldfilter = 0;

% neutral gas measurements - note inconsequential unless using this to get
% out density measurements
ngas = 0.000281;    % Argon    actually this is already ngas-1
%ngas = 0.000132;    % Hydrogen
%ngas = 0.000035;    % Helium
l_probe = .633;     % red HeNe
%l_probe = .532;     % greeny

if nargin<3 || isempty(shotdate)
    shotdate = [];
end

%if plasmaflag
%    lambda_probe = .8; % GEMINI
%    lambda = .8;
%else
    lambda_probe = l_probe;
    lambda = l_probe;
    %end


if plasmaflag
    ngas = 0;
end

% find date and use TCC wire shots to place images in the same coordinate system
if isempty(shotdate)
    I = findstr(file,'\');
    if isempty(I); I=0; end;
    shotdate = file(I(end)+1:I(end)+8);
end
[~, xshift_i, direction, micronperpixel, flippicture] = getimageshift(shotdate);    % get the xaxis right from positions of the needle ...

% ref image
%Aref = double(imread(reffile,'tiff'));
if isstr(reffile)
    Aref = double(imread(reffile));
else
    Aref = double(reffile);
end
% data image
%Adat = double(imread(file,'tiff'));
if isstr(file)
    Adat = double(imread(file));
else
    Adat = double(file);
end

% Adat(Adat>75) = 75;
% Aref(Aref>75) = 75;
% 
% 
% Adat(Adat<35) = 35;
% Aref(Aref<35) = 35;
% %figure,imagesc(Adat)

% use laser direction
if direction == -1
    Adat = flipdim(Adat,2);
    Aref = flipdim(Aref,2);
    xshift_i = size(Aref,2) - xshift_i;
end

if flippicture
    Adat = flipdim(Adat,1);
    Aref = flipdim(Aref,1);
end

% axis and quantities
nc = 8.8e-12*9.1e-31/1.602e-19^2*(2*pi*3e8/(lambda_probe*1e-6))^2/1e6;        % I was lazy

% cut relevant region and rotate the image
if useoldarea==0
    % select relevant region which is later unwrapped
    figure
    imagesc(Adat)
    if rotateflag
        title('select line for rotation, click and double click')
        [xr, yr] = getline;
        xr = xr(1:2); yr = yr(1:2);
        [xr, I] = sort(xr);
        yr = yr(I);
        theta = 180/pi*atan((yr(2)- yr(1)) / (xr(2)- xr(1)) );
        imagesc(imrotate(Adat, theta))
    else
        theta = NaN;
    end
    title('\Delta\Phi, Select region for unwrapping');
    %imagesc(xim,yim,abs(P))
    %title('Transmission, Select region for unwrapping');
    disp('Select region for unwrapping');
    [~,rect]=imcrop;
    close
    rect = round(rect);
    if rect(1)<1; rect(1) = 1; end      % in case the selected rectangle is too big
    if rect(1)+rect(3)>size(Adat,2); rect(3) = size(Adat,2)-rect(1); end
    if rect(2)<1; rect(2) = 1; end
    if rect(2)+rect(4)>size(Adat,1); rect(4) = size(Adat,1)-rect(2); end
    save(['cutarea.mat'],'rect','theta');
elseif useoldarea==1
    load(['cutarea.mat']);
elseif useoldarea == -1
    rect = [1 1 size(Adat,2)-1 size(Adat,1)-1];
    theta = 0;
end


% axis according to the roation (or not)
if isnan(theta) % no rotation
    xim = ((0:size(Adat,2)-1)-xshift_i)*micronperpixel;
    yim = (0:size(Adat,1)-1)*micronperpixel;
else
    tmp = abs(theta/180*pi);
    Dx = (cos(tmp)*size(Adat,2) + sin(tmp)*size(Adat,1));
    Dy = (sin(tmp)*size(Adat,2) + cos(tmp)*size(Adat,1));
    Aref = imrotate(Aref, theta);
    Adat = imrotate(Adat, theta);
    xim = ((0:size(Adat,2)-1)-xshift_i)*micronperpixel* Dx/size(Adat,2);
    yim = (0:size(Adat,1)-1)*micronperpixel* Dy/size(Adat,1);
end

% cut data, rect in pixel
xim = xim(rect(1):rect(1)+rect(3));
yim = yim(rect(2):rect(2)+rect(4));

% put y to start at zero ... this is helpful in most cases (assumes we always cut from the nozzle)
yim = yim - yim(1);

Aref = Aref(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));
Adat = Adat(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));
Adat = Adat - mean(mean(Adat));
Aref = Aref - mean(mean(Aref));

% Fourier Transformation (2D)
% zero padding
mrows = 2^(ceil(log2(size(Aref,1)))  + 0);
mcols = 2^(ceil(log2(size(Aref,2)))  + 0);
rstart = floor((mrows - size(Aref,1))/2);   % number of upper rows of zeros
cstart = floor((mcols - size(Aref,2))/2);

Areftmp = zeros(mrows, mcols);
Areftmp(rstart+1 : rstart+size(Aref,1), cstart+1 : cstart+size(Aref,2) ) = Aref;

Adattmp = zeros(mrows, mcols);
Adattmp(rstart+1 : rstart+size(Adat,1), cstart+1 : cstart+size(Adat,2) ) = Adat;

ref_ftx = fftshift(fft2(Areftmp));
dat_ftx = fftshift(fft2(Adattmp));

% filter fringes
if useoldfilter
    load(['fourierfilter.mat']);
else
    [tmp,Filter] = filterfringes(dat_ftx);
    save fourierfilter Filter
end
dat_ftx = Filter.*dat_ftx;
ref_ftx = Filter.*ref_ftx;

% do the backtransform
ref_iftx = ifft2(fftshift(ref_ftx));
dat_iftx = ifft2(fftshift(dat_ftx));

% get rid of the extra zeros from fft
ref_iftx = ref_iftx( rstart+1 : rstart+size(Aref,1), cstart+1 : cstart+size(Aref,2) );
dat_iftx = dat_iftx( rstart+1 : rstart+size(Aref,1), cstart+1 : cstart+size(Aref,2) );

% % % % % % % % % figure,imagesc(real(ref_iftx))
% % % % % % % % % figure,imagesc(real(dat_iftx))

%dat_iftx=Adat;
%ref_iftx = Aref;

% complex phase factor
I = find(ref_iftx~=0);
P = zeros(size(dat_iftx));
P(I) = dat_iftx(I)./ref_iftx(I);
% figure,imagesc(real(P))
% unwrap phase
if plasmaflag
%    Phi1 = unwrap(angle(dat_iftx),[],1);
%    Phi2 = unwrap(angle(ref_iftx),[],1);
%     Phi1 = angle(dat_iftx);
%     Phi2 = angle(ref_iftx);
%     Phi1 = flipdim(Phi1,2);   
%     Phi2 = flipdim(Phi2,2);   
    
%     figure
%     subplot(2,2,1)
%     imagesc(Phi1)
%     title('Phidat not unw')
%     subplot(2,2,2)
%     imagesc(Phi2)
%     title('Phiref not unw')
% 
%     Phi1 = unwrap(Phi1);
%     Phi2 = unwrap(Phi2);
%     tmp = mean(Phi1(1:5,:)); % this is the baseline
%     Phi1 = Phi1 - tmp(ones(size(Phi1,1),1),:); % substract baseline
%     tmp = mean(Phi2(1:5,:)); % this is the baseline
%     Phi2 = Phi2 - tmp(ones(size(Phi2,1),1),:); % substract baseline
%     
%     subplot(2,2,3)
%     imagesc(Phi1)
%     title('Phidat')
%     subplot(2,2,4)
%     imagesc(Phi2)
%     title('Phiref')
%     
%     return
    
%    Phi = Phi1 - Phi2;
    %%%%% This was "working"
    if ~expensiveunwrapflag
        tmp = angle(P);
        Phi = unwrap( tmp, [], 1 );
    elseif strcmp(expensiveunwrapflag,'Quality')
        Phi = QualityGuidedUnwrap2D(P)
    elseif strcmp(expensiveunwrapflag,'GoldStein')
        Phi = GoldsteinUnwrap2D(P);
    end    
%     figure
%     imagesc(Phi)
%     title('select point to start re-unwrapping');
%     [xtmp, ytmp] = ginput(1);
%     Phi1 = Phi(:,1:round(xtmp));
%     Phi2 = Phi(:,round(xtmp)+1:end);
%     Phi1 = unwrap( flipdim(Phi1,2), [], 2);
%     Phi1 = flipdim(Phi1,2);
%     Phi2 = unwrap(Phi2,[],2);
%     Phi = [Phi1 Phi2];
%    Phi = Phi1 - Phi2;

    nomp = 10;  % number of columns to average over
    % take care of constant gradient
%    constgrad = (mean(Phi(end - nomp+1:end,:), 1) - mean(Phi(1:nomp,:), 1)) / size(Phi,1);    % gradient along the columns
%    xtmp = (1:size(Phi,1))';
%    Phi = Phi - constgrad(ones(size(Phi,1),1),:) .* xtmp(:,ones(size(Phi,2),1));

    tmp = mean(Phi(1:5,:)); % this is the baseline
    Phi = Phi - tmp(ones(size(Phi,1),1),:); % substract baseline
else
    Phi = angle(P);
%    Phi = flipdim(Phi,2);
     Phi = unwrap(Phi,[],2);
%    Phi = flipdim(Phi,2);
    
%    Phi = GoldsteinUnwrap2D(P);
%    Phi = unwrap(Phi,[],1);
%    Phi = GoldsteinUnwrap2D(P);

    nomp = 10;  % number of columns to average over
    % take care of constant gradient
    constgrad = (mean(Phi(:,end - nomp+1:end), 2) - mean(Phi(:,1:nomp), 2)) / size(Phi,2);    % gradient along the columns
    xtmp = 1:size(Phi,2);
    Phi = Phi - constgrad(:, ones(size(Phi,2),1)) .* xtmp(ones(size(Phi,1),1),:);

    tmp = mean(Phi(:,end-nomp+1:end),2); % this is the baseline for neutral density shots
    Phi = Phi - tmp(:,ones(size(Phi,2),1)); % substract baseline
end

if flipphase
    Phi = -Phi;
end

% plot stuff
% % % % % % % % % % % % % % % % % % % % % % % figure
% % % % % % % % % % % % % % % % % % % % % % % subplot(2,2,1)
% % % % % % % % % % % % % % % % % % % % % % % imagesc(xim,yim,Aref)
% % % % % % % % % % % % % % % % % % % % % % % title('reference')
% % % % % % % % % % % % % % % % % % % % % % % axis tight
% % % % % % % % % % % % % % % % % % % % % % % colorbar

% % % % % % % % % % % % % % % % % % % % % % % % % % % subplot(2,2,2)
% % % % % % % % % % % % % % % % % % % % % % % % % % % imagesc(xim,yim,Adat)
% % % % % % % % % % % % % % % % % % % % % % % % % % % title('data')
% % % % % % % % % % % % % % % % % % % % % % % % % % % axis tight
% % % % % % % % % % % % % % % % % % % % % % % % % % % colorbar

mag = (abs(P));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % subplot(2,2,3)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % imagesc(xim,yim,(abs(P)));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %title('log_{10}T')
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % title('Transmission')
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % xlabel('x [\mum]');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % ylabel('y [\mum]');
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % colorbar

% % % % % % % % % % % % % % % % % % % % % % % % % % subplot(2,2,4)
% % % % % % % % % % % % % % % % % % % % % % % % % % imagesc(xim,yim,Phi);
% % % % % % % % % % % % % % % % % % % % % % % % % % title('\Delta\Phi')
% % % % % % % % % % % % % % % % % % % % % % % % % % xlabel('x [\mum]');
% % % % % % % % % % % % % % % % % % % % % % % % % % ylabel('y [\mum]');
% % % % % % % % % % % % % % % % % % % % % % % % % % colorbar

Phicut = Phi;
% calc dens and stuff from lineouts, uncomment below if necessary

% [Pmax1, dchannel1, ne1, Philine, xl] =lineoutcalc(xim,yim,Phicut,nc,lambda, plasmaflag, ngas);
% disp(['max Phaseshift: ' num2str(Pmax1)]);
% disp(['mean diameter [µm]: ' num2str(dchannel1)]);
% disp(['mean ne: ' num2str(ne1/1e18)]);
% 
% figure
% plot(xl,Philine);
% xlabel('r [µm]')
% ylabel('\Delta \Phi [rad]')
% 
% lineoutstruct.Pmax = Pmax1;
% lineoutstruct.dchannel = dchannel1;
% lineoutstruct.ne = ne1;
% lineoutstruct.Phi = Philine;
% lineoutstruct.x = xl;

%if isstr(file)
%    save([file(1:end-4) '.mat'], 'Phi', 'xim', 'yim', 'mag', 'lineoutstruct');
%end

if isstr(file)
    save([file(1:end-4) '.mat'], 'Phi', 'xim', 'yim', 'mag');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pmax, dchannel, ne, Philineout, lineoutpos] =lineoutcalc(xim,yim,Phicut,nc,lambda, plasmaflag, ngas)

% replot
h=imagesc(xim,yim,Phicut);
colorbar
xlabel('y [\mum]')
ylabel('\Delta \Phi')
title('select lineout, click & click')
[xline,yline]=getline;
xl = linspace(xline(1), xline(2), 100);
yl = linspace(yline(1), yline(2), 100);

% calculate phase along the lineout
[Xim, Yim]=meshgrid(xim,yim);
Philineout = interp2(Xim,Yim, Phicut, xl, yl);

dl = sqrt( (xl(2)-xl(1))^2 + (yl(2)-yl(1))^2 );     % step length along lineout in micron
lineoutpos = (0:length(xl)-1)*dl;
Pmax = max(Philineout);
if Pmax>0
    I = find(Philineout>Pmax/2);
else
    I = find(Philineout<Pmax/2);
end
Pfwhm = length(I)*dl;   % this is the full width at half maximum

dchannel = 2/sqrt(3)*Pfwhm;
refind = lambda/2/pi*Pmax/(dchannel);
if plasmaflag
    ne = nc*(1-(1-refind).^2);
else
    ne = 2.69e19/ngas  * (refind);
end






function [Indrow,Indcol] = selectcenterline(Phi,xim,yim,micronperpixel)
figure
imagesc(xim,yim,Phi)
title('Select line along center of channel, Finish with ENTER');
disp('Select line along center of channel, Finish with ENTER');
[Xc,Yc]=ginput;
[Xc,I]=sort(Xc);
Yc = Yc(I);
Xc = Xc/micronperpixel;
Yc = Yc/micronperpixel;
Indcol = ceil(min(Xc)) : floor(max(Xc));
Indrow = zeros(size(Indcol));
for a = 1:length(Xc)-1
    x1 = Xc(a); x2 = Xc(a+1);
    y1 = Yc(a); y2 = Yc(a+1);
    I = find(Indcol>=x1 & Indcol<=x2);
    m = (y2-y1)/(x2-x1);
    Indrow(I) = round(m*(Indcol(I)-x1)+y1);
end
hold on
plot(Indcol*micronperpixel,Indrow*micronperpixel,'k','LineWidth',1)
hold off
disp('Press ENTER to continue')
pause
close
disp('Think, think, think, ...')


function [xn,yn,NEMAT] = phi2dens(Phi, Irow, Icol, resolution, lambda_probe, nc, micronperpixel)
plotflag = 0;
noangles = 18;
if plotflag
    figure
end
NEMAT = NaN*ones(max(Icol)*2, ceil( length(Irow)/resolution ));
b = 0;
for a = 1:resolution:length(Irow)-resolution
    b=b+1;
    
    meanrow = round(mean(Irow(a:a+resolution)));
    
    projected = lambda_probe/2/pi* mean([Phi(1:meanrow,Icol(a):Icol(a)+resolution); Phi(meanrow:-1:1,Icol(a):Icol(a)+resolution)],2);
    projected = projected-min(projected);
    I = iradon(projected(:,ones(noangles,1)), 180/noangles); 
    lineout = I(ceil(end/2),:);
    ne = nc*(1-(1-lineout).^2);
    centershift = size(Phi,1)-ceil(mean(Irow(a : a+resolution)));
    
    NEMAT(ceil(end/2)-centershift-ceil(length(ne)/2)+1 : ceil(end/2)-centershift+floor(length(ne)/2),b) = ne(:);
    if plotflag
        subplot(2,1,1)
        plot(projected,'k')
        hold on
        back = radon(I,0);
        plot(back,'b');
        hold off
    
        subplot(2,1,2)
        plot(ne,'r')
        pause
    end
end
if plotflag
    close
end
% delete NaN's in NEMAT
I = find( sum(~isnan(NEMAT),2) );
NEMAT = NEMAT(I,:);
I = find( sum(~isnan(NEMAT),1) );
NEMAT = NEMAT(:,I);

xn = [1:size(NEMAT,2)]*resolution*micronperpixel;
yn = (1:size(NEMAT,1))*micronperpixel;
xn=xn(:);
yn = yn(:);


function [diameter, nemean, lowB, highB] = meanvalues(NEMATc)
NEMATc(find(isnan(NEMATc))) = 0;
intNEcum = cumtrapz(NEMATc, 1);    % cumulative integral over ne along r
intNE = intNEcum(end,:);       % complete integral over ne along r
[I,J] = find( intNEcum > .25*intNE(ones(size(intNEcum,1),1),:) & intNEcum < .75*intNE(ones(size(intNEcum,1),1),:));
lowB = zeros(size(NEMATc,2),1);
highB = zeros(size(NEMATc,2),1);
for a = 1:size(NEMATc,2)
    Itmp = find(J==a);
    if ~isempty(Itmp)
        lowB(a) = min(I(Itmp));
        highB(a) = max(I(Itmp));
    else
        lowB(a) = NaN;
        highB(a) = NaN;
    end
end
diameter = highB - lowB;
nemean = intNE(:)./diameter/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Phic = myunwrap2(Phi)
% SLOW !!!!!!!!!!!!!
% once again go through the columns and look for jumps larger than pi
Phi1 = unwrap(Phi,[],1);
Phi2 = unwrap(Phi(end:-1:1,:),[],1);
Phi2 = Phi2(end:-1:1,:);
tmp1 = Phi2 - Phi1;
tmp2 = abs(abs(Phi1) - abs(Phi2));

Ierr = find( abs(tmp1(round(end/2),:) > pi));   % problematic regions

tmpdiff = zeros(size(tmp2,1)-1, size(tmp2,2));
tmpdiff(:,Ierr) = diff(tmp2(:,Ierr));   % differences allong rows of tmp2
[Mjump, Ijump] = max(abs(tmpdiff));

%tmp3 = -ones(size(tmp2));
%tmp3(:,Ierr) = tmp2(:,Ierr);
%[Mjump, Ijump] = min(tmp3,[],1);
%Ijump( find(Mjump==-1) ) = [];           % delete indecees where everything is ok already

Phi3 = Phi1;
for a = 1:length(Ierr)
    Phi3(Ijump(a):end, Ierr(a)) = Phi2(Ijump(a):end, Ierr(a));
end

figure
subplot(2,2,1)
imagesc(Phi1);
subplot(2,2,2)
imagesc(Phi2);
subplot(2,2,3)
imagesc(Phi3);
subplot(2,2,4)
imagesc(tmp2);

Phic = Phi3;




function Phic = myunwrap(Phi,cut)
% once again go through the columns and look for jumps larger than 'cut'
ahead = 1;  % if phi jumps over 2*pi within this number of steps, correct 
cut = 3*pi/2;
vertdir = 1;
if vertdir; Phi = Phi'; end
for a = 1:ahead:size(Phi,2)-ahead
    difftmp = Phi(:,a+ahead) - Phi(:,a);
    I = find(abs(difftmp)>cut);
    if ~isempty(I)
        difftmp2 = difftmp(I,ones(ahead,1));
        Phi(I,a+1 : a+ahead) = Phi(I,a+1 : a+ahead) - 2*pi*sign(difftmp2);
    end
end
if vertdir; Phi = Phi'; end
Phic = Phi;


%function Phic = correctjumps(Phi,cut)
%ahead = 5;
%for a = ahead:ahead:size(Phi,2)-ahead
%    difftmp = sum(Phi(:,a+1:a+ahead) - Phi(:,a+1-ahead:a),2)/ahead;
%    I = find(abs(difftmp)>cut);
%    if ~isempty(I)
%        difftmp2 = difftmp(I,ones(ahead,1));
%        Phi(I,a+1:a+ahead) = Phi(I,a+1:a+ahead) - cut*fix(difftmp2/(cut));
%    end
%end
%Phic = Phi;
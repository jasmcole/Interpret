function rho = FitToPlateau(phi, calibdata, gas, lambda, plotflag)

h = calibdata.micperpix*0.000001;
lambda = lambda*1e-9;
rho_atm = (101.325e3*6.022e23)/(8.31*273);

switch gas
    case 'Argon'
        k = 0.000281;
        disp(['Assuming Argon gas with n - 1 = ' num2str(k) ' at atmospheric pressure'])
        disp(['Laser wavelength is ' num2str(lambda*1e9) 'nm'])
    case 'Plasmaold'
        k = (1.6e-19)^2*lambda^2*rho_atm/( 8*pi^2 * 9.11e-31 * 8.85e-12 * 9e16);
        disp(['Assuming Plasma'])
        disp(['Laser wavelength is ' num2str(lambda*1e9) 'nm'])
    case 'Plasma'
        %ncrit for 800nm in m^-3
        k = rho_atm/(2*1748e24);
        disp(['Assuming Plasma'])
        disp(['Laser wavelength is ' num2str(lambda*1e9) 'nm'])
end

phi = fliplr(phi);

ymidautoflag = 1;
[zsize ysize] = size(phi);

if (ymidautoflag == 1)
    % Auto ymid finder
    yaxis = 1:ysize;
    zpoints = 0.1:0.1:0.9;
    zpoints = zpoints*zsize;
    zpoints = round(zpoints);
    for n = 1:length(zpoints)
        yres(n) = trapz(phi(zpoints(n),:).*yaxis)./trapz(phi(zpoints(n),:));
    end
    ymid = round(mean(yres))
    imagesc(phi)
    line([ymid ymid], [0 zsize], 'color', 'white')
    hold on
    scatter(yres, zpoints)
    hold off
    drawnow
else
    % Interactive ymid finder
    ymidfig = figure;
    imagesc(phi)
    [ymid zmid] = ginput(1);
    ymid = round(ymid)
end

phihalf = phi(:, ymid:end);
phihalf = -phihalf*sign(mean(mean(phihalf)));

y = 0:length(phihalf(1,:))-1;
y = y*h;
ymax = max(y);

options = optimset('TolFun',1e-10);

for z_index = 1:zsize,
    
    philine = phihalf(z_index,:);
    philine = philine - mean(philine(end-30:end));
    
    R1 = trapz(philine.*y)/trapz(philine);
    R2 = 2*R1;
    n0 = (lambda*rho_atm*philine(1))/(6*pi*k*R1);
    A = 4*pi*k*n0/(rho_atm*lambda);
    
    params0 = [A R1 R2];
    
    paramsfit = lsqnonlin(@PhiFit,params0,[],[],options,y,philine);
    
    phifit = PhiAnalytic(paramsfit, y);
    
    subplot(1,2,1)
    plot(y,philine,y,phifit)
    drawnow
    
    n0 = paramsfit(1)*rho_atm*lambda/(4*pi*k);
    R1 = paramsfit(2);
    R2 = paramsfit(3);
    
    rholine = RhoAnalytic(y, n0, R1, R2);
    
    subplot(1,2,2)
    plot(y,rholine)
    drawnow
    
    rho(z_index, :) = rholine;
    
end

rho = [fliplr(rho) rho];
rho = imrotate(rho, -90);
[zaxis yaxis] = size(rho);
yaxis = 1:yaxis;
zaxis = 1:zaxis;
yaxis = 1e3*h*yaxis;
zaxis = 1e3*h*zaxis;

if (plotflag == 1)
    figure
    subplot(2,2,1)
    imagesc(imrotate(phi,-90))
    axis image
    subplot(2,2,2)
    imagesc(yaxis, zaxis, rho)
    xlabel('x /mm')
    ylabel('y /mm')
    axis image xy
    subplot(2,2,3)
    plot(zaxis, mean(rho')/1e6)
    
end

end

function phi = PhiAnalytic(params, yvec)

[A R1 R2] = deal(params(1), params(2), params(3));

for n = 1:length(yvec)
    y = yvec(n);
    if (y < R1)
        phi(n) = A*(sqrt(R1^2 - y^2) + ...
            (R2/(R2 - R1))*(sqrt(R2^2 - y^2) - sqrt(R1^2 - y^2)) - ...
            (0.5/(R2 - R1))*( R2*sqrt(R2^2 - y^2) - R1*sqrt(R1^2 - y^2) + y^2*log( (R2 + sqrt(R2^2 - y^2))/(R1 + sqrt(R1^2 - y^2)))));
    else
        phi(n) = A*( ( (R2/(R2 - R1))*(sqrt(R2^2 - y^2)) ) - (0.5/(R2 - R1))*( R2*sqrt(R2^2 - y^2) + y^2*log( (R2 + sqrt(R2^2 - y^2))/y ) ) );
    end
end

phi = real(phi);

end

function phidiff = PhiFit(params, yvec, phi)

[A R1 R2] = deal(params(1), params(2), params(3));

for n = 1:length(yvec)
    y = yvec(n);
    if (y < R1)
        phidiff(n) = A*(sqrt(R1^2 - y^2) + ...
            (R2/(R2 - R1))*(sqrt(R2^2 - y^2) - sqrt(R1^2 - y^2)) - ...
            (0.5/(R2 - R1))*( R2*sqrt(R2^2 - y^2) - R1*sqrt(R1^2 - y^2) + y^2*log( (R2 + sqrt(R2^2 - y^2))/(R1 + sqrt(R1^2 - y^2)))));
    else
        phidiff(n) = A*( ( (R2/(R2 - R1))*(sqrt(R2^2 - y^2)) ) - (0.5/(R2 - R1))*( R2*sqrt(R2^2 - y^2) + y^2*log( (R2 + sqrt(R2^2 - y^2))/y ) ) );
    end
end

phidiff = real(phidiff);

phidiff = phidiff - phi;

end

function rho = RhoAnalytic(y, n0, R1, R2)

for n = 1:length(y)
    if (y(n) < R1)
        rho(n) = n0;
    else
        rho(n) = n0*(R2 - y(n))/(R2 - R1);
    end
end

rho(rho < 0) = 0;

end
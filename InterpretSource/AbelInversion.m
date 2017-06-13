%Performs Abel inversion of phase map and outputs density profile
%micperpix is in microns, gas is a string, lambda is in nm

function [rho, yaxis, zaxis, rhomean] = AbelInversion(Phi, calibdata, gas, lambda, plotflag, handles)

set(handles.StatusBox, 'String', 'Beginning Abel inversion'); drawnow

Phi(isnan(Phi)) = 0;

if (calibdata.ymid == 0)
    
    ymidautoflag = 0;
    
    if (ymidautoflag == 1)
        % Auto ymid finder
        [zsize ysize] = size(Phi);
        yaxis = 1:ysize;
        zpoints = 0.1:0.1:0.9;
        zpoints = zpoints*zsize;
        zpoints = round(zpoints);
        for n = 1:length(zpoints)
            yres(n) = trapz(Phi(zpoints(n),:).*yaxis)./trapz(Phi(zpoints(n),:));
        end
        ymid = round(mean(yres))
        imagesc(Phi)
        line([ymid ymid], [0 zsize], 'color', 'white')
        hold on
        scatter(yres, zpoints)
        hold off
        drawnow
    else
        % Interactive ymid finder
        set(handles.StatusBox, 'String', 'Select symmetry axis'); drawnow
        imagesc(Phi)
        [ymid zmid] = ginput(1);
        ymid = round(ymid);
        set(handles.StatusBox, 'String', ['Selected ymid = ' num2str(ymid)]); drawnow
        line([ymid ymid], [0 1e4], 'color', 'white'); drawnow
        %arrowline([ymid 1.3*ymid], [zmid zmid]); drawnow
    end
    
end

if (calibdata.ymid > 0)
    ymid = calibdata.ymid;
end

ysize = length(Phi(1,:))-1;
zsize = length(Phi(:,1));
h = calibdata.micperpix * 1e-6;
asmooth = calibdata.asmooth;

lambda = lambda*1e-9;
ncrit = 1748e24 * (800e-9/lambda)^2;

%n=1+k(rho/rho_atm)
%rho_atm=(P/RT)*Na
rho_atm=(101.325e3*6.022e23)/(8.31*273); %units m^-3

% From http://www.kayelaby.npl.co.uk/general_physics/2_5/2_5_7.html
switch gas
    case 'Argon'
        k = 0.000281;
    case 'Helium'
        k = 0.000035;
    case 'Nitrogen'
        k = 0.000298;
    case 'Plasma'
        k = rho_atm/(2*ncrit);
end

Ny = ysize - ymid + 1;
Philine = zeros(1, Ny);
Philinemirror = zeros(1, 2*Ny); 
yaxis = 1:Ny;
yaxis2 = 1:2*Ny;
rho = zeros(zsize, ysize-ymid);

t = 0;

for z_index = 1:zsize,
    
    Philine = Phi(z_index, ymid:ysize);
    Philineold = Philine;
    Philinemirror = [fliplr(Philine) Philine];
    Philinemirror = csaps(yaxis2, Philinemirror, asmooth, yaxis2);
    Philine = Philinemirror(end/2 + 1:end);

    if (~mod(z_index, round(zsize/10)))
        % This plots the smoothed and unsmoothed phase for comparison
        axes(handles.DensitydiagAxes)
        plot(Philineold, '.')
        hold on
        plot(Philine, 'red')
        hold off
        percentdone = [num2str(floor(100*z_index/zsize)) '%'];
        set(handles.StatusBox, 'String', ['Doing Abel inversion...' percentdone]); drawnow
    end
    
    dPhidy = -gradient(Philine,1)/h;
    dPhidy(1) = 0.1*dPhidy(2);
        
    [yg rg] = meshgrid(([1:Ny] - 0.9)*h, (0:Ny-2)*h);
    dPhig = repmat(dPhidy, [Ny-1, 1]);
    integrand = dPhig./sqrt(yg.^2 - rg.^2);
    integrand = integrand.*(rg < yg);
    rho(z_index, :) = sum(integrand');
end

rho = rho * ((lambda*rho_atm)/(2*pi^2*k))*h;

rho = [fliplr(rho) rho];
rho = imrotate(rho, -90);
[zaxis yaxis] = size(rho);
zaxis = -zaxis/2:zaxis/2 - 1;
yaxis = 1:yaxis;
yaxis = 1e3*h*yaxis;
zaxis = 1e3*h*zaxis;

if (plotflag == 1)
    axes(handles.DensityAxes)
    imagesc(yaxis, zaxis, rho)
    xlabel('z /mm')
    ylabel('y /mm')
    axis image xy
    colorbar
    axes(handles.DensitydiagAxes)
    rhomean = mean(rho(round(end/2)-10:round(end/2)+10,:));
    plot(yaxis, rhomean)
    xlabel('z /mm')
    ylabel('Density /m^{-3}')
    title('Mean density within central 20 pixels')
end

set(handles.StatusBox, 'String', 'Finished Abel inversion'); drawnow

end


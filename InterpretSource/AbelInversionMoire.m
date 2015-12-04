%Performs Abel inversion of phase map and outputs density profile
%micperpix is in microns, gas is a string, lambda is in nm

function [rho, yaxis, zaxis, rhomean] = AbelInversion(Phi, calibdata, gas, lambda, plotflag, handles)
set(handles.StatusBox, 'String', 'Beginning Abel inversion'); drawnow
%Phi = fliplr(Phi);

Phi(isnan(Phi)) = 0;

[ysize xsize] = size(Phi);
xsize = 1:xsize;
ysize = 1:ysize;
x = {ysize,xsize};
[Phi,p] = csaps(x,Phi,1,x);

if (calibdata.ymid == 0)
    % Interactive ymid finder
    set(handles.StatusBox, 'String', 'Select symmetry axis'); drawnow
    imagesc(Phi)
    [ymid zmid] = ginput(1);
    ymid = round(ymid);
    set(handles.StatusBox, 'String', ['Selected ymid = ' num2str(ymid)]); drawnow
    line([ymid ymid], [0 1e4], 'color', 'white'); drawnow
else
    ymid = calibdata.ymid;
end

ysize = length(Phi(1,:))-1;
zsize = length(Phi(:,1));
h = calibdata.micperpix*0.000001;
asmooth = calibdata.asmooth;

lambda = lambda*1e-9;

% n = 1 + k(rho/rho_atm)
% rho_atm=(P/RT)*Na
rho_atm=(101.325e3*6.022e23)/(8.31*273); %units m^-3

switch gas
    case 'Argon'
        k=0.000281;
    case 'Nitrogen'
        k=0.000265;
    case 'Plasmaold'
        k = (1.6e-19)^2*lambda^2*rho_atm/( 8*pi^2 * 9.11e-31 * 8.85e-12 * 9e16);
    case 'Plasma'
        %ncrit for 800nm in m^-3
        k = rho_atm/(2*1748e24);
end
count = 0;
for z_index = 1:zsize,
    
    Philine = Phi(z_index,:);
    Philine = Philine(ymid:ysize);
    Philineold = Philine;
    Philine = [-fliplr(Philine) Philine];
    yaxis = 1:length(Philine);
    Philine = csaps(yaxis, Philine, asmooth, yaxis);
    Philine = Philine(end/2 + 1:end);
    yaxis = 1:length(Philine);
    
    if (count > zsize/10)
        count = 0;
        % This plots the smoothed and unsmoothed phase for comparison
        axes(handles.DensitydiagAxes)
        plot(Philineold, '.')
        hold on
        plot(Philine, 'red')
        hold off
        percentdone = [num2str(floor(100*z_index/zsize)) '%'];
        set(handles.StatusBox, 'String', ['Doing Abel inversion...' percentdone]); drawnow
    end
    
    dPhidy = Philine;
    dPhidy(1) = 0;
    
    for r_index = 1:ysize-ymid,
        for y_index = r_index:ysize - ymid + 1
            y = (double(y_index) - 0.9)*h;
            r = (double(r_index) - 1.0)*h;
            integrand(y_index - r_index + 1) = dPhidy(y_index)*(y^2 - r^2)^-0.5;
            
        end
        
        rho(z_index, r_index) = (rho_atm/(2*pi^2*k))*(1e-3*calibdata.Moire_p_um/calibdata.Moire_d_mm)*h*trapz(integrand);
        integrand(end) = [];
    end
    
    count = count + 1;
    
end

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
end

set(handles.StatusBox, 'String', 'Finished Abel inversion'); drawnow

end


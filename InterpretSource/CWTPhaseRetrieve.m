function phaser = CWTPhaseRetrieve(handles)

data = double(handles.dataimage);
intref = double(handles.refimage);
calibdata = handles.calibdata;

[Nrows Ncols] = size(data);

rng default % Reset random number generator for reproducible results

t = linspace(0,1,Ncols);
dt = t(2) - t(1);

[ysize xsize] = size(data);
lambda = xsize/calibdata.nfringes;

count = 0;
for row = 1:Nrows
    
    set(handles.StatusBox, 'String', ['CWT ' num2str(round(100*row/Nrows)) '% complete']); drawnow
    
    y = data(row,:);
    yref = intref(row,:);
    y = y - mean(y);
    yref = yref - mean(yref);
    y = y/max(abs(y));
    yref = yref/max(yref);
    sig = struct('val',y,'period',dt);
    sigref = struct('val',yref,'period',dt);
    
    MorletFourierFactor = 2/(6+sqrt(2+6^2));
    
    %     omegamin = 200;
    %     omegamax = 350;
    %     nscales = 100;
    %     omegas = linspace(omegamin, omegamax, nscales);
    %     scales = Ncols./(omegas*MorletFourierFactor)
    
    scales = linspace(lambda*calibdata.CWT_fringelo,calibdata.CWT_fringehi*lambda,calibdata.CWT_nlambda);
    omegas = 1./scales;
    
    coeffs = cwt(y,scales,'cmor1-1');
    refcoeffs = cwt(yref,scales,'cmor1-1');
    
    for n = 1:length(t)
        [val index] = max(coeffs(:,n));
        [valref indexref] = max(refcoeffs(:,n));
        omegar(n) = omegas(index);
        phase(n) = angle(coeffs(index,n));
        phaseref(n) = angle(refcoeffs(indexref,n));
        yr(n) = real(coeffs(index,n));
        yrref(n) = real(refcoeffs(index,n));
    end
    
    phaseref(isnan(phaseref)) = 0;
    datar(row,:) = yr;
    datarref(row,:) = yrref;
    phaser(row,:) = unwrap(phase - phaseref);
    phaser = unwrap(phaser);
    
    yr = yr/max(yr);
    
    count = count + 1;
    if (count/Nrows > 0.1)
        if (row > 1)
            axes(handles.PhasediagAxes)
            pcolor(t,scales,abs(coeffs)); shading flat
            hold on
            scatter(t,1./omegar)
            hold off
            xlabel('Position in Image')
            ylabel('Wavelengths /pixels')
            
            drawnow
        end
        count = 0;
    end
    
end

phaser = -(imrotate(phaser, 90));
phaser = phaser';
phaser = fliplr(phaser);
axes(handles.PhaseAxes)
imagesc(phaser); axis image xy; colorbar

end
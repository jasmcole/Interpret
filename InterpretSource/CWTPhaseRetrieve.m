function [phaser, mask] = CWTPhaseRetrieve(handles)

data = double(handles.dataimage);
intref = double(handles.refimage);
calibdata = handles.calibdata;
maskthresh = calibdata.CWT_maskthresh;

[Nrows, Ncols] = size(data);
mask = ones(size(data));

t = linspace(0,1,Ncols);

[ysize, xsize] = size(data);
lambda = xsize/calibdata.nfringes;
lambdalo = calibdata.CWT_lambdalo;
lambdahi = calibdata.CWT_lambdahi;

count = 0;

for row = 1:Nrows
    
    set(handles.StatusBox, 'String', ['CWT ' num2str(round(100*row/Nrows)) '% complete']); drawnow
    
    y = data(row,:);
    yref = intref(row,:);
    y = y - mean(y);
    yref = yref - mean(yref);
    y = y/max(abs(y));
    yref = yref/max(yref);

    scales = linspace(lambdalo, lambdahi, calibdata.CWT_nlambda);
    omegas = 1./scales;
    
    coeffs = cwt(y,scales,'cmor1-1');
    refcoeffs = cwt(yref,scales,'cmor1-1');
    
    for n = 1:length(t)
        [val, index] = max(coeffs(:,n));
        [valref, indexref] = max(refcoeffs(:,n));
        omegar(n) = omegas(index);
        
        mask(row, n) = abs(val)*sqrt(omegar(n));
        phase(n) = angle(coeffs(index,n));
        phaseref(n) = angle(refcoeffs(indexref,n));
        yr(n) = real(coeffs(index,n));
        yrref(n) = real(refcoeffs(index,n));
    end
    
    phaseref(isnan(phaseref)) = 0;
    datar(row,:) = yr;
    datarref(row,:) = yrref;
    phaser(row,:) = (phase - phaseref);
    
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


mask = mask/max(mask(:));

axes(handles.PhasediagAxes)
imagesc(mask); axis image xy; colorbar; title('Phase mask')
hold on
contour(mask, [maskthresh, maskthresh], 'Color', 'white')
hold off

mask(mask < maskthresh) = 0;
mask(mask > 0.0) = 1;

phaser = phaser.*(-mask);
axes(handles.PhaseAxes)
imagesc(phaser); axis image xy; colorbar

end
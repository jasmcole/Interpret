function phase = HilbertPhaseRetrieve(handles)

data = double(handles.dataimage);
intref = double(handles.refimage);
calibdata = handles.calibdata;

[nrows ncols] = size(data);

phase = zeros(size(data));
mask = zeros(size(data));
xaxis = 1:ncols;

set(handles.StatusBox, 'String', 'Doing Hilbert transform...'); drawnow
count = 0;
for n = 1:nrows
    y = data(n,:);
    yplot = y;
    yref = intref(n,:);
    movav = csaps(xaxis, y, calibdata.hsmooth1, xaxis);
    movavref = csaps(xaxis, yref, calibdata.hsmooth1, xaxis);
    y = y - movav;
    yref = yref - movavref;
    y = csaps(xaxis, y, calibdata.hsmooth2, xaxis);
    ysmooth = y;
    yref = csaps(xaxis, yref, calibdata.hsmooth2, xaxis);
    Hy = hilbert(y);
    Hyref = hilbert(yref);
    mag(n,:) = abs(Hy);
    magref(n,:) = abs(Hyref);
    y = y./abs(Hy);
    yref = yref./abs(Hyref);
    Hy = hilbert(y);
    Hyref = hilbert(yref);
    angleHyref = angle(Hyref);
    angleHyref(isnan(angleHyref)) = 0;
    phase(n,:) = (angle(Hy) - angleHyref);
    imorig(n,:) = y;
    imref(n,:) = yref;
    mask(n,:) = (yplot - (ysmooth+movav));
    count = count + 1;
    if (count > nrows/10)
        count = 0;
        percentdone = [num2str(floor(100*n/nrows)) '%'];
        axes(handles.PhasediagAxes)
        x = 1:length(y);
        plot(x,yplot,x,ysmooth+movav, x,movav)
        UpdateStatus(['Doing Hilbert transform...' percentdone], handles);
    end
    
end

axes(handles.PhaseAxes)
imagesc(phase)
axis image xy
title('Phase')
colorbar
axes(handles.PhasediagAxes)
imagesc(imorig)
axis image xy
title('Filtered image')
end

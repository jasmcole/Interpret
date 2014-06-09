function phase = HilbertPhaseRetrieve(handles)

data = double(handles.dataimage);
intref = double(handles.refimage);
calibdata = handles.calibdata;

[nrows ncols] = size(data);

phase = zeros(size(data));
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
    phase(n,:) = unwrap(angle(Hy) - angleHyref);
    imorig(n,:) = y;
    imref(n,:) = yref;
    count = count + 1;
    if (count > nrows/10)
        count = 0;
        percentdone = [num2str(floor(100*n/nrows)) '%'];
        axes(handles.PhasediagAxes)
        x = 1:length(y);
        plot(x,yplot,x,movav,x,ysmooth+movav)
        set(handles.StatusBox, 'String', ['Doing Hilbert transform...' percentdone]); drawnow
    end
    
end

%Remove edges to get rid of problems
phase(:,1:10) = [];
phase(:,end-9:end) = [];

phase = flipud(phase);
phase = flipud(unwrap(flipud(phase)));
phase = imrotate(phase, 90);

phase = unwrap(phase);

phase = imrotate(phase, 90);
phase = fliplr(phase);
imorig = imorig;

axes(handles.PhaseAxes)
imagesc(phase)
axis image xy
colorbar
axes(handles.PhasediagAxes)
imagesc(imorig)
axis image xy

end

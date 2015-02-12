function phase = PhaseRetrieve(handles)

data = double(handles.dataimage);
intref = double(handles.refimage);
calibdata = handles.calibdata;

fftim = fft2(data);
fftref = fft2(intref);

if (max([calibdata.xfft calibdata.yfft calibdata.wfft calibdata.hfft]) == 0)
    axes(handles.PhasediagAxes)
    set(handles.StatusBox, 'Value',0,'String','Select Fourier region'); drawnow
    imagesc(log(abs(fftshift(fftim))))
    %axis image xy
    xlim([0.4*length(fftim(1,:)) 0.6*length(fftim(1,:))])
    ylim([0.2*length(fftim(:,1)) 0.8*length(fftim(:,1))])
    rect = getrect;
    x = round(rect(1));
    y = round(rect(2));
    w = round(rect(3));
    h = round(rect(4));
    set(handles.StatusBox, 'Value',0,'String',['You picked xfft = ' num2str(x) ', yfft = ' num2str(y) ', wfft = ' num2str(w) ', hfft = ' num2str(h)]); drawnow
else
    x = calibdata.xfft;
    y = calibdata.yfft;
    w = calibdata.wfft;
    h = calibdata.hfft;
end

fftim = fftshift(fftim);
fftref = fftshift(fftref);
newfftim = zeros(size(data));
newfftref = zeros(size(data));
newfftim(y:y+h, x:x+w) = fftim(y:y+h, x:x+w);
newfftref(y:y+h, x:x+w) = fftref(y:y+h, x:x+w);
newfftim = double(newfftim);
newfftref = double(newfftref);
newfftim = fftshift(newfftim);
newfftref = fftshift(newfftref);

datafiltered = ifft2(newfftim);
reffiltered = ifft2(newfftref);

if(sum(sum(isnan(reffiltered))) > 0)
    reffiltered = ones(size(datafiltered));
end
if (max(max(intref)) == min(min(intref)))
    reffiltered = ones(size(datafiltered));
end
phase = angle(datafiltered./reffiltered);
mag = abs(datafiltered./reffiltered);

axes(handles.PhaseAxes)
imagesc(phase); axis image xy; colorbar
axes(handles.PhasediagAxes)
imagesc(sin(phase)); axis image xy

phase = phase;

end

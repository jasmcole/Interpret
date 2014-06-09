function UpdateCalibrationDatabase(handles)

set(handles.StatusBox,'String','Saving calibration'); drawnow

outdata{1} = get(handles.MicperpixBox,'String'); 
outdata{2} = get(handles.RotationBox,'String');
outdata{3} = get(handles.xBox,'String');
outdata{4} = get(handles.yBox,'String'); 
outdata{5} = get(handles.wBox,'String'); 
outdata{6} = get(handles.hBox,'String'); 
outdata{7} = get(handles.xfftBox,'String'); 
outdata{8} = get(handles.yfftBox,'String'); 
outdata{9} = get(handles.wfftBox,'String'); 
outdata{10} = get(handles.hfftBox,'String'); 
outdata{11} = handles.calibdata.reference;
outdata{12} = get(handles.hsmooth1Box,'String'); 
outdata{13} = get(handles.hsmooth2Box,'String');
outdata{14} = get(handles.asmoothBox,'String');
outdata{15} = get(handles.ymidBox,'String'); 


calibration = get(handles.CalibBox,'String');

I = importdata('CalibrationDatabase.xls'); 

[nfields ncols] = size(I.textdata);

column = 0;
number = '2';

for n = 2:ncols
    if (strcmp(calibration,I.textdata{1,n}) == 1)
        column = n;
    end
end

if (column == 0)
    outdata = [calibration, outdata];
    number = '1';
    column = ncols + 1;
end

alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

if (column > 26)
   letter1 = floor(column/26);
   letter2 = column - 26*letter1;
   letter1 = alphabet(letter1);
   letter2 = alphabet(letter2);
   letter = [letter1 letter2];
end

if (column < 27)
   letter = alphabet(column); 
end

startcell = [letter number];
endcell = [letter num2str(nfields)];
range = [startcell ':' endcell];

xlwrite('CalibrationDatabase.xls',outdata','Sheet1', startcell);

set(handles.StatusBox,'String','Calibration saved'); drawnow

end
function UpdateCalibrationDatabase(handles)

copyfile('CalibrationDatabase.csv', 'CalibrationDatabase_backup.csv')

set(handles.StatusBox,'String','Saving calibration'); drawnow

newdata = get(handles.CalibrationTable, 'Data');
assignin('base','newdata',newdata)
calibration = newdata{1,1};
if(strcmp(calibration, 'Default'))
    set(handles.StatusBox,'String','Cannot overwrite defaults - save new calibration'); drawnow
    pause(2)
    return;
end
assignin('base','calibration',calibration)

data = read_mixed_csv('CalibrationDatabase.csv', ',');
[nparams nrecords] = size(data);

found = 0;
i = 0;

while ~found
    i = i + 1;
    if(strcmp(calibration, data{1,i}))
        found = 1;
    end
end

if (found > 0)
    for n = 1:nparams
        data{n,i} = newdata{n,1};
    end
    cell2csv('CalibrationDatabase.csv', data, ',')
end

set(handles.StatusBox,'String','Calibration saved'); drawnow

end
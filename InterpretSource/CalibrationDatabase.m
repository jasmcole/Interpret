function calibdata = CalibrationDatabase(calibration)

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
    
    for n = 1:nparams
        if(isempty(str2num(data{n,i})))
            calibdata.(data{n,1}) = data{n,i};
        else
            calibdata.(data{n,1}) = str2num(data{n,i});
        end
    end
            
end
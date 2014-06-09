function calibdata = CalibrationDatabase(calibration)

I = importdata('CalibrationDatabase.xls');

[nfields ncols] = size(I.textdata);

calibdb(1:ncols -1) = struct;

for n = 1:nfields
    fieldnames{n} = I.textdata{n,1};
end

%Skip first column - this contains names of fields
for n = 2:ncols
    for m = 1:length(fieldnames)
        
        if (isempty(I.textdata{m,n}))
            calibdb(n-1).(fieldnames{m}) = I.data(m,n-1);
        else
            calibdb(n-1).(fieldnames{m}) = I.textdata{m,n};
        end
        
    end
end

%Assume name field is 'Name'
for n = 1:ncols - 1
    if (strcmp(calibration, calibdb(n).Name) == 1)
        for m = 1:nfields
            if (isstr(calibdb(n).(fieldnames{m})))
                if (~isempty(str2num(calibdb(n).(fieldnames{m}))))
                    calibdata.(fieldnames{m}) = str2num(calibdb(n).(fieldnames{m}))
                else
                    calibdata.(fieldnames{m}) = calibdb(n).(fieldnames{m})
                end
            else
                calibdata.(fieldnames{m}) = calibdb(n).(fieldnames{m})
            end
        end
    end
end

calibdata.fieldnames = fieldnames;

end
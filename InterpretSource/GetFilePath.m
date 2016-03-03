function filepath = GetFilePath(fileexpression, year, month, day, run, shot)

newfile = fileexpression;
%Year/Month/Day/Run/Shot/Shot with no zeros/wildcard
%YMDRS^*
letters = 'YMDRS';

for s = 1:length(letters)
    for m = 2:4
        str = repmat(letters(s),1,m);
        inds = strfind(fileexpression, str);
        for n = 1:length(inds)
            try
                formatstring = ['%0' num2str(m) 'd'];
                switch letters(s)
                    case 'Y'
                        newfile(inds(n):inds(n)+m-1) = num2str(year, formatstring);
                    case 'M'
                        newfile(inds(n):inds(n)+m-1) = num2str(month, formatstring);
                    case 'D'
                        newfile(inds(n):inds(n)+m-1) = num2str(day, formatstring);
                    case 'R'
                        newfile(inds(n):inds(n)+m-1) = num2str(run, formatstring);
                    case 'S'
                        newfile(inds(n):inds(n)+m-1) = num2str(shot, formatstring);
                end
            end
        end
    end
end

filepath = fileparts(newfile);

end
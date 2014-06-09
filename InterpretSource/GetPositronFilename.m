function fname = GetPositronFilename(month, day, run, shot, diagnostic)

    prefixdir  = '/Volumes/Drobo/Experimental Data/2013GeminiSarri/';

    datefolder = ['2013' num2str(month, '%02d') num2str(day, '%02d')];
    runfolder  = [datefolder 'r' num2str(run, '%03d')];
    shotprefix = [runfolder '_s' num2str(shot, '%03d') '_F20_' diagnostic];
    fname      = [prefixdir datefolder '/' runfolder '/' shotprefix];

end
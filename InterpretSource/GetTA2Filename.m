function fname = GetTA2Filename(month, day, run, shot, diagnostic)
    
    prefixdir  = '/Volumes/Drobo/Experimental Data/TA2 Data Complete/';
    datefolder = ['2012' num2str(month, '%02d') num2str(day, '%02d')];
    runfolder  = [datefolder 'r' num2str(run, '%03d')];
    shotprefix = [runfolder '_s' num2str(shot, '%03d') '_' diagnostic];
    fname      = [prefixdir datefolder '/' runfolder '/' shotprefix];

end
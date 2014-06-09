function fname = GetGeminiFilename(month, day, run, shot, diagnostic)

prefixdir  = '/Volumes/Drobo/Experimental Data/2012_Oct_Gemini/ServerData/';

datefolder = ['2012' num2str(month) num2str(day, '%02d')];
runfolder  = [datefolder 'r' num2str(run, '%03d')];
shotprefix = [runfolder '_s' num2str(shot, '%03d') '_F20_' diagnostic];
fname      = [prefixdir datefolder '/' runfolder '/' shotprefix];

end
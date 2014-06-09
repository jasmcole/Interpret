shots = 1:1;

for n = 1:length(shots)
    datafile = ['/Volumes/G2012_BK1/2012Gemini/ServerData/20121106/20121106r001/20121106r001_s0' num2str(shots(n), '%02d') '_F20_Interferometer.fit' ];
    q = AnalyseInterferogram(datafile);
end
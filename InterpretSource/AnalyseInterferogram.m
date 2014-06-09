function output = AnalyseInterferogram(datafile)

RetrievalType = 'Hilbert';
RhoRetrieval  = 'Abel';
GoldsteinFlag = 0;

lambda        = 800;
medium        = 'Plasma';
calib         = 'Gemini Positron 20130207'; 
rhoplotflag   = 1;

disp(['Image calibration ' calib])

calibdata = CalibrationDatabase(calib)

intref = calibdata.reference;

close all

switch RetrievalType
    case 'FFT'
        phase = PhaseRetrieve(datafile, intref, GoldsteinFlag, calibdata);
    case 'CWT'
        phase = CWTPhaseRetrieve(datafile, intref, GoldsteinFlag, calibdata);
    case 'Hilbert'
        phase = HilbertPhaseRetrieve(datafile, intref, GoldsteinFlag, calibdata);
end

switch RhoRetrieval
    case 'Abel'
        rho = AbelInversion(phase, calibdata, medium, lambda, rhoplotflag);
    case 'Fit'
        rho = FitToPlateau(phase, calibdata, medium, lambda, rhoplotflag);
end

output.rho = rho;
output.phase = phase;

[y x] = size(rho);
x = 1:x;
y = 1:y;
y = y - mean(y);
x = x*calibdata.micperpix/1000;
y = y*calibdata.micperpix/1000;
output.x = x;
output.y = y;


peakdensity = GetPeakDensity(rho);

end
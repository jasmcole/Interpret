function peakdensity = GetPeakDensity(rho)

    [ylen xlen] = size(rho);
    
    yavradius = 30;
    
    for n = 1:xlen
        peakdensity(n) = mean(rho(ylen/2 - yavradius:ylen/2 + yavradius, n));
    end
    
    mean(peakdensity)/1e6
    std(peakdensity)/1e6
    
    subplot(2,2,4)
    plot(peakdensity/1e6)
    line([0 xlen], [mean(peakdensity)/1e6 mean(peakdensity)/1e6])
    line([0 xlen], [mean(peakdensity)/1e6 + std(peakdensity)/1e6 mean(peakdensity)/1e6 + std(peakdensity)/1e6])
    line([0 xlen], [mean(peakdensity)/1e6 - std(peakdensity)/1e6 mean(peakdensity)/1e6 - std(peakdensity)/1e6])
    drawnow
    
    

end
function I2 = RemoveHotPixels(I, diam, nstd)

%This function removes pixels that are 'hotter' than the background image
%It works by first finding a mean image, that is where every pixel becomes
%the average of those around it, in a square region of side 'diam'
%The total standard deviation of pixels in the image is calculated, and if
%any pixel is more than nstd deviations away from its local mean then it is
%set to the local mean
%Raising diam gives a more smoothed image, changing nstd changes the
%threshold for hot pixel removal
%
%Jason Cole 2014 j.cole11@imperial.ac.uk

I2 = I;

[ysize xsize] = size(I);

filter = ones(diam, diam)/(diam^2);
Imean = conv2(I,filter,'same');
Istd = std(I(:));

I2(I > Imean + nstd*Istd) = Imean(I > Imean + nstd*Istd);

end
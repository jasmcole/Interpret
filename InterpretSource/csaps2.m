function I2 = csaps2(I, smooth)

[ysize xsize] = size(I);
xaxis = 1:xsize;
yaxis = 1:ysize;
I2 = csaps({yaxis, xaxis}, I, smooth, {yaxis, xaxis});

end
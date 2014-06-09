function [data intref] = SmoothInterferogram(data, intref)

[nrows ncols] = size(data);

phase = zeros(size(data));
xaxis = 1:ncols;

for n = 1:nrows
    y = data(n,:);
    yref = intref(n,:);
    movav = csaps(xaxis, y, 1e-4, xaxis);
    movavref = csaps(xaxis, yref, 1e-4, xaxis);
    y = y - movav;
    yref = yref - movavref;
    y = csaps(xaxis, y, 1e-4, xaxis);
    yref = csaps(xaxis, yref, 1e-4, xaxis);
    Hy = hilbert(y);
    Hyref = hilbert(yref);
    y = y./abs(Hy);
    yref = yref./abs(Hyref);
    data(n,:) = y;
    intref(n,:) = yref;
    n/nrows
end

end
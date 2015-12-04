function phase = VolkovUnwrap(handles)

I = handles.phase;

I = [I fliplr(I); flipud(I) flipud(fliplr(I))];

I = mod(I,2*pi);

[Ny, Nx] = size(I);

Z = exp(1i*I);

qx = ([0:Nx/2 (-Nx/2+1):-1] + 0.01)*2*pi/Nx;
qy = ([0:Ny/2 (-Ny/2+1):-1] + 0.01)*2*pi/Ny;

[qxg, qyg] = meshgrid(qx,qy);
[dIdx, dIdy] = gradient(I);
[dZdx, dZdy] = gradient(Z);

dkdy = (real(dZdy./(1i*Z)) - dIdy)/(2*pi);
dkdx = (real(dZdx./(1i*Z)) - dIdx)/(2*pi);

k = real( (1/(1i))*(ifft2((fft2(dkdx).*qxg + fft2(dkdy).*qyg)./(qxg.^2 + qyg.^2)) ));
k = k - min(k(:));
k = round(k);
k = k(1:end/2, 1:end/2);
I = I(1:end/2, 1:end/2);

phase = I + 2*pi*k;
phase = phase - min(phase(:));

end
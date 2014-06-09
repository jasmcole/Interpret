% Assume Matlab dimensions, y is first index, x is second

function fnew = ReconstructFromDerivativeTest(f, x, y)

dx = x(2) - x(1);
dy = y(2) - y(1);

Nx = length(x);
Ny = length(y);

qxmax = 2*pi/dx;
qymax = 2*pi/dy;

qx = linspace(0,qxmax,Nx);
qy = linspace(0,qymax,Ny);

dfdx = zeros(Ny,Nx);
dfdy = zeros(Ny,Nx);

for n = 1:length(qx)
    for m = 1:length(qy)
        qxg(m,n) = qx(n);
        qyg(m,n) = qy(m);
        oneoverq2(m,n) = 1/(qxg(m,n)^2 + qyg(m,n)^2);
        
        if ((n > 1) && (m > 1) && (n < Nx) && (m < Ny))
            dfdx(m,n) = (f(m,n+1) - f(m,n-1))/(2*dx);
            dfdy(m,n) = (f(m+1,n) - f(m-1,n))/(2*dy);
        end
    end
end
   
oneoverq2(1,1) = 0;
fnew = real((1/(2*pi*1i))*ifft2( (fft2(dfdx).*qxg + fft2(dfdy).*qyg).*oneoverq2 ));


% [qyg, qxg] = meshgrid(qy, qx);
% 
% oneoverq2 = 1./(qxg.^2 + qyg.^2);
% oneoverq2(1,1) = 0;
% 
% dfdx = gradient(f, dx);
% dfdy = gradient(f' , dy)';
% 
% fnew = real((1/(2*pi*1i))*ifft2( (fft2(dfdx).*qxg + fft2(dfdy).*qyg).*oneoverq2 ));

end
function U = Fresnel_propagate_CT(U0,lambda,z,x1)

% Written by Peter
% Idea: Use the Fresnel-Diffraction integral written as a convolution, in
% conjunction with the Fourier-Convolution Theorem.
% see https://en.wikipedia.org/wiki/Fresnel_diffraction
k=2*pi/lambda;	%Wavenumber [1/um]
Nx1=length(x1);
x1range=max(x1)-min(x1);
y1=x1;
[X1,Y1]=meshgrid(x1,y1);

%Fresnel number (rule of thumb: if FN<0.01, screen is in the far field)
FN=x1range^2/lambda/z;
disp(['Fresnel number = ' num2str(FN) ' (if <0.01 --> screen is in the far field)']);

h=1/lambda/z*exp(1i*(k*z+k/2/z*(X1.^2+Y1.^2)));
U=ifftshift(ifft2(fft2(U0).*fft2(h)));
end
function omegaout=FourierAxis(direct_axis)

% this function automatically generates x-axis for the Fourier space,
% given the x-axis in direct space. (e.g.: direct space --> timedelay; Fourier space --> omega)
% axis will be aligned to the output of FFT, which means first element will
% be 0 frequency, and not -Nyquist.

N=length(direct_axis);
isodd=logical(mod(N,2));
direct_step=direct_axis(2)-direct_axis(1);
direct_range=max(direct_axis)-min(direct_axis);
% omega_max=1/2*2*pi*1/direct_step;
% omega_min=-omega_max;
omega_step=2*pi/direct_range;
if isodd
    omegaout=[0:omega_step:(N-1)/2*omega_step -(N-1)/2*omega_step:omega_step:-omega_step];
else
    omegaout=[0:omega_step:(N/2-1)*omega_step -N/2*omega_step:omega_step:-omega_step];
end
end
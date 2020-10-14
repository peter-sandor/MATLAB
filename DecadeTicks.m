function tickarray = DecadeTicks(A)
% This function generates arrays containing the Xtick positions for line
% plots. The positions are placed equidistantly, such that they are integer
% multiples of some power of 10. It is also a requirement that there are neither not too many, nor not too few tickmarks.
% (This makes the plot look cleaner.)

Ntickmax=5;
Ntickmin=3;
corrfactor=5;
Amin=min(A);
Amax=max(A);
decade_Amin=floor(log10(abs(Amin)));
decade_Amax=floor(log10(abs(Amax)));
% decade_diff=decade_Amax-decade_Amin;
decade_diff=floor(log10(abs(Amax-Amin)));
decade_max=max(decade_Amin,decade_Amax);

if decade_Amin==decade_Amax
    tickstep=10^(decade_diff);
    tickmin=floor(sigdig(Amin,(decade_Amin-decade_diff)+1)*10^(decade_Amin-decade_diff))*10^(decade_diff);
    tickmax=ceil(sigdig(Amax,(decade_Amax-decade_diff)+1)*10^(decade_Amax-decade_diff))*10^(decade_diff);
else
    tickstep=10^(decade_max);
    tickmin=floor(sigdig(Amin,1)*10^(decade_Amin-decade_max))*10^(decade_Amin);
    tickmax=ceil(sigdig(Amax,1)*10^(decade_Amax-decade_max))*10^(decade_Amax);
end
Ntick=(tickmax-tickmin)/tickstep;
if round(Ntick)<Ntickmin
    tickstep=10^floor(log10(Amax-Amin)-1)*corrfactor;
elseif Ntick>Ntickmax
    tickstep=tickstep*corrfactor;
end
tickarray=tickmin:tickstep:tickmax;
end


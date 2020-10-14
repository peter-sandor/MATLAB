function intfactor = amp2int(ampin)
% AOM amplitude to intensity calibration based on the data and fit that can
% be found in:
% [l:\2011 09 27 E\diodeAmplitude\intensity_calibration.mat]
if iscell(ampin)==1
    intfactor=-1.6103*ampin.^3+2.3518*ampin.^2+0.25858*ampin;
else
    intfactor=-1.6103*ampin.^3+2.3518*ampin.^2+0.25858*ampin;
end
end
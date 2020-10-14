function [energies,psiSquared]=PEspectrum(s_ion,ladderind,eMax,varargin)

% This function calculates the PhotoeElectron spectrum from the output of
% the 0D StrongFieldIonization code based on the population of the chosen
% ladder states.

if nargin==3
    timeind=size(s_ion,1);
elseif nargin==4
    timeind=varargin{1};
else
    disp('invalid number of arguments.')
    energies=[];
    psiSquared=[];
    return;
end
ionPlot = s_ion(timeind,ladderind(1):ladderind(2));
eMaxReduced = eMax/(1.602e-19);
resolution = 0.01; %in eV
energies=(0:resolution:eMaxReduced); 
maxL = size(ionPlot,2);
runningSum = zeros(size(energies));
for l=1:maxL
    legendAll = legendre(l-1,2*energies/eMaxReduced-1);
    polynomial = legendAll(1,:);
    runningSum(:) = ionPlot(l)*sqrt((2*l+1)/eMaxReduced)*polynomial' + runningSum(:);
end
energies=energies-eMaxReduced/2;
psiSquared = abs(runningSum).^2;
end
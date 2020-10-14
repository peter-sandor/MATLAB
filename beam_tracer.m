function output=beam_tracer(input,system)

% This code calculates the parameters of a gaussian beam after propagation
% of a series of lenses, separated by free space. The geometry of the
% optical setup and the input beam parameters are to be specified in the
% structure "input". Reference for the calculations:
% Saleh B.E.A., Teich M.C. Fundamentals of photonics (Wiley, 1991)
% page 80 and onward

% system: an array of structs, where each struct has two fields:
%       system(n).element: a string,either 'lns' for thin lens or 'fsp'
%           for freespace propagation, etc.
%       system(n).param: a numerical array, with the relevant parameters for that particular element (focal length in [mm] for lens, distance also in [mm] for freespace propagation, etc.)

% the structure "input" should contain the following fields:
% input.W: initial beam radius (in millimeters)
% input.R: initial wavefront curvature (in millimeters)
% input.lambda: operating wavelength (in millimeters)

N=length(system);
q0=1/(1/input.R-1i*input.lambda/(pi*input.w^2));
vec=[q0;1];
output.lambda=input.lambda;

for indS=1:N
    if system(indS).element=='lns'    % focusing by thin lens
        vec=[1,0;-1/system(indS).param,1]*vec;
    elseif system(indS).element=='fsp'    % propagation in free space
        vec=[1,system(indS).param;0,1]*vec;
    elseif system(indS).element=='rff'    % refraction at flat interface; system(k).param(1:2): refractive index of first and second medium, respectively.
        vec=[1,0;0,system(indS).param(1)/system(indS).param(2)]*vec;
        output.lambda=system(indS).param(1)/system(indS).param(2)*output.lambda;
    elseif system(indS).element=='rfc'    % refraction at curved interface; system(k).param(1:3): refractive index of first and second medium, and radius of curvature, respectively.
        vec=[1,0;(system(indS).param(1)-system(indS).param(2))/(system(indS).param(2)*system(indS).param(3)),system(indS).param(1)/system(indS).param(2)]*vec;
        output.lambda=system(indS).param(1)/system(indS).param(2)*output.lambda;
    else 'unknown element'
    end
    vec=vec/vec(2);
end

output.R=1/real(1/vec(1)); % output wavefront curvature
output.w=sqrt(-1/imag(1/vec(1))*output.lambda/pi); % output beam radius
output.BFL=output.R/(1+(output.lambda*output.R/(pi*output.w^2))^2); % back focal length of the system
output.beam_waist=output.w./sqrt(1+(pi*output.w^2./(output.lambda*output.R)).^2); % calculate beam waist (minimum radius) using the output wavefront curvature and beam radius
output.z_waist=output.R./(1+(output.lambda*output.R./(pi*output.w^2)).^2); % calculate beam waist (minimum radius) using the output wavefront curvature and beam radius
end
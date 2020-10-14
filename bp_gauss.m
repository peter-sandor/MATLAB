function output=bp_gauss(input)

% This code calculates the parameters of a gaussian beam after propagation
% of a series of lenses, separated by free space. The geometry of the
% optical setup and the input beam parameters are to be specified in the
% structure "input". Reference for the calculations:
% Saleh B.E.A., Teich M.C. Fundamentals of photonics (Wiley, 1991)
% page 80 and onward

% the structure "input" should contain the following fields:
%
% input.element: an array of strings, each string being either 'lns' for
%               lens or 'fsp' for freespace propagation
% input.s: an array of cells, where each cell contains an array with the relevant parameters for the element of the same array index (focal length in [mm] for lens, distance also in [mm] for freespace propagation, etc.)
% input.W: initial beam radius (in millimeters)
% input.R: initial wavefront curvature (in millimeters)
% input.lambda: operating wavelength (in millimeters)

N=length(input.s);
q0=1/(1/input.R-i*input.lambda/(pi*input.w^2));
vec=[q0;1];
output.lambda=input.lambda;

for k=1:N
    if input.element(k,:)=='lns'    % focusing by thin lens
        vec=[1,0;-1/input.s{k}(1),1]*vec;
    elseif input.element(k,:)=='fsp'    % propagation in free space
        vec=[1,input.s{k}(1);0,1]*vec;
    elseif input.element(k,:)=='rff'    % refraction at flat interface; input.s{}(1:2): refractive index of first and second medium, respectively.
        vec=[1,0;0,input.s{k}(1)/input.s{k}(2)]*vec;
        output.lambda=input.s{k}(1)/input.s{k}(2)*output.lambda;
    elseif input.element(k,:)=='rfc'    % refraction at curved interface; input.s{}(1:3): refractive index of first and second medium, and radius of curvature, respectively.
        vec=[1,0;(input.s{k}(1)-input.s{k}(2))/(input.s{k}(2)*input.s{k}(3)),input.s{k}(1)/input.s{k}(2)]*vec;
        output.lambda=input.s{k}(1)/input.s{k}(2)*output.lambda;
    else 'unknown element'
    end
    vec=vec/vec(2);
end

output.R=1/real(1/vec(1)); % output wavefront curvature
output.w=sqrt(-1/imag(1/vec(1))*output.lambda/pi); % output beam radius
output.BFL=output.R/(1+(output.lambda*output.R/(pi*output.w^2))^2); % back focal length of the system
output.beam_waist=output.w./sqrt(1+(pi*output.w^2./(output.lambda*output.R)).^2); % calculate beam waist (minimum radius) using the output wavefront curvature and beam radius
end
function out = LegendreP(x,Jmax,varargin)
% This function generates the associated Legendre functions up to order
% Jmax.
% Input:    x: double array of size (N,1) or (1,N), -1<=x<=1
%           Jmax: integer value, Jmax>=0
%           optional input: 'norm' (string) to normalize the functions
% Output:   out = double array of size (N,Jmax+1,2*Jmax+1)
%                   first index: x, second index: L (0 to Jmax), third
%                   index: m (-Jmax to Jmax)) /m = 0 is at index Jmax+1/

if strcmp(varargin,'norm')
    flag_norm = 1;
else
    flag_norm = 0;
end
if ~isempty(vec2ind(x<-1)) || ~isempty(vec2ind(x>1))
    disp('x has to be: -1<=x<=1 !')
    out=[];
    return;
end
x=map2colvec(x);
Nx=length(x);
out = zeros([Nx, Jmax+1 2*Jmax+1]);
out(:,1,Jmax+1)=1;
if Jmax>0
    J = 1;
    out(:,J+1,Jmax+1) = x;
    out(:,J+1,Jmax+1+J) = -sqrt(1-x.^2);
    out(:,J+1,Jmax+1-J) = sqrt(1-x.^2)/2;
    while J<Jmax
        J=J+1;
        % calculate positive M-branch
        out(:,J+1,Jmax+1+J) = -(2*J-1)*sqrt(1-x.^2).*out(:,J,Jmax+J);
        out(:,J+1,Jmax+J) = (2*J-1)*x.*out(:,J,Jmax+J);
        
        for M=0:J-2
            out(:,J+1,Jmax+1+M) = 1/(J-M)*((2*J-1)*x.*out(:,J,Jmax+1+M) - (J+M-1)*out(:,J-1,Jmax+1+M));
        end
        % calculate negative M-branch
        for M=1:J
            factorM=1/prod([J-M+1:1:J+M]);
            out(:,J+1,Jmax+1-M) = (-1)^M*factorM*out(:,J+1,Jmax+1+M);
        end
    end
end
if flag_norm
    for ind1 = 1:Jmax+1
        out(:,ind1,:) = out(:,ind1,:)*sqrt((2*(ind1-1)+1)/2);
    end
end
end
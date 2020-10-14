function [U,temp_min,temp_max] = Fresnel_propagate_basic(U0,lambda,z,x1,x2)

% Written by Spencer, modified by Peter

k=2*pi/lambda;	%Wavenumber [1/um]

Nx1=length(x1);
x1range=max(x1)-min(x1);
x1step=x1range/(Nx1-1);
y1step=x1step;
% x1min=-x1range/2; % [micron]
% x1max=-x1min; % [micron]
% x1step=2*x1max/(Nx1-1);
% x1=x1min:x1step:x1max;
% y1step=x1step;
y1=x1;
[X1,Y1]=meshgrid(x1,y1);

Nx2=length(x2);
x2range=max(x2)-min(x2);
x2step=x2range/(Nx2-1);
% x2min=-x2range/2; % [micron]
% x2max=-x2min;
% x2step=2*x2max/(Nx2-1);
% x2=x2min:x2step:x2max;
y2=x2;
% [X2,Y2]=meshgrid(x2,y2);

 %Fresnel number (rule of thumb: if FN<0.01, screen is in the far field)
FN=(x1range*x2range)^2/lambda/z;
disp(['Fresnel number = ' num2str(FN) ' (if <0.01 --> screen is in the far field)']);
for ind1=1:Nx2
    for ind2=1:Nx2
        U(ind1,ind2)=-1i.*exp(1i.*k.*z)./(lambda.*z).*sum(sum(U0.*exp(1i.*k.*((x2(ind1)-X1).^2+(y2(ind2)-Y1).^2)./(2.*z)),1),2).*x1step.*y1step;
        temp_min(ind1,ind2)=min(min((x2(ind1)-X1).^2+(y2(ind2)-Y1).^2))/(10*z^3*lambda);
        temp_max(ind1,ind2)=max(max((x2(ind1)-X1).^2+(y2(ind2)-Y1).^2))/(10*z^3*lambda); % this
%     value has to be <<1 for the Fresnel approximation to be valid
    end
    if mod(ind1,20)==0
        disp([num2str(ind1) '/' num2str(Nx2)]);
    end
end
end
% clear;
N2=180;
x_min=0;
x_max=99;
x_step=(x_max-x_min)/(N2-1);
x=(x_min:x_step:x_max)';

eps=x_step/10;
eps2=1e-4;
sigma=20;
omega=0.5;
a=39;

noise_factor=0;
id=1;
switch id
    case 1
        func.name='step';
        func.projd=zeros(size(x));
        func.projd(abs(x)<=a)=2*sqrt(a^2-x(abs(x)<=a).^2);
    case 2
        func.name='gauss';
        func.projd=sigma*sqrt(pi)*exp(-x.^2/sigma^2);
    case 3
        func.name='bessel';
        func.projd=2*cos(omega*x)/omega;
    case 4
        func.name='a2-r2';
        func.projd=zeros(size(x));
        func.projd(abs(x)<=a)=4/3*(a^2-x(abs(x)<=a).^2).^(3/2);
    case 5
        x=x(abs(x)>=eps2);
        func.name='acosh';
        func.projd=zeros(size(x));
        func.projd(abs(x)<=a)=a*sqrt(a^2-x(abs(x)<=a).^2)-x(abs(x)<=a).^2.*acosh(a./x(abs(x)<=a));
    case 6
        sigma=.6;
        func.projd=exp(-(x-x(ceil(N2/2))).^2/sigma^2);
end
func1=func.projd+noise_factor*max(func.projd)*rand(size(x));% add noise

func_invtd=invAbel1(func1);
func_invtd=func_invtd/max(func_invtd);
[rho,func_invtd2]=invAbel5(x,func1',length(x));

switch id
    case 1
        func.orig=double(abs(rho)<=a);
    case 2
        func.orig=exp(-rho.^2/sigma^2);
    case 3
        func.orig=besselj(0,omega*rho);
    case 4
        func.orig=zeros(size(rho));
        func.orig(abs(rho)<=a)=(a^2-rho(abs(rho)<=a).^2);
    case 5
        func.orig=zeros(size(rho));
        func.orig(abs(rho)<=a)=a-rho(abs(rho)<=a);
    case 6
        func.orig=zeros(size(rho));
end
func2=func.orig;
[x2,func_projd]=Abel(rho,func.orig,N2);

figure;plot(rho,func_invtd2,'r')
hold on;plot(rho,func.orig/max(func.orig),'k')
% plot(rho,func_invtd,'r--')
plot(rho,func1/max(func1),'m')
plot(1:100,func_projd/max(func_projd),'b')

%%
% clear;
N2=180;
x_min=0;
x_max=99;
x_step=(x_max-x_min)/(N2-1);
x=(x_min:x_step:x_max)';

eps=x_step/10;
eps2=1e-4;
sigma=20;
omega=0.5;
a=39;
a_max=39;
id=1;
func.projd=zeros(size([N2 N2]));
for ind1=1:N2
    a(ind1)=a_max*max(cos(x(ind1)/a_max),0);
    switch id
        case 1
            func.name='step';
%             func.projd=zeros(size([N2 N2]));
            func.projd(abs(x)<=a(ind1),ind1)=2*sqrt(a(ind1)^2-x(abs(x)<=a(ind1)).^2);
        case 2
            func.name='gauss';
            func.projd=sigma*sqrt(pi)*exp(-x.^2/sigma^2);
        case 3
            func.name='bessel';
            func.projd=2*cos(omega*x)/omega;
        case 4
            func.name='a2-r2';
            func.projd=zeros(size(x));
            func.projd(abs(x)<=a)=4/3*(a^2-x(abs(x)<=a).^2).^(3/2);
        case 5
            x=x(abs(x)>=eps2);
            func.name='acosh';
            func.projd=zeros(size(x));
            func.projd(abs(x)<=a)=a*sqrt(a^2-x(abs(x)<=a).^2)-x(abs(x)<=a).^2.*acosh(a./x(abs(x)<=a));
        case 6
            sigma=.6;
            func.projd=exp(-(x-x(ceil(N2/2))).^2/sigma^2);
    end 
    
    
end



%%

flagB = 1; %lattice-centered
% flagB = 0; %pixel-centered

%set up a toolbox, containing weights distributions, which only depends on
%number of rings, 'NumRing', in a slice (x = constant).

%Given x and AOI(area of interest, which is defined by 'Radius'), 'NumRing'
%is determined by NumRing =(Radius^2-x^2)^0.5. And so is the linear maps
%between projected pixel values, 'Proj(y,x)', and ring densities, 'Ring(y,x)'.
% Proj(NumRing,y) = Area(NumRing,y,k)*Ring(NumRing,k), summed over k>=y.
% Ring(NumRing,y) = Weight(NumRing,y,k)*Proj(NumRing,k), summed over k<=y.
% 'Area' and 'Weight' are pre-define matrix-valued functions of NumRing st. Area*Weight = 1.

% Here to define 'Area' and its inverse 'Weight' in one structure array 'Dist_Area':
SizePixel = 150; %size of 'Dist_Area' in pixels
SizeFactor = 1; %refinement for better resolution
Size = SizePixel*SizeFactor;
Dim = 360*SizeFactor;

if flagB == 1    % areas first, lattice-centered
    Dist_Area = struct('Area',{},'Weight',{});
    for NumRing = 1:Size
        Dist_Area(NumRing).Area = zeros(NumRing,NumRing);
        Dist_Area(NumRing).Weight = zeros(NumRing,NumRing);
        for y = 1:NumRing
            Sum = 0;
            for k = y:NumRing
                Dist_Area(NumRing).Area(y,k) = k^2*(acos((y-1)/k)-acos(y/k))+y*(k^2-y^2)^.5-(y-1)*(k^2-(y-1)^2)^.5-Sum;
                Sum = Sum + Dist_Area(NumRing).Area(y,k);
            end
        end
        Dist_Area(NumRing).Weight = inv(Dist_Area(NumRing).Area);
    end
end
if flagB == 0     % areas first, pixel-centered
    Dist_Area = struct('Area',{},'Weight',{});
    for NumRing = 1:Size
        Dist_Area(NumRing).Area = zeros(NumRing,NumRing);
        Dist_Area(NumRing).Weight = zeros(NumRing,NumRing);
        for y = 1:NumRing
            Sum = 0;
            for k = y:NumRing
                Dist_Area(NumRing).Area(y,k) = (k-.5)^2*(acos((y-1.5)/(k-.5))-acos((y-.5)/(k-.5)))+(y-.5)*((k-.5)^2-(y-.5)^2)^.5-(y-1.5)*((k-.5)^2-(y-1.5)^2)^.5-Sum;
                Sum = Sum + Dist_Area(NumRing).Area(y,k);
            end
        end
        Dist_Area(NumRing).Area(1,:) = Dist_Area(NumRing).Area(1,:)/2; %centre ring
        Dist_Area(NumRing).Weight = inv(Dist_Area(NumRing).Area);
    end
end
Inv = Dist_Area(99).Weight*func1;
plot(rho,Inv/max(Inv),'g');
% hold on;plot(rho,func_invtd/max(func_invtd),'g')
% 
% figure;plot(x2,func_projd,'m')
% hold on;plot(x,func.projd,'b')
function COLOR_out = generate_colormap(varargin)

if nargin == 1
    Ncolor = varargin{1};
else
    Ncolor = 8;
end
if Ncolor>=8
    Nmin = floor(Ncolor^(1/3));
    Nmax = Ncolor/Nmin^2;
    
    A{1} = 0:1/(Nmin-1):1;
    A{2} = A{1};
    A{3} = 0:1/(Nmax-1):1;
    color_gray = [1/(Nmin-1)/2 1/(Nmin-1)/2 1/(Nmin-1)/2];
else
    A{1} = [0 1];
    A{2} = A{1};
    A{3} = A{1};
    color_gray = [0.5 0.5 0.5];
end

COLOR_out = combine(A);
COLOR_out(:,2:3) = circshift(COLOR_out(:,2:3),[0 1]);
COLOR_out(end,:) = color_gray; % replace white with grey
COLOR_out(Ncolor+1:end,:) = [];
end
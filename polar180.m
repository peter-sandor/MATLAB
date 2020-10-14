function hpol = polar180(varargin)
    %POLAR  Polar coordinate plot.
    %   POLAR(THETA, RHO) makes a plot using polar coordinates of
    %   the angle THETA, in radians, versus the radius RHO.
    %   POLAR(THETA, RHO, S) uses the linestyle specified in string S.
    %   See PLOT for a description of legal linestyles.
    %
    %   POLAR(AX, ...) plots into AX instead of GCA.
    %
    %   H = POLAR(...) returns a handle to the plotted object in H.
    %
    %   Example:
    %      t = 0 : .01 : 2 * pi;
    %      polar(t, sin(2 * t) .* cos(2 * t), '--r');
    %
    %   See also PLOT, LOGLOG, SEMILOGX, SEMILOGY.
    
    %   Copyright 1984-2015 MathWorks, Inc.
    
    % Parse possible Axes input
    [cax, args, nargs] = axescheck(varargin{:});
    if nargs < 1
        error(message('MATLAB:narginchk:notEnoughInputs'));
    elseif nargs > 3
        error(message('MATLAB:narginchk:tooManyInputs'));
    end
    
    if nargs < 1 || nargs > 3
        error(message('MATLAB:polar:InvalidDataInputs'));
    elseif nargs == 2
        theta = args{1};
        rho = args{2};
        if ischar(rho)
            line_style = rho;
            rho = theta;
            [mr, nr] = size(rho);
            if mr == 1
                theta = 1 : nr;
            else
                th = (1 : mr)';
                theta = th(:, ones(1, nr));
            end
        else
            line_style = 'auto';
        end
    elseif nargs == 1
        theta = args{1};
        line_style = 'auto';
        rho = theta;
        [mr, nr] = size(rho);
        if mr == 1
            theta = 1 : nr;
        else
            th = (1 : mr)';
            theta = th(:, ones(1, nr));
        end
    else % nargs == 3
        [theta, rho, line_style] = deal(args{1 : 3});
    end
    if ischar(theta) || ischar(rho)
        error(message('MATLAB:polar:InvalidInputType'));
    end
    if ~isequal(size(theta), size(rho))
        error(message('MATLAB:polar:InvalidInputDimensions'));
    end
    
    % get hold state
    cax = newplot(cax);
    
    next = lower(get(cax, 'NextPlot'));
    hold_state = ishold(cax);
    
    % get x-axis text color so grid is in same color
    % get the axis gridColor
    axColor = get(cax, 'Color');
    gridAlpha = get(cax, 'GridAlpha');
    axGridColor = get(cax,'GridColor').*gridAlpha + axColor.*(1-gridAlpha);
    tc = axGridColor;
    ls = get(cax, 'GridLineStyle');
    
    % Hold on to current Text defaults, reset them to the
    % Axes' font attributes so tick marks use them.
    fAngle = get(cax, 'DefaultTextFontAngle');
    fName = get(cax, 'DefaultTextFontName');
    fSize = get(cax, 'DefaultTextFontSize');
    fWeight = get(cax, 'DefaultTextFontWeight');
    fUnits = get(cax, 'DefaultTextUnits');
    set(cax, ...
        'DefaultTextFontAngle', get(cax, 'FontAngle'), ...
        'DefaultTextFontName', get(cax, 'FontName'), ...
        'DefaultTextFontSize', get(cax, 'FontSize'), ...
        'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
        'DefaultTextUnits', 'data');
    
    % only do grids if hold is off
    if ~hold_state
        
        % make a radial grid
        hold(cax, 'on');
        % ensure that Inf values don't enter into the limit calculation.
        arho = abs(rho);
%         maxrho = max(arho(arho ~= Inf));
        maxrho=4;
        hhh = line([-maxrho, -maxrho, maxrho, maxrho], [-maxrho, maxrho, maxrho, -maxrho], 'Parent', cax);
        set(cax, 'DataAspectRatio', [1, 1, 1], 'PlotBoxAspectRatioMode', 'auto');
        v = [get(cax, 'XLim') get(cax, 'YLim')];
        ticks = sum(get(cax, 'YTick') >= 0);
        delete(hhh);
        % check radial limits and ticks
        rmin = 0;
%         rmax = v(4);
%         rmax=1.2*max(rho);
        rmax=4;
        rticks=2;
%         rticks = max(ticks - 1, 2);
%         if rticks > 5   % see if we can reduce the number
%             if rem(rticks, 2) == 0
%                 rticks = rticks / 2;
%             elseif rem(rticks, 3) == 0
%                 rticks = rticks / 3;
%             end
%         end
        
        % define a circle
        th = 0 : pi / 50 : pi;
        xunit = -sin(th);
        yunit = cos(th);
        % now really force points on x/y axes to lie on them exactly
        inds = 1 : (length(th) - 1) / 2 : length(th);
%         xunit(inds(2 : 2 : length(inds))) = zeros(length(2 : 2 : length(inds)), 1);
%         yunit(inds(1 : 2 : length(inds))) = zeros(length(1 : 2 : length(inds)), 1);
        % plot background if necessary
        if ~ischar(get(cax, 'Color'))
            patch('XData', xunit * rmax, 'YData', yunit * rmax, ...
                'EdgeColor', tc, 'FaceColor', get(cax, 'Color'), ...
                'HandleVisibility', 'off', 'Parent', cax);
        end
        
        % draw radial circles
        rinc = (rmax - rmin) / rticks;
%          for indS = 0:1:floor(rmax)% (rmin + rinc) : rinc : rmax
%             hhh = line(xunit * indS, yunit * indS, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1, ...
%                 'HandleVisibility', 'off', 'Parent', cax);
%             text(-indS, 0, ...
%                 ['' num2str(indS)], 'VerticalAlignment', 'bottom', ...
%                 'HandleVisibility', 'off', 'Parent', cax);
%         end
%         set(hhh, 'LineStyle', '-'); % Make outer circle solid
        
        % plot spokes
        Nspoke=4;
        th = (1 : Nspoke) * pi / (Nspoke);
        cst = -sin(th);
        snt = cos(th);
        cs = [cst; zeros(size(cst))];
        sn = [snt; zeros(size(snt))];
        line(rmax * cs, rmax * sn, 'LineStyle', ls, 'Color', tc, 'LineWidth', 1, ...
            'HandleVisibility', 'off', 'Parent', cax);
        line(rmax * cs, -rmax * sn, 'LineStyle', ls, 'Color', tc, 'LineWidth', 1, ...
            'HandleVisibility', 'off', 'Parent', cax);
        % annotate spokes in degrees
        rt = 1.2 * rmax;
%         text(0, rt, '0',...
%                 'HorizontalAlignment', 'center', ...
%                 'HandleVisibility', 'off', 'Parent', cax);
        for indS = 1 : Nspoke
            text(rt * (cst(indS)), rt * snt(indS), int2str(indS * 180/Nspoke),...
                'HorizontalAlignment', 'center', ...
                'HandleVisibility', 'off', 'Parent', cax);
            if indS == length(th)
%                 loc = int2str(0);
                text(rt * 0, rt * 1, '0', 'HorizontalAlignment', 'center', ...
                'HandleVisibility', 'off', 'Parent', cax);
                text(rt * 0, rt * -1, '180', 'HorizontalAlignment', 'center', ...
                'HandleVisibility', 'off', 'Parent', cax);
            end

        end
        
        % set view to 2-D
        view(cax, 2);
        % set axis limits
        axis(cax, rmax * [-1, 1, -1.15, 1.15]);
    end
    
    % Reset defaults.
    set(cax, ...
        'DefaultTextFontAngle', fAngle , ...
        'DefaultTextFontName', fName , ...
        'DefaultTextFontSize', fSize, ...
        'DefaultTextFontWeight', fWeight, ...
        'DefaultTextUnits', fUnits );
    
    % transform data to Cartesian coordinates.
    xx = -rho .* sin(theta);
    yy = rho .* cos(theta);
    
    % plot data on top of grid
    if strcmp(line_style, 'auto')
        q = plot(xx, yy, 'Parent', cax);
    else
        q = plot(xx, yy, line_style, 'Parent', cax);
    end
    
    if nargout == 1
        hpol = q;
    end
    
    if ~hold_state
        set(cax, 'DataAspectRatio', [1, 1, 1]), axis(cax, 'off');
        set(cax, 'NextPlot', next);
    end
    set(get(cax, 'XLabel'), 'Visible', 'on');
    set(get(cax, 'YLabel'), 'Visible', 'on');
    
    if ~isempty(q) && ~isdeployed
        makemcode('RegisterHandle', cax, 'IgnoreHandle', q, 'FunctionName', 'polar');
    end
end

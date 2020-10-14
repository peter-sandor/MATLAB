function varargout = GVD_Sanjay(lambda0,material,varargin)

N=length(lambda0);
c=0.3; % [um/fs]
dlambda=0.00001;

%If trying to examine 400nm in a bbo uncomment the function below to get
%the correct index of refraction for the 400nm light. 
%n0 = IndexEllipsoid(lambda0,'BBO_o','BBO_e',29.2);
A = ismember(lambda0,0.4);
if (A(1) == 1) && (strncmp(material, 'BBO', 3) == 1)
    n0 = IndexEllipsoid(lambda0,'BBO_o','BBO_e',29.2);
else
    n0 = sellmeier(lambda0,material);
end
method = 'pchip';
if ~isempty(n0)
    n_minus = interp1(lambda0,n0,lambda0-dlambda/2,method,'extrap');
    n_plus = interp1(lambda0,n0,lambda0+dlambda/2,method,'extrap');
    dn_dlambda=(n_plus-n_minus)./dlambda;
    n_group = n0 - lambda0 .* dn_dlambda;
    v_group =[map2colvec(c ./ n_group)];
    
    dn_minus = interp1(lambda0,dn_dlambda,lambda0-dlambda/2,method,'extrap');
    dn_plus = interp1(lambda0,dn_dlambda,lambda0+dlambda/2,method,'extrap');
    d2n_dlambda2 = (dn_plus-dn_minus)./dlambda;
	GVD = map2colvec(((map2colvec(lambda0).^3)/(2*pi*c^2)).*map2colvec(d2n_dlambda2)); %[fs^2/micron/rad]
    
    d2n_minus = interp1(lambda0,d2n_dlambda2,lambda0-dlambda/2,method,'extrap');
    d2n_plus = interp1(lambda0,d2n_dlambda2,lambda0+dlambda/2,method,'extrap');
    d3n_dlambda3 = (d2n_plus-d2n_minus)./dlambda;
	TOD=[-map2colvec(map2colvec(lambda0).^4/4/pi^2/c^3.*(3*(map2colvec(d2n_dlambda2))+map2colvec(lambda0.*d3n_dlambda3)))]; %[fs^3/micron]
    
    % Phi(omega)=Phi0 + Phi1*(omega-omega0) +  Phi2*(omega-omega0)^2 +
    % Phi3*(omega-omega0)^3 + ...
    % here Phi2=GVD/2, Phi3=TOD/6

    varargout{3}=TOD;
    varargout{2}=GVD;
    varargout{1}=v_group;
    
    if nargin==3 && strcmp(varargin{1},'plot')
        figure;
        subplot(311)
        plot(lambda0,v_group,'k')
        title(['material: ' material])
        % xlabel('lambda [\mum]')
        ylabel('group velocity [\mum/fs]')
        grid on;

        subplot(312)
        plot(lambda0,GVD*1e3,'k')
        % xlabel('lambda [\mum]')
        ylabel('GVD [fs^2/mm/rad]')
        grid on;

        subplot(313)
        plot(lambda0,TOD*1e3,'k')
        xlabel('lambda [\mum]')
        ylabel('TOD [fs^3/mm/rad^2]')
        grid on;
    end
end
    
end
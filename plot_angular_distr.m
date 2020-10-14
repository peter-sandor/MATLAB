function varargout=plot_angular_distr(varargin)
% coeffs has to be an Nx1 vector if only the angular dependence is to be
% calculated. It should be an Nx2 matrix with the uncertainties in the
% coefficients in the second column if the histogram is to be plotted.
coeffs=varargin{1};
basetype=varargin{2};
plottype=varargin{3};
if nargin==4
    hax=varargin{4};
else
    figure;hax(1)=axes;
    figure;hax(2)=axes;
end
Ncoeff=size(coeffs,1);
term_index=1:Ncoeff;

basisnames{1}='1';
if strcmp(basetype,'sin')
    basisnames{2}='sin(\theta)^2';
    basisnames{3}='sin(2*\theta)^2';
    basisnames{4}='sin(3*\theta)^2';
    basisnames{5}='sin(4*\theta)^2';
    coeff_int=(cos(2*pi*(term_index-1))-1+8*(term_index-1).^2)./(8*(term_index-1).^2-2);
%     coeff_int=[2; 4/3; 16/15; 36/35];
elseif strcmp(basetype,'cos')
    basisnames{2}='cos(\theta)^2';
    basisnames{3}='cos(2*\theta)^2';
    basisnames{4}='cos(3*\theta)^2';
    basisnames{5}='cos(4*\theta)^2';
    coeff_int=-(cos(2*pi*(term_index-1))+3-8*(term_index-1).^2)./(8*(term_index-1).^2-2);
%     coeff_int=[2; 2/3; 14/15; 34/35];
else
    disp('invalid basis function ID!');
    varargout{1}=[];
    varargout{2}=[];
    varargout{3}=[];
    return;
end
coeff_int(1)=2;
flaghist=size(coeffs,2);
if flaghist==2
    coeffs_sigma=coeffs(:,2);
    coeffs=coeffs(:,1); 
end
Ntheta=(Ncoeff-1)*200 + 1;
theta=map2colvec(0:2*pi/(Ntheta-1):2*pi);

normfactor=0;
for ind1=1:Ncoeff
    normfactor=normfactor+coeffs(ind1)*coeff_int(ind1)/2/pi;
end
coeffs=coeffs/normfactor; % normalize coefficients so that the angular distribution has unit area
if flaghist==2
    coeffs_sigma=coeffs_sigma/normfactor;
end

rho=coeffs(1);
if strcmp(basetype,'sin')
    for ind1=2:Ncoeff
        newterm=eval(['sin(' num2str((ind1-1)) '*theta).^2']);
        rho=rho+coeffs(ind1)*newterm;
    end
elseif strcmp(basetype,'cos')
    for ind1=2:Ncoeff
        newterm=eval(['cos(' num2str((ind1-1)) '*theta).^2']);
        rho=rho+coeffs(ind1)*newterm;
    end
end
% normfactor2=trapz(theta(1:floor(Ntheta/2)+1),sin(theta(1:floor(Ntheta/2)+1)).*rho(1:floor(Ntheta/2)+1))/2/pi
% rho=rho/normfactor; % normalize angular distribution to have unit area
varargout{1}=theta;
varargout{2}=rho;
varargout{3}=normfactor;

if strcmp(plottype,'newplot')
    str_title=[];
    for ind1=1:length(coeffs)
        str_title=[str_title '(' num2str(roundP(coeffs(ind1),3)) '\pm' num2str(roundP(coeffs_sigma(ind1),3)) ')*' basisnames{ind1} ' + '];
    end
    str_title=str_title(1:end-3);
    
    pax=polarP(hax(1),theta,rho,'r.-');
    set(pax,'linewidth',1);
    title(hax(1),str_title)
elseif strcmp(plottype,'addplot')
    hold on;
    [x,y]=pol2cart(hax(1),theta,rho);
    line(y,x,'linewidth',2);
    hold off;
end
%%
if flaghist==2
    Ns=5000;
    Nr=100;
    Nsigma=2.4;
    ax_limit=1.4*max(rho);
    coeffs2=zeros([Ncoeff,Ns]);
    rho_ens=zeros([Ntheta Ns]);
    for ind1=1:Ns
        rho_ens(:,ind1)=zeros([Ntheta 1]);
        for ind2=1:Ncoeff
            coeffs2(ind2,ind1) = randn*coeffs_sigma(ind2) + coeffs(ind2);
            if ind2 == 1
                rho_ens(:,ind1) = rho_ens(:,ind1) + coeffs2(ind2,ind1)*ones([Ntheta 1]);
            else
                if strcmp(basetype,'sin')
                    rho_ens(:,ind1) = rho_ens(:,ind1) + map2colvec(coeffs2(ind2,ind1)*sin((ind2-1)*theta).^2);
                elseif strcmp(basetype,'cos')
                    rho_ens(:,ind1) = rho_ens(:,ind1) + map2colvec(coeffs2(ind2,ind1)*cos((ind2-1)*theta).^2);
                end
            end
        end
    end

    [Rs,phist]=polar_binner(theta,rho_ens,Nr);
    % figure;imagescP(Rs,theta/pi*180,phist)
    % ylabel('\theta [degrees]')
    % xlabel('R')

    str_title=[];
    for ind1=1:length(coeffs)
        str_title=[str_title '(' num2str(roundP(coeffs(ind1),3)) '\pm' num2str(roundP(coeffs_sigma(ind1),3)) ')*' basisnames{ind1} ' + '];
    end
    str_title=str_title(1:end-3);
    
    pax=polarP(hax(2),theta,rho);
    set(pax,'linewidth',2,'color','r')
    [x0,y0]=pol2cart(theta,rho);
    hold on;
    for ind1=round(Nr/20):Nr
        [x(:,ind1),y(:,ind1)]=pol2cart(theta-(theta(2)+theta(1))/2,Rs(ind1)*ones([Ntheta 1]));
        color_line3(y(:,ind1),x(:,ind1),zeros([Ntheta 1]),phist(:,ind1),'parent',hax(2),'linewidth',4);
    end
    line(y0,x0,'parent',hax(2),'linewidth',2,'color','red');
    line([0 0],[-0.8*ax_limit 0.8*ax_limit],[0 0],'parent',hax(2),'color','k','linestyle','--','linewidth',2);
    colormap(hax(2),'gray');
    colormap(hax(2),flipud(colormap(hax(2))));
%     colorbar;
    axis square;
    view([0 90])
    xlim([-ax_limit ax_limit])
    ylim([-ax_limit ax_limit])
    zlim([0 1])
    title(hax(2),str_title)
end
end
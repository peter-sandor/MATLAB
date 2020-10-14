function s=POP(im, bParams, qParams, cartImage)
% Apply polar onion peeling method for analyzing VMI images. The code
% implements the method shown at G.M. Roberts et-al Rev. Sci. Instr. 80, 053104 (2009)
% The code is limited to images no larger than 1024x1024 pixels, and uses
% arbitrary (even) beta parameters.
%
% Inputs:
%  im       - a Velocity map image, the image is assumed to be square
%             matrix, centered and corrected for elliptisity, tilt, etc.
%  bParams  - a vector which each beta parameter to fit for (excluding 0
%             which is always included), e.g. [2 4] for beta2 & beta4
%  qParams  - a vector which states the quadrants of the image to include
%             in the analysis, e.g. [2 3] will treat only the left image
%             side. This is sometimes needed if there is a burned spot, or
%             images are only symmetric in half of the plane, as in two
%             color (w,2w) experiments.
%  cartImage- the kind of Cartesian image to return, acceptable value are
%             'sim' for the simulated (fit) image
%             'exp' for the experimental image
%              0   for no Cartesian image (this speeds up run time)
%
% Outputs:
%  s            - A Matlab structure that contains the following:
%  s.iraraw     - a 2d triangular array (polar) of the raw data folded
%                 to (0,pi/2) - G(R,alpha)
%  s.iradecon   - a 2d triangular array (polar) of the simulated
%                 deconvolved data folded to (0,pi/2) - g_fit(r;R,alpha)
%  s.iradeconExp- a 2d triangular array (polar) of the experimental
%                 deconvolved data folded to (0,pi/2) - g_fit(r;R,alpha)
%  s.sqp        - a 2d rectangular array (polar) derived from iraraw
%  s.PESId      - a 1d vector of the radial projection of the (simulated)
%                 deconvolved data
%  s.PESIdExp   - a 1d vector of the radial projection of
%                 the (experimental) deconvolved data
%  s.Betas      - a matrix containing each beta parameters (as specified in
%                 bParams, e.g. if bParams=[2 4] then s.Betas(1,:) gives b2
%                 and s.Betas(2,:) gives b4
%
% And depending on the option selected for cartImage:
%  s.simImage    - a 2d reconstruction after the onion peeling in Cartesian
%                  coordinates of the simulated image
%  s.expImage    - a 2d reconstruction after the onion peeling in Cartesian
%                  coordinates of the experimental image
%  s.sqpdecon    - a 2d rectangular array (polar) derived from iradecon
%  s.sqpdeconExp - a 2d rectangular array (polar) derived from iradeconExp
%
% Example:
%
% load('testimg.mat');
% s = POP(im,[2 4], 1:4, 'sim');
% figure;
% subplot(2,2,1);imagesc(im);        axis square;title('raw image')
% subplot(2,2,2);imagesc(s.iraraw);  axis square;title('ira raw')
% subplot(2,2,3);imagesc(s.simImage);axis square;title('simulated image')
% subplot(2,2,4);imagesc(s.iradecon);axis square;title('ira deconvoluted')
%
%   Comments \ improvements are welcomed
%   Adi Natan (natan@stanford.edu)
%   Adam Chatterley (a.s.chatterley@dur.ac.uk)

%% defaults
if (nargin < 4);                 cartImage = 0; end ;
if (nargin < 3);  qParams=1:4 ;  cartImage = 0; end ;
if (nargin < 2);  bParams=[2 4];  qParams=1:4 ; cartImage = 0; end ;

% Check that beta Params are in the scope of the code
for i=bParams
    if (mod(i, 2) ~= 0 || i <= 0 || i > 24)
        error('Only even positive beta parameters <=24 supported!');
    end
end

% Trapping specific warning ID in case of nearly singular matrix (line 161)
temperr = warning('error', 'MATLAB:nearlySingularMatrix');

%% inital steps
x0=ceil(size(im,1)/2);y0=ceil(size(im,2)/2); % Assuming image is centered
% load lut - impulse response basis set  (delta functions)
load('delta_lut.mat');
RR=(0:size(im)/2);
PPR=single(floor(0.5*pi*(RR+1))-1); % calc the  # of pixels per radius
AngleInc = single(0.5*pi./PPR'); % angle increment per radius
AngleInc(1)=0; % avoid inf at origin

%set radius of square matrix around center
L = min([x0,y0,(size(im) - [x0,y0])]);

%create  image quadrants
Q(:,:,1) =     im(x0:-1:x0-L+1,y0:y0+L-1);
Q(:,:,2) =     im(x0:-1:x0-L+1,y0:-1:y0-L+1);
Q(:,:,3) =     im(x0:x0+L-1,y0:-1:y0-L+1);
Q(:,:,4) =     im(x0:x0+L-1,y0:y0+L-1);

%add image quadrants
a4=zeros(size(Q,1),size(Q,2));
for ii=qParams;
    a4=a4+Q(:,:,ii);
end

% normalize by total image intensity
a4=a4*4/numel(qParams);

ira = zeros(L-2,PPR(L)); % initialize the  matrix
ira(1,1) = a4(1,1);      % origin pixel remains the same
PESR = zeros(L-2,PPR(L));% initialize the reconstructed matrix
PESR(1,1) = a4(1,1);

%% creating the 2d triangular array polar image
for r=2:L-2
    npr=PPR(r); % determine # polar pix in radius
    angincr=AngleInc(r);    % the angular increment per radius
    qp=0:npr;
    xp=r*sin(angincr*qp)+1;  % polar x-coordinate
    yp=r*cos(angincr*qp)+1;  % polar y-coordinate
    % define scale fractional weight of cart pixels in polar pixels
    xc=round(xp);yc=round(yp);
    xd=1-abs(xc-xp);
    yd=1-abs(yc-yp);
    % gather intensity per radius per angle (ira)
    ira(r,1:npr+1) = xd.*yd.*a4(xc+(yc-1)*L);
    ira(r,2:npr) = ira(r,2:npr) +xd(2:npr).*(1-yd(2:npr)).*a4(xc(2:npr)+...
        (yc(2:npr)-1+(-1).^(yp(2:npr)<yc(2:npr)))*L) + ...
        (1-xd(2:npr)).*yd(2:npr).*a4(xc(2:npr)+...
        (-1).^(xp(2:npr)<xc(2:npr))+yc(2:npr)*L) + ...
        (1-xd(2:npr)).*(1-yd(2:npr)).*a4(xc(2:npr)+...
        (-1).^(xp(2:npr)<xc(2:npr))+L*(yc(2:npr)+...
        (-1).^(yp(2:npr)<yc(2:npr))));
    
    PESR(r,1:npr+1) = r - 1 + ira(r,1:npr+1);
end

iraraw = ira; % save to compare later onion peeling

%initialize some more  parameters
betas = zeros(length(bParams),L-2);
PESId = zeros(1,L-2);
iradecon = zeros(L-2,PPR(L));
iradeconExp = zeros(L-2,PPR(L));
PESIdExp = zeros(1,L-2);

%% Onion peel:
for r = L-2:-1:2
    B = zeros(1,numel(bParams+1));
    npr=PPR(r); % # of polar pixels in radius -1
    qp=0:npr;
    y = ira(r,qp+1); % assign row of all pixels in radius r
    B(1)=sum(y);
    
    % one fit coefficient for each B param
    fitCoefs = ones(numel(bParams + 1), PPR(r)+1);
    for ii=(1:numel(bParams))
        % assign relevant Legendre polynomials to fitCoef
        fitCoefs(ii+1,:) = leg(bParams(ii), cos(AngleInc(r)*qp));
        B(ii+1) = y*fitCoefs(ii+1,:)'; % fill B matrix
    end
    
    A = fitCoefs * fitCoefs'; % matrix A for least square fitting
    A(1,1)=npr+1;
    
    try
        Ain = inv(A);
    catch % switch to Moore-Penrose pseudoinverse in case of
        % nearly singular matrix
        Ain =  pinv(A);
    end
    
    Beta = zeros(1,length(bParams+1));
    Beta(1) = B*Ain(:,1);
    
    for ii=1:numel(bParams)
        Beta(ii+1)  = B*Ain(:,ii+1)/Beta(1); % generate beta matrix
        betas(ii,r) = Beta(ii+1); % copy for output betas
    end
    
    if 0==Beta(1)
        PESId(r)=0;
        continue;
    end
    
    % generate matrices for alpha; R/rp * cos(alpha); and the basis set
    % scaled by pixels per radius
    alphaMat = (AngleInc(1:r) * (0:PPR(r)));
    rrpCosAlphaMat = repmat((1:r)'/r, 1, PPR(r)+1).*cos(alphaMat);
    
    itMat = repmat((lut(r,1:r)*(npr+1)./ (PPR(1:r)+1))', 1, PPR(r)+1);
    
    bContrib = ones(r, PPR(r)+1);
    for ii=1:numel(bParams) %  add each beta contribution
        bContrib = bContrib + Beta(ii+1)*leg(bParams(ii), rrpCosAlphaMat);
    end
    
    % generate the simulated image for this radius
    factMat = Beta(1).*itMat.*bContrib;
    % save the simulated data
    PESId(r) = PESId(r) + sqrt(r)*factMat(end,:)*sin(alphaMat(end,:))';
    iradecon(r,1:PPR(r)+1)=factMat(end,1:PPR(r)+1)/sqrt(r);
    % save the experimental data
    PESIdExp(r) = PESId(r) + sqrt(r)*ira(r,1:PPR(r)+1)*sin(alphaMat(end,:))';
    iradeconExp(r,1:PPR(r)+1)=ira(r,1:PPR(r)+1)/sqrt(r);
    % subtract away the simulated image
    ira(1:r,1:PPR(r)+1) = ira(1:r,1:PPR(r)+1)-factMat;
end

%% assign zero value to infs and NaNs if present.
PESId(~isfinite(PESId))=0;
iradecon(~isfinite(iradecon))=0;
iradeconExp(~isfinite(iradeconExp))=0;

%%  Create a square polar projection from the triangular one
[M,N] = size(iraraw);
% # of non-zero polar pixels per radius in triangular polar plot
nzpr = sum(iraraw>0,2);
% create a matrix that "squares" the triangular polar matrix
sqp=zeros(M,N);
sqp(1,:)=ones(1,N).*iraraw(1,1);

for ii=2:M
    sqp(ii,:)=interp1(1:nzpr(ii),iraraw(ii,1:nzpr(ii)),linspace(1,nzpr(ii),N),'spline');
end

%% 2d transform to Cartesian coordinates of the simulated\exp image
if strcmp(cartImage, 'sim')
    x = [];
    y = [];
    z = [];
    for r = 1:L-5
        qp = 0:PPR(r);
        x = [x r*cos(qp*AngleInc(r))];
        y = [y r*sin(qp*AngleInc(r))];
        z = [z iradecon(r,qp+1)];
    end
    
    % Use proper interpolator according to the Matlab version used
    if verLessThan('matlab', '8.1')
        F = TriScatteredInterp(double(x(:)),double(y(:)),double(z(:)),'natural');
    else
        F = scatteredInterpolant(double(x(:)),double(y(:)),double(z(:)),'natural');
    end
    
    [xx,yy] = meshgrid(0:L-2,0:L-2);
    zz = F(xx,yy);
    zz = max(zz,zeros(size(zz)));
    zz = [zz(end:-1:1,end:-1:1) zz(end:-1:1,:); zz(:,end:-1:1) zz];
    
    s=struct('iraraw',iraraw,'iradecon',iradecon, 'iradeconExp', iradeconExp, 'sqp',sqp,'PESId',PESId,'PESIdExp',PESIdExp,'Betas',betas,'simImage',zz);
    return
    
elseif (strcmp(cartImage,'exp'))
    x = [];
    y = [];
    z = [];
    for r = 1:L-5
        qp = 0:PPR(r);
        x = [x r*cos(qp*AngleInc(r))];
        y = [y r*sin(qp*AngleInc(r))];
        z = [z iradeconExp(r,qp+1)];
    end
    
    % Use proper interpolator according to the Matlab version used
    if verLessThan('matlab', '8.1')
        F = TriScatteredInterp(double(x(:)),double(y(:)),double(z(:)),'natural');
    else
        F = scatteredInterpolant(double(x(:)),double(y(:)),double(z(:)),'natural');
    end
    
    [xx,yy] = meshgrid(0:L-2,0:L-2);
    zz = F(xx,yy);
    zz = max(zz,zeros(size(zz)));
    zz = [zz(end:-1:1,end:-1:1) zz(end:-1:1,:); zz(:,end:-1:1) zz];
    
    s=struct('iraraw',iraraw,'iradecon',iradecon, 'iradeconExp', iradeconExp,'sqp',sqp, 'PESId',PESId,'PESIdExp',PESIdExp,'Betas',betas,'expImage',zz);
    
    return
end
s=struct('iraraw',iraraw,'iradecon',iradecon,'iradeconExp', iradeconExp, 'sqp',sqp,'PESId',PESId,'PESIdExp',PESIdExp, 'Betas',betas);


function p=leg(m,x)
%  This function returns Legendre polynomial P_m(x) where m is the degree
%  of polynomial and X is the variable.
%  The x2=x.*x is a optimized to minimize the # of operations.
switch m
    case 0
        p=ones(size(x));
        return
    case 1
        p=x;
        return
    case 2
        p=(3*x.*x -1)/2;
        return
    case 4
        x2=x.*x;
        p = ((35.*x2-30).*x2+3)/8;
        return
    case 6
        x2=x.*x;
        p = (((231.*x2-315).*x2+105).*x2-5)/16;
        return
    case 8
        x2=x.*x;
        p = ((((6435.*x2-12012).*x2+6930).*x2-1260).*x2+35)/128;
        return
    case 10
        x2=x.*x;
        p = (((((46189.*x2-109395).*x2+90090).*x2-30030).*x2+3465).*x2-63)/256;
        return
    case 12
        x2=x.*x;
        p = ((((((676039.*x2-1939938).*x2+2078505).*x2-1021020).*x2+225225).*x2-18018).*x2+231)/1024;
        return
    case 14
        x2=x.*x;
        p = (((((((5014575.*x2-16900975).*x2+22309287).*x2-14549535).*x2+4849845).*x2-765765).*x2+45045).*x2-429)/2048;
        return
    case 16
        x2=x.*x;
        p = ((((((((300540195.*x2-1163381400).*x2+1825305300).*x2-1487285800).*x2+669278610).*x2-162954792).*x2+19399380).*x2-875160).*x2+6435)/32768;
        return
    case 18
        x2=x.*x;
        p = (((((((((2268783825.*x2-9917826435).*x2+18032411700).*x2-17644617900).*x2+10039179150).*x2-3346393050).*x2+624660036).*x2-58198140).*x2+2078505).*x2-12155)/65536;
        return
    case 20
        x2=x.*x;
        p = ((((((((((34461632205.*x2-167890003050).*x2+347123925225).*x2-396713057400).*x2+273491577450).*x2-116454478140).*x2+30117537450).*x2-4461857400).*x2+334639305).*x2-9699690).*x2+46189)/262144;
        return
    case 22
        x2=x.*x;
        p = (((((((((((263012370465.*x2-1412926920405).*x2+3273855059475).*x2-4281195077775).*x2+3471239252250).*x2-1805044411170).*x2+601681470390).*x2-124772655150).*x2+15058768725).*x2-929553625).*x2+22309287).*x2-88179)/524288;
        return
    case 24
        x2=x.*x;
        p = ((((((((((((8061900920775.*x2-47342226683700).*x2+121511715154830).*x2-178970743251300).*x2+166966608033225).*x2-102748681866600).*x2+42117702927300).*x2-11345993441640).*x2+1933976154825).*x2-194090796900).*x2+10039179150).*x2-202811700).*x2+676039)/4194304;
        return
end


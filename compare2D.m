function varargout=compare2D(frq1a,frq3a,data2Da,frq1b,frq3b,data2Db,varargin)

% frqvec=1.135:0.001:1.170;
if nargin==6
    frq_min=max([min(frq1a) min(frq1b) min(frq3a) min(frq3b)]);
    frq_max=min([max(frq1a) max(frq1b) max(frq3a) max(frq3b)]);
    frq_step=frq1a(2)-frq1a(1);
    frqvec1=frq_min:frq_step:frq_max;
    frqvec3=frqvec1;
elseif nargin==8
    frqvec1=varargin{1};
    frqvec3=varargin{2};
end
spa=real(data2Da);
spb=real(data2Db);
[FRQ1a FRQ3a]=ndgrid(frq1a,frq3a);
[FRQ1b FRQ3b]=ndgrid(frq1b,frq3b);
[FRQ1ip FRQ3ip]=ndgrid(frqvec1,frqvec3);
spipa=interp2(FRQ3a,FRQ1a,spa,FRQ3ip,FRQ1ip);
spipb=interp2(FRQ3b,FRQ1b,spb,FRQ3ip,FRQ1ip);

spsuma=sum(sum(spipa));
spsumb=sum(sum(spipb));
% spsuma=max(max(spipa));
% spsumb=max(max(spipb));
varargout{1}=sign(spipa/spsuma-spipb/spsumb).*(spipa/spsuma-spipb/spsumb).^2;
varargout{2}=frqvec1;
varargout{3}=spipa/spsuma;
varargout{4}=spipb/spsumb;
end
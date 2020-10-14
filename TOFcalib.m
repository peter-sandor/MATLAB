function varargout = TOFcalib(t1,m1,t2,m2,varargin)
% t=A*sqrt(m/q)+t0
% m=((t-t0)/A)^2*q
% Arguments:
%   delay       - arrival time axis in [ns] or [pixel]
%   q           - charge of fragments in AU
%   t1          - arrival time of fragment1 in [ns]
%   m1          - mass of fragment1 in [AMU]
%   t2          - arrival time of fragment2 in [ns]
%   m2          - mass of fragment2 in [AMU]
% Returned variable:
%   mass        - calibrated mass axis in [AMU]

if nargin==4
    delay=1:8000;
    q=1;
elseif nargin==5
    delay=varargin{1};
    q=1;
elseif nargin==6
    delay=varargin{1};
    q=varargin{2};
end
A=sqrt(q)*(t2-t1)/(sqrt(m2)-sqrt(m1));
t0=t1-A*sqrt(m1/q);
mass=((delay-t0)/A).^2;
varargout{1}=mass;
varargout{2}=t0;
varargout{3}=A;
end
function varargout = Morse(varargin)
r=varargin{1};
re=varargin{2};
a=varargin{3};
De=varargin{4};
n=varargin{5};

hbar=1;
mass=1;
x=a*r;
xe=a*re;
lambda=sqrt(2*mass*De)/a/hbar;
z=2*lambda*exp(-x+xe);
N=factorial(n)*sqrt(a*(2*lambda-2*n-1)/gamma(n+1)/gamma(2*lambda-n));

potential=De*(1-exp(-a*(r-re))).^2-De;
if nargout>=2;
    nrg=[];
    for n=0:floor(lambda)
        nrg=[nrg; -a^2*hbar^2/2/mass*(lambda-n-1/2)^2];
    end
    varargout{2}=nrg;
    if nargout==3
        PSI=[];
        for n=0:floor(lambda)
            PSI=cat(1,PSI,N.*z.^(lambda-n-1/2).*exp(-z/2).*polyval(LaguerreGen(n,2*lambda-2*n-1),z));
        end
        varargout{3}=PSI;
    end
end
varargout{1}=potential;
end
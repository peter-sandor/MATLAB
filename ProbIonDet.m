function P = ProbIonDet(lambda,Ngen,Nelectron,Nion,p,q,varargin)
% This function calculates the probability of generating 'Ngen'
% ion+electron pairs in the laser field and detecting exactly 'Nelectron'
% electrons and 'Nion' ions of those, given the electron detection
% efficiency 'p' and ion detection efficiency 'q'.
% lambda is the expectation value of the number of pairs generated,
% assuming Poissonian statistics.
% varargin adds two optional inputs if the different ion fragments need to
% be treated separately: the first of the two are the occurences (xi), the
% second of the two (pi) are the probabilites of generating the individual
% fragments (both arguments are vectors, with lengths equal to the number of fragment types.)
if nargin==6
    P = nchoosek(Ngen,Nelectron).*p.^Nelectron*(1-p).^(Ngen-Nelectron).*nchoosek(Ngen,Nion).*q^Nion.*(1-q).^(Ngen-Nion).*lambda.^Ngen./factorial(Ngen).*exp(-lambda);
elseif nargin==8
    % xi = varargin{1};
    % pi = varargin{2};
	probfactor=1;
    for ind1=1:length(varargin{1})
        probfactor=probfactor*varargin{2}(ind1)^varargin{1}(ind1);
    end
    P = nchoosek(Ngen,Nelectron).*p.^Nelectron*(1-p).^(Ngen-Nelectron).*nchoosek(Ngen,Nion).*q^Nion.*(1-q).^(Ngen-Nion).*nchooseks(Nion,varargin{1}).*probfactor.*lambda.^Ngen./factorial(Ngen).*exp(-lambda);
end
end
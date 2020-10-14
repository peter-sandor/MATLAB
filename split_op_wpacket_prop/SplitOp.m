function varargout = SplitOp(input)
% This is the core solver for the 1D wavepacket propagation code
% Based on David Tannor's book, section 11.7.1
% input: struct with the following fields:
%       Psi0: Initial wavefunction as a function of position 'x', 1D array
%		V: potential energy as a function of position 'x' and possibly of time 't'; either 1D (Nx-by-1) or 2D (Nx-by-Nt) array
%		x: 1D coordinate array for the position
%		p: 1D coordinate array for the momentum
%		t: 1D array containing the times at which the wavefunction should be evaluated
% Psi0, V, x and p should have the same number of elements!
%
% The output is a struct with the following fields:
% 		x,p,V,t: same as for the input
% 		F: struct which can be used to play the result of the propagation as a movie. See the MATLAB command 'movie'.
%		autocorr: 1D array, the autocorrelation of the wavefunction at time 't' with that at the start time
%		spectrum: 1D array, the Fourier-transform of the autocorrelation
%		exp_value: 1D array, <Psi*|x|Psi>
%		std_dev: 1D array, <Psi*|x^2|Psi>-(<Psi*|x|Psi>)^2
%		Psi: 2D array, Psi(x,t) (position-space wavefunction)
%		Psi_F: 2D array, Psi(p,t) (momentum-space wavefunction)
	
tic;
Nt=length(input.t);
Nx=length(input.x);
if size(input.V,2)==1
    input.V=extend(input.V,Nt);
end
t_step=mean(diff(input.t));
Psi=zeros([length(input.x) length(input.t)]);
Psi_F=zeros([length(input.x) length(input.t)]);
Psi(:,1)=input.Psi0;
exp_value=zeros(size(input.t));
std_dev=zeros(size(input.t));
autocorr=zeros(size(input.t));
fprintf(1, 'Calculating:\n');
for ind1=1:Nt
    Psi(:,ind1)=squeeze(Psi(:,ind1)).*map2colvec(exp(-i*input.V(:,ind1)*t_step/2));
    Psi_F(:,ind1)=fft(squeeze(Psi(:,ind1)));
    Psi_F(:,ind1)=map2colvec(Psi_F(:,ind1)).*map2colvec(exp(-i*1/2*input.p.^2*t_step));
    Psi(:,ind1)=ifft(Psi_F(:,ind1));
    Psi(:,ind1)=squeeze(Psi(:,ind1)).*map2colvec(exp(-i*input.V(:,ind1)*t_step/2));
    if ind1<Nt
        Psi(:,ind1+1)=Psi(:,ind1);
    end
    exp_value(ind1)=sum(conj(Psi(:,ind1)).*input.x.*Psi(:,ind1));
    std_dev(ind1)=sqrt(sum(conj(Psi(:,ind1)).*input.x.^2.*Psi(:,ind1))-exp_value(ind1)^2);
    autocorr(ind1)=map2rowvec(conj(Psi(:,ind1)))*map2colvec(input.Psi0);
    fprintf(1, '.');
    if mod(ind1,100)==0
        fprintf(1, '\n');
    end
end
fprintf(1, '\n');
disp('Done.');
spectrum=fft(autocorr);

out.x=input.x;
out.p=input.p;
out.V=input.V;
out.t=input.t;
out.frq=FourierAxis(input.t)/2/pi;
out.autocorr=autocorr;
out.spectrum=spectrum;
out.exp_value=exp_value;
out.std_dev=std_dev;
out.Psi=Psi;
out.Psi_F=Psi_F;
varargout{1}=out;
toc
end
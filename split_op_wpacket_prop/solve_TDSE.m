% solve TDSE with split-operator method (see David Tannor's book at section 11.7.1)
% use Hartree atomic units
% written by Peter in 2013 fall
% All calculations are 1D, on a single potential energy surface
% set integer variable 'id' to select which problem you want to solve:
% id=1: free-particle wavepacket
% id=2: wavepacket in harmonic oscillator potential
% id=3: wavepacket in Morse potential
% id=4: wavepacket in predissociative potential
% id=5: wavepacket scattering from ramp-up potential
% id=6: free electron in pulsed oscillating electric field
% id=7: electron in Morse potential + pulsed oscillating electric field
% id=8: Coulomb potential with absorbing boundary conditions
% id=9: Coulomb potential with absorbing boundary conditions + pulsed
%   oscillating electric field

id=9;
switch id
    case 1
        Nt=2^8;
        Nx=2^10;
        x_min=0;
        x_max=512;
        x_step=(x_max-x_min)/(Nx-1);
        input.x=map2colvec(x_min:x_step:x_max);
        input.potential=zeros([Nx 1])-1e-1*1i*[zeros([Nx-40 1]); map2colvec(0:39)]; % free particle with absorbing boundary condition on the right
        t_min=0;
        t_max=40;
        ind0=200;
        m0=40;
    case 2
        Nt=2^8;
        Nx=2^10;
        x_min=0;
        x_max=512;
        x_step=(x_max-x_min)/(Nx-1);
        input.x=map2colvec(x_min:x_step:x_max);
        input.potential=5e-6*(input.x-256).^2; % harmonic
        t_min=0;
        t_max=1e4;
        ind0=500;
        m0=1;
    case 3
        Nt=2^8;
        Nx=2^10;
        x_min=0;
        x_max=512;
        x_step=(x_max-x_min)/(Nx-1);
        input.x=map2colvec(x_min:x_step:x_max);
        [a,b]=Morse(map2rowvec(input.x),120,0.009,1e0,1); % Morse
        input.potential=a;
        t_min=0;
        t_max=5e3;
        ind0=160;
        m0=1;
    case 4
        Nt=2^10;
        Nx=2^10;
        x_min=-100;
        x_max=20;
        x_step=(x_max-x_min)/(Nx-1);
        input.x=map2colvec(x_min:x_step:x_max);
        input.potential=1/10*(1/2*input.x.^2+0.03*(input.x+10).^3).*1./(1+exp(-1/2*(input.x+40)))-10./(1+exp(1/2*(input.x+40)))-5e-1*1i*[map2colvec(39:-1:0); zeros([Nx-40 1])]; % predissociative with absorbing boundary condition on the left
        t_min=0;
        t_max=80;
        ind0=914;
        m0=1;
    case 5
        Nt=2^8;
        Nx=2^10;
        x_min=0;
        x_max=512;
        x_step=(x_max-x_min)/(Nx-1);
        input.x=map2colvec(x_min:x_step:x_max);
        input.potential=5e-2*[zeros([512 1]); map2colvec(0:199); zeros([Nx-512-200 1])]; % ramp up potential barrier
        t_min=0;
        t_max=100;
        ind0=400;
        m0=400;
    case 6
        Nt=2^8;
        Nx=2^10;
        x_min=-Nx;
        x_max=+Nx;
        x_step=(x_max-x_min)/(Nx-1);
        input.x=map2colvec(x_min:x_step:x_max);
        t_min=0;
        t_max=2e3;
        t_step=(t_max-t_min)/(Nt-1);
        input.t=t_min:t_step:t_max;
        t0=1e3;
        phi=pi;
        omega=2*pi/(5*2.7/24.2e-3);
%         omega=2*pi/(10/24.2e-3);
        tau=15/24.2e-3; % temporal extent of field [24 as]
        field_amp=3*0.0053;
        envelope=field_amp/omega*exp(-(input.t-t0).^2/tau^2);
        A=envelope.*cos(omega*input.t+phi);
        field=diff(A)/t_step;
        field=[0 field];
        input.potential=map2colvec(input.x)*map2rowvec(field); % free electron in pulsed oscillating electric field
        ind0=512;
        m0=1;   
    case 7
        Nt=2^9;
        Nx=2^10;
        x_min=0;
        x_max=+2*Nx;
        x_step=(x_max-x_min)/(Nx-1);
        input.x=map2colvec(x_min:x_step:x_max);
        t_min=0;
        t_max=5e3;
        t_step=(t_max-t_min)/(Nt-1);
        input.t=t_min:t_step:t_max;
        t0=2.5e3;
        [a,b]=Morse(map2rowvec(input.x),120,0.009,1e0,1);
        phi=pi;
%         omega=2*pi/(2.7/24.2e-3);
        omega=2*pi/(27/24.2e-3);
        tau=15/24.2e-3; % temporal extent of field [24 as]
        field_amp=1*0.0053;
        envelope=field_amp/omega*exp(-(input.t-t0).^2/tau^2);
        A=envelope.*cos(omega*input.t+phi);
        field=diff(A)/t_step;
        field=[0 field];
        input.potential=map2colvec(input.x)*map2rowvec(field)+extend(map2colvec(a),Nt); % electron in Morse potential + pulsed oscillating electric field
        ind0=160;
        m0=1;
    case 8
        Nt=2^8;
        Nx=2^10;
        x_min=-Nx/4;
        x_max=+Nx/4;
        x_step=(x_max-x_min)/(Nx-1);
        input.x=map2colvec(x_min:x_step:x_max);
        t_min=0;
        t_max=0.5e3;
        t_step=(t_max-t_min)/(Nt-1);
        input.t=t_min:t_step:t_max;
        t0=2.5e3;
        % 1D Coulomb potential with absorbing boundary conditions at +Inf
        % and -Inf:
        input.potential=-abs(1./input.x)-1e-1*1i*[zeros([Nx-40 1]); map2colvec(0:39)]-1e-1*1i*[map2colvec(0:39); zeros([Nx-40 1])];
%         input.potential(input.potential<1000)=-1000;
        ind0=Nx/2;
        m0=1;
    case 9
        Nt=2^9;
        Nx=2^10;
        x_min=-Nx/4;
        x_max=+Nx/4;
        x_step=(x_max-x_min)/(Nx-1);
        input.x=map2colvec(x_min:x_step:x_max);
        t_min=0;
        t_max=5e3;
        t_step=(t_max-t_min)/(Nt-1);
        input.t=t_min:t_step:t_max;
        t0=2.5e3;
        phi=pi;
%         omega=2*pi/(2.7/24.2e-3);
        omega=2*pi/(27/24.2e-3);
        tau=15/24.2e-3; % temporal extent of field [24 as]
        field_amp=1*0.0053;
        envelope=field_amp/omega*exp(-(input.t-t0).^2/tau^2);
        A=envelope.*cos(omega*input.t+phi);
        field=diff(A)/t_step;
        field=[0 field];
        % 1D Coulomb potential with absorbing boundary conditions at +Inf
        % and -Inf:
        input.potential = map2colvec(input.x)*map2rowvec(field)-1*abs(1./input.x)-1e-1*1i*[zeros([Nx-40 1]); map2colvec(0:39)]-1e-1*1i*[map2colvec(0:39); zeros([Nx-40 1])];
%         input.potential(input.potential<1000)=-1000;
        ind0=Nx/2;
        m0=1;
end
t_step=(t_max-t_min)/(Nt-1);
input.t=t_min:t_step:t_max;
input.p=map2colvec([0:(Nx/2-1) -(Nx/2):-1]/Nx/x_step*2*pi);
% frq=[0:Nt/2-1 -(Nt/2):-1]/Nt/t_step;

x0=input.x(ind0);
sigma=8;
% sigma=0.0016;
sigma_mom=1/4/sigma;
input.Psi0=zeros(size(input.x));

% construct initial wavepacket in momentum space
for ind1=1:Nx
    input.Psi0=input.Psi0+exp(-(input.p(ind1)-input.p(m0)).^2/sigma_mom^2).*exp(1i*input.p(ind1)*(map2colvec(input.x-input.x(ind0))));
end
input.Psi0=input.Psi0./sqrt(sum(abs(input.Psi0).^2));
input.V=map2colvec(input.potential);
% this is where magic happens!
solution = SplitOp(input); % propagator code
%% retrieve an eigenfunction corresponding to a specific energy 'E'
omega=2*pi*(-0.0168);
eigenfnct=sum(permute(extend(exp(1i*omega*solution.t),Nx),[2 1]).*solution.Psi,2);
dfigure;
plot(solution.x,[abs(eigenfnct) real(eigenfnct)])
hold on;plot(solution.x,real(solution.V-solution.V(ind0)),'k-')
hold off;
title(['E= ' num2str(omega/27.2) ' eV'])
%% create movie
F=TDSE_movie(solution);
fig2=figure;
movie(fig2,F,1,15)
% movie2avi(aviobj,'test.avi')
%%
figure;imagesc(solution.t,solution.x,abs(solution.Psi).^2)
xlabel('delay [at.u.]')
ylabel('distance [at.u.]');
title('|\Psi(x,t)|^2')
%%
figure(2)
subplot(222)
plot(solution.t,abs(solution.autocorr),'k-')
xlabel('delay [at.u.]')
ylabel('C(t)')
subplot(221)
plot(solution.t,real(solution.exp_value),'k-')
xlabel('delay [at.u.]')
ylabel('expectation value')
subplot(223)
plot(solution.t,abs(solution.std_dev),'k-')
xlabel('delay [at.u.]')
ylabel('standard deviation')
subplot(224)
plot(fftshift(solution.frq),fftshift(abs(solution.spectrum)),'k-')
xlabel('\nu [at.u.]')
ylabel('power spectrum')
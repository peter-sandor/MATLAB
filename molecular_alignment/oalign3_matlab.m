function [output] = oalign3_matlab(input)
%     CALCULATES TIME DEPENDENT EVOLUTION OF <COS(THETA)>,
%     <COS(THETA).^2) AND THE TIME
%     DEPENDENT ANGULAR DISTRIBUTION (IN THETA) FOR A LINEAR RIGID
%	  ROTOR MOLECULE THAT HAS BEEN SUBJECTED TO A Single Cycle HCP	
%	  AND AN INTENSE LASER PULSE
%	  INTEGRATION STEPSIZE IS VARIED TO KEEP THE STATE NORMALIZED TO
%	  UNITY WITHIN DESIRED ACCURACY
%	
%	  INITIAL BOLTZMANN STATE DISTRIBUTION ASSUMED
%	  PEAK HCP FIELD AND DURATION, LASER INTENSITY AND DURATION,
%	  LASER PULSE DELAY, MOLECULAR DIPOLE MOMENT, ROTATIONAL
%	  CONSTANT AND ROTATIONAL TEMPERATURE ARE INPUTS
%     
%       REAL*8 FHCP, TAU, R0, B, E(103), V(103), TSTEP, RP(203,2503)
%       REAL*8 W(103), TCOS(2503),TMAX, RC(103),IC(103),THETASTEP,W2(103)
%       REAL*8 F, T, RCI(103), ICI(103), RT(103), IT(103),NORM,P(103,203)
%       REAL*8 ACCUR, OLDNORM, TT, RRC(103), IIC(103), AMP(103), TK, AD
%       REAL*8 BIGAMP,JNORM,PIECE,TOTNORM,SUMAMP,THETA(203),PROB(203,2503)
%       REAL*8 SOMX2, PLGNDR, PMMP1, PLL, FACT, PMM, FAC, X, IP(203,2503)
%       REAL*8 TCOS2(2503), I0, ILASER, TAUL, VJJ(103), VJ2(103),DELAY
%       REAL*8 PIECE2, D, TTEST1, OLDTTEST1, TTEST2, OLDTTEST2, EVEN, ODD
%       REAL*8 TSTART, TSTART1, TSTART2
%       INTEGER*4 JMAX, M, MMAX, JJ, J, MAXJ, JJMAX, MAXT, BIGJ, BIGT
%       INTEGER*4 MMAXJ, AMPJ, THETANUM, I, BIGTHETA
%       INTEGER*4 L, LL, II,III, K
%       CHARACTER*12 FIL, FIL2

%     CHARACTER*8 FILO2, JUNK
%     CHARACTER*1 SUF3,SUF2,SUF1
%     CHARACTER*4 SUF
%       PI=4.0*atan(1.0)	
BIGJ =100;
BIGT= 2500;
BIGTHETA = 200;
TOTNORM = 0;
MMAXJ = 0;
%       OPEN(1,STATUS='OLD',FILE='oalign-3.dat')
%    'ENTER THE COS(THETA) AND PSI(THETA) FILENAME-
%             MUST BE 8 CHAR. LONG-SUFFIX WILL BE '.<ANGLE>'
%      READ(1,7)FILO2
%       READ(1,6)FIL
% 	  READ(1,6)FIL2
% 7     FORMAT(A8)
% 6     FORMAT(A12)    
FHCP=-input.FHCP/5.142E6; % 'ENTER THE PEAK HCP FIELD (1 KV/CM).'
TAU=input.TAU/2.42E-5; % 'ENTER THE HCP DURATION (GAUSSIAN FWHM) IN PSEC.'   
I0=-input.I0/3.5E16/4.0; % 'ENTER THE PEAK LASER INTENSITY IN (W/CM.^2)' and CONVERT I0 TO -E**2/4
TAUL=input.TAUL/2.42E-5; % 'ENTER THE LASER PULSE DURATION (GAUSSIAN FWHM) IN PSEC.'
DELAY = input.DELAY/2.42E-5; % 'ENTER THE RELATIVE DELAY OF THE LASER PULSE IN PSEC.'
AD=input.AD; % 'ENTER THE POLARIZABILITY ANISOTROPY IN (A.U.)'
R0=input.R0; % 'ENTER THE MOLECULAR DIPOLE MOMENT IN (A.U.)'
B=input.B/(2.0*109737.0); % 'ENTER ROTATIONAL CONSTANT IN 1/CM.'
D=input.D/(2.0*109737.0); % 'ENTER CENTRIFUGAL DISTORTION CONSTANT IN 1/CM.'
% 'ENTER THE RELATIVE ABUNDANCE FACTORS FOR EVEN AND ODD J STATES'
EVEN=input.EVEN;
ODD=input.ODD;
JMAX=input.JMAX; % 'ENTER THE MAXIMUM NUMBER OF J STATES TO CONSIDER.  '
if JMAX > BIGJ
    disp('MAX J EXCEEDS CURRENT ARRAY LIMITS')
end    

%    	READ(1,*)TSTEP
TSTEP0=1.0/(20.0*B*JMAX); % 'ENTER THE INITIAL TIME STEP IN PSEC.  '
TSTEP = TSTEP0;
TSTART1=-3.5*TAU+TSTEP/2.0;
TSTART2=DELAY-3.5*TAUL+TSTEP/2.0;
if TSTART2 < TSTART1
  TSTART=TSTART2;
else
  TSTART=TSTART1;
end
%	  B =0.037 (I2); 2.0 (N2); 1.43 (O2); 30.44 (D2); 0.89 (F2)
%	  B=1.34 (LiF); 0.71(LiCl); 7.5 (LiH); 1.93(CO); 1.7(NO); 0.72(SO)
%     B=0.203 (COS)     	
TK=input.TK; % 'ENTER THE ROTATIONAL TEMPERATURE IN K.  '
TMAX=input.TMAX/2.42E-5; % 'ENTER THE MAX DELAY IN PSEC.
DT=input.DT/2.42E-5; % 'ENTER DELAY BETWEEN DATA POINTS IN PSEC.'
NT=round((TMAX-TSTART)/DT);
if (NT > BIGT)
    disp('MAX TIME EXCEEDS CURRENT ARRAY LIMITS.');
end
ACCUR=input.ACCUR; % 'ENTER THE SINGLE STEP NORMALIZATION ACCURACY.  '
THETANUM=input.THETANUM; % 'ENTER THE NUMBER OF THETA STEPS FOR DENSTY PLOT.'  
THETASTEP=pi/(THETANUM-1);

TCOS=zeros([1 NT]);
TCOS2=zeros([1 NT]);
PROB=zeros([THETANUM NT]);
% Pleg=zeros([THETANUM+1 JMAX 2*JMAX+1]);
THETA=(0:THETANUM)*THETASTEP;
TK=TK*1.38D-23;
TK=TK/1.6E-19;
TK=TK/27.2;
TK=TK/B;
%	  TEMPERATURE IN UNITS OF B
%	  CHOOSE MAX M VALUE
%	  PROB = exp(-B*x*(x+1)/(k*T))/(k*T/B)
%	  MMMAX=-LOG(TK/1000.0)*TK
%	  IF (MMMAX .LT. 2) MMAX =2
%	  IF (MMMAX .GT. 2) MMAX=INT(sqrt(MMMAX)+1)
%	  V(1)=0.0
%	  E(1)=0.0
%	  W(1)=0.0
SUMAMP = 0.0;
W(1)=0.0;
E(1)=0.0;
W2(1)=0.0;
W2(2)=0.0;
for J1=2:JMAX+2 % 10
    E(J1)=(J1-1)*(J1)*(B-D*J1*(J1-1));
    W(J1)=E(J1)-E(J1-1);
    if (J1~=2)
        W2(J1)=E(J1)-E(J1-2);
    end
    AMP(J1)=(2.0*(J1-1)+1.0)*exp(-(J1-1)*J1/TK)/TK;
    SUMAMP = SUMAMP+AMP(J1);
end
AMP(1) = 1.0-SUMAMP;
if AMP(1) < 0
    AMP(1) = 0.0;
end
if AMP(1)> 1
    AMP(1) = 1.0;
end
BIGAMP=AMP(1);
AMPJ=1;
for J1=1:JMAX+1 % 13
    AMP(J1) = AMP(J1)/(2.0*(J1-1)+1.0);
    if ((J1-1)/2.0-round((J1-1)/2.0)) < 0.1
        AMP(J1)=2.0*EVEN*AMP(J1)/(EVEN+ODD);
%         WRITE(*,*)'EVEN',J1-1
    else
        AMP(J1)=2.0*ODD*AMP(J1)/(EVEN+ODD);
%         WRITE(*,*)'ODD', J1-1
    end
    if (AMP(J1) > BIGAMP)
        BIGAMP = AMP(J1);
        AMPJ = J1;
    end
%		WRITE(*,*) AMP(J)
%		PAUSE
end
%	  PAUSE
for J=AMPJ:JMAX+1 % 11
    if AMP(J)<(BIGAMP/1000)
        break;
    end
end
MMAX = J; % 12
disp(['MAX THERMALLY POPULATED J STATE IS ' num2str(MMAX)]);
%	  M LOOP
RC(1)=0.0;
IC(1)=0.0;
RT(1)=0.0;
IT(1)=0.0;
ICI(1)=0.0;
RCI(1)=0.0;
RRC(1)=0.0;
IIC(1)=0.0;
M=0;
flag_exit1=0;
flag_exit2=0;

while (M < MMAX) && (flag_exit1==0)

    for J=M+1:(JMAX+2) % loop 20
        % CALCULATE COUPLINGS
        V(J)=sqrt(((J-1).^2-M^2)/((2.0*(J-1)+1.0)*(2.0*(J-1)-1.0)));
        VJJ(J)=(2.0*(J-1)*J-2.0*M^2-1.0)/((2.0*(J-1)+3.0)*(2.0*(J-1)-1.0));
        VJ2(J)=sqrt(((J-1).^2-M^2)*((J-2).^2-M^2)/((2.0*(J-1)+1.0)*(2.0*(J-1)-1.0)*(2.0*(J-1)-1.0)*(2.0*(J-1)-3.0)));

        % CALCULATE LEGENDRE POLYNOMIALS
        L=J-1; % site 5000
        for I2=1:THETANUM+1 % loop 5100
            X=cos(THETA(I2));
            P(J,I2)=0.0;
            if((M<0) || (M>L) || (abs(X)>1))
                break;
            end
            PMM=1.0;
            if(M>0)
                SOMX2=sqrt((1.0-X)*(1.0+X));
                FACT=1.0;
                for I1=1:M % loop 5011
                    PMM=-PMM*FACT*SOMX2;
                    FACT=FACT+2.0;
                end
            end
            if(L==M)
                PLGNDR=PMM;
            else
                PMMP1=X*(2*M+1)*PMM;
                if(L==M+1)
                    PLGNDR=PMMP1;
                else
                    for LL=M+2:L % loop 5012
                        PLL=(X*(2*LL-1)*PMMP1-(LL+M-1)*PMM)/(LL-M);
                        PMM=PMMP1;
                        PMMP1=PLL;
                    end
                    PLGNDR=PLL;
                end
            end  
            P(J,I2)=sqrt((2.0*L+1.0)/2.0)*PLGNDR;
            for  I1=(L-M+1):(L+M) % loop 5150
                FAC = I1;
                P(J,I2)=P(J,I2)/sqrt(FAC);
            end
        end
    end
    PLeg{M+1}=P;

%     JSTATE=0;
    % INITIALIZE AMPLITUDES
    for JSTATE=M+1:MMAX % loop 1000, site 30
    % Main loop for solving TDSE with population only on a single J-state
    % (--> JSTATE) as the initial condition. Evolution of the J-state
    % amplitudes is calculated for times within the laser pulse envelope
    % (--> T). For times after the envelope (--> TT), PROB and TCOS2 are
    % calculated for the angular probability distribution and its second
    % moment, respectively.
    %       WRITE(*,*) '|JM> =', (JSTATE-1),M, MMAX
    disp(['|J,M,MMAX>=' num2str(JSTATE-1) ' ' num2str(M) ' ' num2str(MMAX)]);
        MAXJ=0;
        for J=1:JMAX+2 % 50
            RC(J)=0.0;
            RCI(J) =0.0;
            IC(J)=0.0;
            ICI(J)=0.0;
            RT(J)=0.0;
            IT(J)=0.0;
            RRC(J)=0.0;
            IIC(J)=0.0;
        end
        RP=zeros([THETANUM,NT]);
        IP=zeros([THETANUM,NT]);
%         for  I1=1:THETANUM+1 % loop 51
%             for I3=1:NT+1 % loop 52
%                 RP(I1,I3) = 0.0;
%                 IP(I1,I3) = 0.0;
%             end
%         end
        RC(JSTATE) = 1.0;
        RCI(JSTATE)=1.0;
        RT(JSTATE)=1.0;
        TSTEP=TSTEP0;
        T=TSTART;
        TT=T;
        indTT=1;
%         disp([num2str(M) ' / ' num2str(T) ' / ' num2str(TT)]);
        FLAG = 0;
        JJ=1;
        OLDNORM=1.0;
        TTEST1 = 0.0;
        TTEST2 = 0.0;
        OLDTTEST1 = 0.0;
        OLDTTEST2 = 0.0;
        flag_exit2=0;

        while (flag_exit2==0) && (flag_exit1==0)
%             disp([num2str(M) ' / ' num2str(T) ' / ' num2str(TT)]);
            F=-FHCP*(exp(-(1.67*T/TAU)^2)-exp(-(1.67*(T/TAU-1.0)).^2))/0.95; % site 98
            ILASER=I0*exp(-(1.67*(T-DELAY)/TAUL).^2);
            TTEST1 = (atan(2.0*(T-DELAY)/TAUL)+pi/2.0)/pi;
            TTEST2 = (atan(1.0*T/TAU)+pi/2.0)/pi;
        % 99    CONTINUE   % site 99
            for J=M+1:JMAX % loop 100
                RC(J)= R0*F*V(J+1)*(-sin(W(J+1)*T)*RT(J+1)+cos(W(J+1)*T)*IT(J+1))*TSTEP+RCI(J) + ...
                	+ AD*ILASER*VJJ(J)*IT(J)*TSTEP + AD*ILASER*VJ2(J+2)*(-sin(W2(J+2)*T)*RT(J+2)+cos(W2(J+2)*T)*IT(J+2))*TSTEP;
                IC(J)= (-R0)*F*V(J+1)*(sin(W(J+1)*T)*IT(J+1)+cos(W(J+1)*T)*RT(J+1))*TSTEP+ICI(J) + ...
                    + (-AD)*ILASER*VJJ(J)*RT(J)*TSTEP + (-AD)*ILASER*VJ2(J+2)*(sin(W2(J+2)*T)*IT(J+2)+cos(W2(J+2)*T)*RT(J+2))*TSTEP;
                if J~=1
                    RC(J)= RC(J) + R0*F*V(J)*(sin(W(J)*T)*RT(J-1)+cos(W(J)*T)*IT(J-1))*TSTEP;
                    IC(J)= IC(J) + (-R0)*F*V(J)*(-sin(W(J)*T)*IT(J-1)+cos(W(J)*T)*RT(J-1))*TSTEP;
                end
                if (J~=1) && (J~=2)
                    RC(J)=RC(J) + AD*ILASER*VJ2(J)*(sin(W2(J)*T)*RT(J-2)+cos(W2(J)*T)*IT(J-2))*TSTEP;
                    IC(J)=IC(J) + (-AD)*ILASER*VJ2(J)*(-sin(W2(J)*T)*IT(J-2)+cos(W2(J)*T)*RT(J-2))*TSTEP;
                end      
            end % end loop 100
        %       WRITE(*,*)((RC(J).^2+IC(J).^2), J=0,JMAX)
        %     CHECK NORMALIZATION
            NORM = 0.0;
            for J=M+1:JMAX % loop 200
                NORM=NORM+RC(J).^2+IC(J).^2;
            end
            coeff1_out(:,M+1)=map2colvec(RC+1i*IC);
            if T>-5e3
                1;
            end
            
        %      write(*,*)oldttest1,ttest1,oldttest2,ttest2
            if (abs(TTEST1-OLDTTEST1)>0.1) || (abs(TTEST2-OLDTTEST2)> 0.1)
                TSTEP=TSTEP/2.0; % site 203
                T=T-TSTEP/2.0;
%                 disp([num2str(M) ' / ' num2str(T) ' / ' num2str(TT)]);
                for K=M+1:JMAX+1 % loop 201
                     RT(K)=(RT(K)+RCI(K))/2.0;
                     IT(K)=(IT(K)+ICI(K))/2.0;
                end
            %	  WRITE(*,*) TSTEP
                if (TSTEP < 1E-9)
                     disp( 'TIME STEP IS TOO SMALL');
                     flag_exit1=1;
                end
                continue;
            else
                if (abs(NORM-OLDNORM)>= ACCUR)
                    TSTEP=TSTEP/2.0; % site 203
                    T=T-TSTEP/2.0;
%                     disp([num2str(M) ' / ' num2str(T) ' / ' num2str(TT)]);
                    for K=M+1:JMAX+1 % loop 201
                         RT(K)=(RT(K)+RCI(K))/2.0;
                         IT(K)=(IT(K)+ICI(K))/2.0;
                    end
                %	  WRITE(*,*) TSTEP
                    if (TSTEP < 1E-9)
                         disp( 'TIME STEP IS TOO SMALL');
                         flag_exit1=1;
                    end
                    continue;
                end
            end
            
        %	  PREPARTE FOR NEXT STEP
            OLDNORM=NORM; % site 250
            OLDTTEST1=TTEST1;
            OLDTTEST2=TTEST2;
            for J=M+1:JMAX+1 % loop 300
                RT(J)=(2.0*RC(J)-RCI(J));
                IT(J)=(2.0*IC(J)-ICI(J));
                RCI(J)=RC(J);
                ICI(J)=IC(J);
            end
            %       WRITE(*,*)((RCI(J).^2+ICI(J).^2), J=0,JMAX),(TT-T)  
            %	  CALCULATE CONTRIBUTION TO cos(THETA)
            if ~((T > 4.5*TAU) && (T > 3.5*TAUL+DELAY))
                %       WRITE(*,*)((RCI(J).^2+ICI(J).^2), J=M,JMAX),(TT-T)
                if (T >= TT)
                    for J=M+1:JMAX+1 % loop 305
                %	     	  WRITE(*,*) (RCI(J).^2+ICI(J).^2),(TT-T)
                %	     	  RRC(J)=(RT(J)-RCI(J))*(TT-T)/TSTEP+RCI(J)
                %	     	  IIC(J)=(IT(J)-ICI(J))*(TT-T)/TSTEP+ICI(J)
                    RRC(J)=2.0*RCI(J)-RT(J);
                    IIC(J)=2.0*ICI(J)-IT(J);
                %	     	  WRITE(*,*) (RCI(J).^2+ICI(J).^2),(TT-T)
                %	  CALCULATE PSI(THETA)	 
                    end
                    statecoeffs(:,M+1)=map2colvec(RRC+1i*IIC);
                %	  WRITE(*,*)((RRC(J).^2+IIC(J).^2),J=M,JMAX),JJ
                    for J=M+1:JMAX+1 % loop 310
                        PIECE2=VJJ(J)*AMP(JSTATE)*(RRC(J).^2+IIC(J).^2);         
                        if (J==0) || (J>2)
                            PIECE=2.0*V(J)*AMP(JSTATE)*((RRC(J)*RRC(J-1)+IIC(J)*IIC(J-1))*cos(W(J)*TT)-(IIC(J-1)*RRC(J)-IIC(J)*RRC(J-1))*sin(W(J)*TT));
                            TCOS(JJ)=TCOS(JJ)+2.0*PIECE;
                            if (M == 0)
                                TCOS(JJ)=TCOS(JJ)-PIECE;
                            end
                            PIECE2=2.0*VJ2(J)*AMP(JSTATE)*((RRC(J)*RRC(J-2)+IIC(J)*IIC(J-2))*cos(W2(J)*TT)-(IIC(J-2)*RRC(J)-IIC(J)*RRC(J-2))*sin(W2(J)*TT))+PIECE2;
                        end
                        TCOS2(JJ)=TCOS2(JJ)+2.0*PIECE2; % site 306
                        if (M == 0)
                            TCOS2(JJ)=TCOS2(JJ)-PIECE2;	                       
                        end
                        JNORM=RRC(J).^2+IIC(J).^2;
                        if (JNORM > 1E-6)
                          MAXJ=J;
                        end
                        for I1=1:THETANUM % loop 333
                            RP(I1,JJ)=(RRC(J)*cos(E(J)*TT)+IIC(J)*sin(E(J)*TT))*P(J,I1)+RP(I1,JJ);
                            IP(I1,JJ)=(-RRC(J)*sin(E(J)*TT)+IIC(J)*cos(E(J)*TT))*P(J,I1)+IP(I1,JJ);
%                             RP(I1,JJ,M)=(RRC(J)*cos(E(J)*TT)+IIC(J)*sin(E(J)*TT))*P(J,I1);
%                             IP(I1,JJ,M)=(-RRC(J)*sin(E(J)*TT)+IIC(J)*cos(E(J)*TT))*P(J,I1);

                        end   	      		  
                    end
                    if (MAXJ == JMAX)
                        disp('J VALUES POPULATED BEYOND JMAX.');
                        flag_exit1=1;
                    end
                    if (MAXJ > MMAXJ)
                      MMAXJ = MAXJ;
                    end
                    TT=TT+DT;
                    indTT=indTT+1;
%                     disp([num2str(M) ' / ' num2str(T) ' / ' num2str(TT)]);
                    JJ=JJ+1;
                end

                if (DELAY >= 0.0)	% site 350
                    if (T < DELAY)
                        if (TSTEP > TAUL/100.)
                            T = T+TSTEP;
%                             disp([num2str(M) ' / ' num2str(T) ' / ' num2str(TT)]);
                            continue;
                        end
                    end
                else
                    if (T < 0.0)
                        if (TSTEP > TAU/100)
                            T = T+TSTEP;
%                             disp([num2str(M) ' / ' num2str(T) ' / ' num2str(TT)]);
                            continue;
                        end
                    end
                end
                T=T+1.5*TSTEP;
                TSTEP=2.0*TSTEP;
%                 disp([num2str(M) ' / ' num2str(T) ' / ' num2str(TT)]);
            %	  WRITE(*,*)T,TSTEP
            else
                flag_exit2=1; % get out of while loop
            end
        end
        JJMAX=JJ; % site 500
        for JJ=JJMAX:NT % loop 550
            for J=M+1:JMAX+1 % loop 540
                PIECE2=VJJ(J)*AMP(JSTATE)*(RRC(J).^2+IIC(J).^2);
                if (J==0) || (J>2)
                    PIECE=2.0*V(J)*AMP(JSTATE)*((RRC(J)*RRC(J-1)+IIC(J)*IIC(J-1))*cos(W(J)*TT)-(IIC(J-1)*RRC(J)-IIC(J)*RRC(J-1))*sin(W(J)*TT));
                    %              if (J.EQ.1) WRITE(*,*) PIECE,W(J),TT*2.42E-5
                    TCOS(JJ)=TCOS(JJ)+2.0*PIECE;
                    if (M == 0)
                        TCOS(JJ)=TCOS(JJ)-PIECE;
                    end
                    PIECE2=2.0*VJ2(J)*AMP(JSTATE)*((RRC(J)*RRC(J-2)+IIC(J)*IIC(J-2))*cos(W2(J)*TT)-(IIC(J-2)*RRC(J)-IIC(J)*RRC(J-2))*sin(W2(J)*TT))+PIECE2;
                end

                TCOS2(JJ)=TCOS2(JJ)+2.0*PIECE2; % site 501
                if (M == 0)
                    TCOS2(JJ)=TCOS2(JJ)-PIECE2;
                end
                for I1=1:THETANUM % loop 555
                    RP(I1,JJ)=(RRC(J)*cos(E(J)*TT)+IIC(J)*sin(E(J)*TT))*P(J,I1)+RP(I1,JJ);
                    IP(I1,JJ)=(-RRC(J)*sin(E(J)*TT)+IIC(J)*cos(E(J)*TT))*P(J,I1)+IP(I1,JJ);
%                     RP(I1,JJ,M)=(RRC(J)*cos(E(J)*TT)+IIC(J)*sin(E(J)*TT))*P(J,I1);
%                     IP(I1,JJ,M)=(-RRC(J)*sin(E(J)*TT)+IIC(J)*cos(E(J)*TT))*P(J,I1);
                end                
            end
        TT = TT+DT;
        indTT=indTT+1;
%         disp([num2str(M) ' / ' num2str(T) ' / ' num2str(TT)]);
        %     	  WRITE(*,*)TT*2.42E-5     	  
        end
        for J=1:JMAX+1 % loop 560
            if (J <= MAXJ)
%                 disp([num2str(RRC(J).^2+IIC(J).^2) ' ' num2str(J-1) ' ' num2str(MAXJ)]);
            end
            TOTNORM = AMP(JSTATE)*(RRC(J).^2+IIC(J).^2)+TOTNORM;
            if (M > 0)
                TOTNORM = TOTNORM+AMP(JSTATE)*(RRC(J).^2+IIC(J).^2);
            end
        end
        for I1=1:THETANUM % loop 770
            for I3=1:NT % loop 780
                PROB(I1,I3)=2.*AMP(JSTATE)*(RP(I1,I3).^2+IP(I1,I3).^2)+PROB(I1,I3);
                if M == 0
                    PROB(I1,I3)=PROB(I1,I3)-AMP(JSTATE)*(RP(I1,I3).^2+IP(I1,I3).^2);
                end
            end 
        end     
    end  % end of loop 1000
    M=M+1;
%     disp([num2str(M) '/' num2str(MMAX)]);
end
T=(TSTART:DT:(TSTART+(NT-1)*DT))*2.42e-5;
% for JJ=1:NT % loop 1500
%     T(JJ)=((JJ-1)*DT+TSTART)*2.42E-5;
% %     WRITE(1,*) T,TCOS(JJ), TCOS2(JJ)
% %     WRITE(2,*) T,(PROB(I,JJ),I=1,THETANUM+1)
% end

disp(['The final normalization is: ' num2str(TOTNORM)]);
output.T=T;
output.theta=THETA;
output.TCOS=TCOS;
output.TCOS2=TCOS2;
output.PROB=PROB;
output.PLeg=PLeg;
output.AMP=AMP;
output.RP=RP;
output.IP=IP;
output.RRC=RRC;
output.IIC=IIC;
output.statecoeffs=statecoeffs;
output.coeff1_out=coeff1_out;
end


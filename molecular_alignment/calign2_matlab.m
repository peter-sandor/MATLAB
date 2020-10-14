function [output] = calign2_matlab(input)
    
i1=input.FHCP/5.142E6; % 'ENTER THE PEAK HCP FIELD (1 KV/CM).'
taul1=input.TAU/2.42E-5; % 'ENTER THE HCP DURATION (GAUSSIAN FWHM) IN PSEC.'   
i0=input.I0/3.5E16/4.0; % 'ENTER THE PEAK LASER INTENSITY IN (W/CM.^2)' and CONVERT I0 TO -E**2/4
taul=input.TAUL/2.42E-5; % 'ENTER THE LASER PULSE DURATION (GAUSSIAN FWHM) IN PSEC.'
delay = input.DELAY/2.42E-5; % 'ENTER THE RELATIVE DELAY OF THE LASER PULSE IN PSEC.'
ad=input.AD; % 'ENTER THE POLARIZABILITY ANISOTROPY IN (A.U.)'
R0=input.R0; % 'ENTER THE MOLECULAR DIPOLE MOMENT IN (A.U.)'
b=input.B/(2.0*109737.0); % 'ENTER ROTATIONAL CONSTANT IN 1/CM.'
d=input.D/(2.0*109737.0); % 'ENTER CENTRIFUGAL DISTORTION CONSTANT IN 1/CM.'
% 'ENTER THE RELATIVE ABUNDANCE FACTORS FOR EVEN AND ODD J STATES'
even=input.EVEN;
odd=input.ODD;
jmax=input.JMAX; % 'ENTER THE MAXIMUM NUMBER OF J STATES TO CONSIDER.  '
if jmax > bigj
    disp('MAX J EXCEEDS CURRENT ARRAY LIMITS')
    break;
end    

%    	READ(1,*)TSTEP
tstep0=1.0/(20.0*b*jmax); % 'ENTER THE INITIAL TIME STEP IN PSEC.  '
tstep = tstep0;
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
tk=input.TK; % 'ENTER THE ROTATIONAL TEMPERATURE IN K.  '
tmax=input.TMAX/2.42E-5; % 'ENTER THE MAX DELAY IN PSEC.
dt=input.DT/2.42E-5; % 'ENTER DELAY BETWEEN DATA POINTS IN PSEC.'
NT=round((tmax-tstart)/dt);
if (NT > BIGT)
    disp('MAX TIME EXCEEDS CURRENT ARRAY LIMITS.');
end
accur=input.ACCUR; % 'ENTER THE SINGLE STEP NORMALIZATION ACCURACY.  '
thetanum=input.THETANUM; % 'ENTER THE NUMBER OF THETA STEPS FOR DENSTY PLOT.'  
thetastep=pi/(THETANUM-1);

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


/* This is a direct translation from Bob Jones' Align-2 code */




	bigj=100;
	bigt=1500;
	bigtheta=200;
	totnorm=0.;
	mmaxj=0;
% /*******************************************************\
% *		Beginning of File Input							*
% \*******************************************************/

	fscanf(fin,"%s", &filename1);
	fscanf(fin,"%s", &dummy);  /* ignore the descriptive string in the DAT-file */
	fscanf(fin,"%s", &filename2);
	fscanf(fin,"%s", &dummy);
	
	fscanf(fin,"%lf", &i0);
	fscanf(fin,"%lf", &i1);
	i0=-i0/3.5e16/4.0; /* convert to -E^2 */
	i1=-i1/3.5e16/4.0;
	fscanf(fin,"%s", &dummy);

	fscanf(fin,"%lf", &taul);
	fscanf(fin,"%lf", &taul1);

	taul/=2.42e-5;
	taul1/=2.42e-5;
	fscanf(fin,"%s", &dummy);

	fscanf(fin,"%lf", &delay);

	delay/=2.42e-5;
	fscanf(fin,"%s", &dummy);

	fscanf(fin,"%lf", &ad);

	fscanf(fin,"%s", &dummy);

	fscanf(fin,"%lf", &b);
	fscanf(fin,"%lf", &d);

	b/=(2.0*109737);
	d/=(2.0*109737);
	fscanf(fin,"%s", &dummy);

	fscanf(fin,"%lf", &even);
	fscanf(fin,"%lf", &odd);

	fscanf(fin,"%s", &dummy);

	fscanf(fin,"%d", &jmax);

	if(jmax>bigj){
		printf("Max J exceeds current array limits.\n");
		return 0;
	}
	tstep0=1.0/(40.0*b*jmax);
	tstep=tstep0;
	fscanf(fin,"%s", &dummy);

	fscanf(fin,"%lf", &tk);

	tk*=1.38e-23;
	tk/=1.6e-19;
	tk/=27.2;
	tk/=b;
	fscanf(fin,"%s", &dummy);
	
	fscanf(fin,"%lf", &tmax);
	fscanf(fin,"%lf", &dt);

	tmax/=2.42e-5;
	dt/=2.42e-5;
	maxt= (int)((tmax+3.5*taul)/dt);
	if(maxt>bigt){
		printf("Max time exceeds current array limits.\n");
		return 0;
	}


	thetanum=bigtheta;
	thetastep=PI/thetanum;
	for(jj=0;j<=maxt;j++){
		tcos2[jj]=0.0;
		for(i=0;i<=thetanum;i++){
			prob[i][jj]=0.0;
		}
	}
	theta[0]=0.0;
	for(i=1;i<=thetanum;i++){
		theta[i]=theta[i-1]+thetastep;
	}
	sumamp=0.0;
	e[0]=0.0;
	w2[0]=0.0;
	w2[1]=0.0;

	for(j=1;j<=jmax;j++){
		e[j]=(j+1)*j*(b-d*j*(j+1));
		if(j!=1){
			w2[j]=e[j]-e[j-2];
		}
		amp[j]=(2.0*j+1.0)*exp(-j*(j+1)/tk)/tk;
		sumamp+=amp[j];
	}
	amp[0]=1.0-sumamp;
	if(amp[0]<0.0) amp[0]=0.0;
	if(amp[0]>1.0) amp[0]=1.0;
	printf("%lf   %lf\n", amp[0], sumamp);
	bigamp=amp[0];
	ampj=0;
	for(j=0;j<=jmax;j++){
		amp[j]= amp[j]/(2.0*j+1.0);
		if((j/2.0)-(int)(j/2.0) < 0.1){
			amp[j]=2.0*even*amp[j]/(even+odd);
			printf("Even %d \n", j);
		}
		else{
			amp[j]=2.0*odd*amp[j]/(even+odd);
			printf("Odd %d \n", j);
		}
		if(amp[j]> bigamp){
			bigamp=amp[j];
			ampj=j;
		}
	}
	mmax=jmax; /* initialize mmax value before looping */
	for(j=ampj; j<=jmax;j++){
		if(amp[j] < (bigamp/1000.0)){
			mmax=j;
			break;
		}
	}
	printf("Max Thermally Populated J State is %d\n", mmax);

	rc[0]=0.0;
	ic[0]=0.0;
	rt[0]=0.0;
	it[0]=0.0;
	ici[0]=0.0;
	rci[0]=0.0;
	rrc[0]=0.0;
	iic[0]=0.0;
	m=0;
% /*******************************************************\
% *		Start of big do() loop							*
% \*******************************************************/
	do{
		printf("M= %d\n", m);
		for(j=m;j<=jmax+1;j++){
			vjj[j]=(2.0*j*(j+1)-2.0*m*m-1.0)/((2.0*j+3.0)*(2.0*j-1.0));
			vj2[j]=sqrt((j*j-m*m)*(pow((j-1), 2)-m*m)/((2.0*j+1.0)*(2.0*j-1.0)*(2.0*j-1.0)*(2.0*j-3.0)));
			/* insert legendre stuff */ 
			l=j;
			for(ii=0; ii<=thetanum;ii++){
				x=cos(theta[ii]);
				p[j][ii]=0.0;
				if(m<0 || m>l || fabs(x)>1) break;
				pmm=1.;
				if(m>0){
					somx2=sqrt((1.-x)*(1.+x));
					fact=1.;
					for(i=0; i<m;i++){
						pmm=-pmm*fact*somx2;
						fact=fact+2.;
					}
				}
				if(l==m){
					plgndr=pmm;
				}
				else{
					pmmp1=x*(2*m+1)*pmm;
					if(l == (m+1)){
						plgndr=pmmp1;
					}
					else{
						for(ll=m+2;ll<=l; ll++){
							pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
							pmm=pmmp1;
							pmmp1=pll;
						}
						plgndr=pll;
					}
				}
				p[j][ii]=sqrt((2.0*l+1.0)/2.0)*plgndr;
				for(i=(l-m+1);i<=(l+m);i++){
					fac=i;
					p[j][ii]/=sqrt(fac);
				}
			}
			/* end legendre stuff */
		}
% 	/*******************************************************\
% 	*		Start of big for() loop							*
% 	\*******************************************************/
		for(jstate=m;jstate<mmax;jstate++){
			printf("|JM> = %d %d %d\n", jstate, m, mmax);
			maxj=0;
			for(j=0;j<=jmax+1;j++){
				rc[j]=0.0;
				rci[j]=0.0;
				ic[j]=0.0;
				ici[j]=0.0;
				rt[j]=0.0;
				it[j]=0.0;
				rrc[j]=0.0;
				iic[j]=0.0;
			}
			for(i=0;i<=thetanum;i++){
				for(iii=0;iii<=maxt;iii++){
					rp[i][iii]=0.0;
					ip[i][iii]=0.0;
				}
			}
			rc[jstate]=1.0;
			rci[jstate]=1.0;
			rt[jstate]=1.0;
			tstep=tstep0;
			t=-3.5*taul+tstep/2.0;
			tt=t;
			jj=0;
			oldnorm=1.0;
% 		/**************** Start of Part A ******************/
			while(1){  /*infinite loop...break from this manually */
				check=fabs((t-delay)/taul1);
				ilaser=0.0;
				if ((t/taul) < 5.0) ilaser=i0*exp(-(1.67*t/taul)*(1.67*t/taul));
				if (check < 5.0) ilaser=ilaser+i1*exp(-(1.67*check)*(1.67*check));
				for(j=m;j<jmax;j++){
					rc[j]=ad*ilaser*vjj[j]*it[j]*tstep+rci[j];
					rc[j]=ad*ilaser*vj2[j+2]*(-sin(w2[j+2]*t)*rt[j+2]+cos(w2[j+2]*t)*it[j+2])*tstep+rc[j];
					ic[j]=-ad*ilaser*vjj[j]*rt[j]*tstep+ici[j];
					ic[j]=-ad*ilaser*vj2[j+2]*(sin(w2[j+2]*t)*it[j+2]+cos(w2[j+2]*t)*rt[j+2])*tstep+ic[j];
					if(j!=0){
						if(j!=1){
							rc[j]=ad*ilaser*vj2[j]*(sin(w2[j]*t)*rt[j-2]+cos(w2[j]*t)*it[j-2])*tstep+rc[j];
							ic[j]=-ad*ilaser*vj2[j]*(-sin(w2[j]*t)*it[j-2]+cos(w2[j]*t)*rt[j-2])*tstep+ic[j];
						}
					}
% /*100 */
				}
				norm=0.0;
				for(j=m;j<jmax;j++){
					norm=norm+rc[j]*rc[j]+ic[j]*ic[j];
				}
		
				if(jj==840)printf("norm is %le \n");
				if(fabs(norm-oldnorm) >= accur){
					tstep=tstep/2.0;
					t=t-tstep/2.0;
					for(k=m;k<=jmax;k++){
						rt[k]=(rt[k]+rci[k])/2.0;
						it[k]=(it[k]+ici[k])/2.0;
					}
					if(tstep < 1e-10){
						printf("Time step is too small %lf %lf\n", tstep, norm-oldnorm);
						return 0;
					}
					continue;
				}
				oldnorm=norm;
				for(j=m;j<=jmax;j++){
					rt[j]=(2.0*rc[j]-rci[j]);
					it[j]=(2.0*ic[j]-ici[j]);
					rci[j]=rc[j];
					ici[j]=ic[j];
				}
				if (t > (3.5*taul1+delay)) break;
				if (t >= tt){
					for(j=m;j<=jmax;j++){
						rrc[j]=2.0*rci[j]-rt[j];
						iic[j]=2.0*ici[j]-it[j];
% 	/* debug */ /*				if(j==m) printf("sum squared is %.8lf\n", rrc[j]*rrc[j]+iic[j]*iic[j]); */
					}
					
					for(j=m;j<=jmax;j++){
						piece2=vjj[j]*amp[jstate]*(rrc[j]*rrc[j]+iic[j]*iic[j]);
						if(j != 0){
							if(j != 1){
								piece2=2.0*vj2[j]*amp[jstate]*((rrc[j]*rrc[j-2]+iic[j]*iic[j-2])*cos(w2[j]*tt)-(iic[j-2]*rrc[j]-iic[j]*rrc[j-2])*sin(w2[j]*tt))+piece2;
							}
						}
						tcos2[jj]=tcos2[jj]+2.0*piece2;
						if(m==0) tcos2[jj]-=piece2;
						jnorm=rrc[j]*rrc[j]+iic[j]*iic[j];
						if(jnorm>1e-6) {
							maxj=j;
% 	/*						printf("new max j is %d\n", maxj); */
						}
						for(i=0;i<=thetanum;i++){
							rp[i][jj]=(rrc[j]*cos(e[j]*tt)+iic[j]*sin(e[j]*tt))*p[j][i]+rp[i][jj];
							ip[i][jj]=(-rrc[j]*sin(e[j]*tt)+iic[j]*cos(e[j]*tt))*p[j][i]+ip[i][jj];
						}
					}
% 			/*		printf("maxj is %d\n", maxj); */
					if(maxj == jmax){
						printf("J values populated beyond jmax\n");
						return 0;
					}
					if(maxj > mmaxj){
						 mmaxj=maxj;
% 	/*Debug*/ /*				printf("mmaxj is %d\n", mmaxj); */
					}
					tt+=dt;
					jj+=1;
				}
				if(t < delay){
					if(tstep > taul1/100.){
						t+=tstep;
						continue;
					}
				}
				t+=1.5*tstep;
				tstep*=2.0;		
			}
		/****************  End of Part A  ******************/
		/**************** Start of Part B ******************/
			jjmax=jj;
			for(jj=jjmax;jj<maxt+1;jj++){
				if(jj==840){
% /*					 printf("tcos2[jj] initial is %lf\t jmax is %d\t m is %d\ttt is %lf\n", tcos2[jj], jmax,m,tt); */
				}
				counter=0;
				for(j=m;j<=jmax;j++){
					piece2=vjj[j]*amp[jstate]*(rrc[j]*rrc[j]+iic[j]*iic[j]);
					if(j!=0){
						if(j!=1){
							piece2=2.0*vj2[j]*amp[jstate]*((rrc[j]*rrc[j-2]+iic[j]*iic[j-2])*cos(w2[j]*tt)-(iic[j-2]*rrc[j]-iic[j]*rrc[j-2])*sin(w2[j]*tt))+piece2;
						}
					}
% 	/*				if(j==2&& jj==840) printf("piece2 is %le\n", piece2); */
					
					tcos2[jj]=tcos2[jj]+2.0*piece2;
					if(m == 0){
						 tcos2[jj]=tcos2[jj]-piece2;
						counter++;					
					}
					for(i=0;i<=thetanum;i++){
						rp[i][jj]=(rrc[j]*cos(e[j]*tt)+iic[j]*sin(e[j]*tt))*p[j][i]+rp[i][jj];
						ip[i][jj]=(-rrc[j]*sin(e[j]*tt)+iic[j]*cos(e[j]*tt))*p[j][i]+ip[i][jj];
					}
				}
% /*				if(jj==840) printf("counter is %d\n", counter); */
% /*				if(jj==840) printf("tcos2[jj] is %lf\n", tcos2[jj]); */
				tt+=dt;
			}
			for(j=0;j<=jmax;j++){
% /*				if(j==0) printf("maxj is %d\n", maxj);   */
				if(j < maxj) printf("%lf\t%d\t%d\n", (rrc[j]*rrc[j]+iic[j]*iic[j]), j, maxj);
				totnorm=amp[jstate]*(rrc[j]*rrc[j]+iic[j]*iic[j])+totnorm;
				if(m>0) totnorm=totnorm+amp[jstate]*(rrc[j]*rrc[j]+iic[j]*iic[j]);
			}
% /*Debug*/ /*		printf("totnorm is %.9lf\n", totnorm); */
			for(i=0;i<=thetanum;i++){
				for(iii=0; iii<=maxt;iii++){
					prob[i][iii]=2.*amp[jstate]*(rp[i][iii]*rp[i][iii]+ip[i][iii]*ip[i][iii])+prob[i][iii];
					if(m==0) prob[i][iii]=prob[i][iii]-amp[jstate]*(rp[i][iii]*rp[i][iii]+ip[i][iii]*ip[i][iii]);
				}	
			}
% 		/****************  End of Part B  ******************/
		}
% 	/*******************************************************\
% 	*		End of big for() loop							*
% 	\*******************************************************/
		m++;

	}while(m<mmax);
% /*******************************************************\
% *		End of big do() loop							*
% \*******************************************************/
	fout1=fopen(filename1, "w");
	fout2=fopen(filename2, "w");
	for(jj=0;jj<=maxt;jj++){
		t=(jj*dt-3.5*taul)*2.42e-5;
		fprintf(fout1, "%lf\t%lf\n", t, tcos2[jj]);
		fprintf(fout2, "%lf\t", t);
		for i=0:thetanum;i++{
			fprintf(fout2, "%lf\t", prob[i][jj]);
		}
		fprintf(fout2, "\n");
	}









output.T=T;
output.TCOS=TCOS;
output.TCOS2=TCOS2;
output.PROB=PROB;
output.PLeg=PLeg;
output.AMP=AMP;
output.RP=RP;
output.IP=IP;
output.RRC=RRC;
end


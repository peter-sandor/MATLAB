function [COSVJJ,COSVJ1,COS2VJJ,COS2VJ1,COS2VJ2] = generate_D_coupling_Bob(maxJ)

NJ=maxJ+1;
COSVJJ = zeros([NJ 2*NJ-1 2*NJ-1]); % cos(theta) coupling, J to J
COSVJ1 = zeros([NJ 2*NJ-1 2*NJ-1]); % cos(theta) coupling, J to J+1
COS2VJJ = zeros([NJ 2*NJ-1 2*NJ-1]); % cos^2(theta) coupling, J to J
COS2VJ1 = zeros([NJ 2*NJ-1 2*NJ-1]); % cos^2(theta) coupling, J to J+1
COS2VJ2 = zeros([NJ 2*NJ-1 2*NJ-1]); % cos^2(theta) coupling, J to J+2
Jvec=0:maxJ;
Mvec=-maxJ:maxJ;
Kvec=Mvec;
for indJ=1:NJ
    J=Jvec(indJ);
    for indM=1:length(Mvec)
        M=Mvec(indM);
        for indK=1:length(Kvec)
            K=Mvec(indK);
            if (J>=abs(M)) && (J>=abs(K))
                COSVJJ(indJ,indM,indK)=M*K/((J+1)*(J+1E-20));
                COS2VJJ(indJ,indM,indK)=((2.0)*(3*K*K-(J+1)*J)*(3*M*M - (J+1)*J)/(3*(J+1E-20)*(J+1)*(2*J+3)*(2*J-1))) + (1.0/3.0);
%                 J=Jvec(indJ)+1;
                COSVJ1(indJ,indM,indK)=sqrt(((J+1)^2-M*M)*((J+1)^2-K*K)/((2.0*(J+1)+1.0)*(2.0*(J+1)-1.0)))/(J+1);
                COS2VJ1(indJ,indM,indK)=(2*M*K)/((J+2)*((J+1E-20)*(J+1)))*sqrt(((J+1)^2-K*K)*((J+1)^2-M*M)/((2.0*(J+1)+1.0)*(2.0*(J+1)-1.0)));
%                 J=Jvec(indJ)+2;
                FACTOR1 = sqrt(((J+2.0)^2-M*M)*((J+1)^2-M*M));
                FACTOR2 = sqrt(((J+2.0)^2-K*K)*((J+1)^2-K*K));
                FACTOR3 = sqrt(1.0/((2.0*(J+2)+1.0)*(2.0*(J+2)-1.0)*(2.0*(J+2)-1.0)*(2.0*(J+2)-3.0)));
                COS2VJ2(indJ,indM,indK)=FACTOR1*FACTOR2*FACTOR3/((J+2.0+1E-20)*(J+1+1E-20));
            end
        end
    end
end
end
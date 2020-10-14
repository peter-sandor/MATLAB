function [cJ1J,cJ1Jp1,cJ2J,cJ2Jp1,cJ2Jp2,varargout] = generate_D_coupling(maxJ)

cJ1J=zeros([maxJ+1 2*maxJ+1 2*maxJ+1]); % <J,M,K|1,0,0|J,M,K>
cJ2J=zeros([maxJ+1 2*maxJ+1 2*maxJ+1]); % <J,M,K|2,0,0|J,M,K>
cJ1Jp1=zeros([maxJ+1 2*maxJ+1 2*maxJ+1]); % <J,M,K|1,0,0|J+1,M,K>
cJ2Jp1=zeros([maxJ+1 2*maxJ+1 2*maxJ+1]); % <J,M,K|2,0,0|J+1,M,K>
cJ2Jp2=zeros([maxJ+1 2*maxJ+1 2*maxJ+1]); % <J,M,K|2,0,0|J+2,M,K>
c0=zeros([maxJ+1 2*maxJ+1 2*maxJ+1]); 
Jvals=0:maxJ;
Mvals=-maxJ:maxJ;
Kvals=Mvals;
for indJ=1:maxJ+1
    J=Jvals(indJ);
    for indK=1:2*maxJ+1
        K=Kvals(indK);
        for indM=1:2*maxJ+1 %Look to see if bi-section search algorithm can speed this up. 
            M=Mvals(indM);
    %         ind_offset=0;
            if (J>=abs(M)) && (J>=abs(K))
                c0(indJ,indM,indK)=1; %Creates a matrix of 1,0's where 1 represents a state combination that exist and 0 for one that doesn't exist. 
                Jp=J;
                if J==0 %populating J = 0. Missing a term of 1/3 for cos2(theta) for cJ2J?
                    cJ1J(indJ,indM,indK)=(-1)^(M-K)*sqrt(3/8*(2*J+1)*(2*Jp+1))/pi*Wigner3j([J 1 Jp],[-K 0 K])*Wigner3j([J 1 Jp],[-M 0 M]);
                    cJ2J(indJ,indM,indK)=(-1)^(M-K)*sqrt(5/8*(2*J+1)*(2*Jp+1))/pi*Wigner3j([J 2 Jp],[-K 0 K])*Wigner3j([J 2 Jp],[-M 0 M]);
                else
                    cJ1J(indJ,indM,indK)=(-1)^(M-K)*sqrt(3/8*(2*J+1)*(2*Jp+1))/pi*approx_Wigner3j([J 1 Jp],[-K 0 K])*approx_Wigner3j([J 1 Jp],[-M 0 M]);
                    cJ2J(indJ,indM,indK)=(-1)^(M-K)*sqrt(5/8*(2*J+1)*(2*Jp+1))/pi*approx_Wigner3j([J 2 Jp],[-K 0 K])*approx_Wigner3j([J 2 Jp],[-M 0 M]);
                end
                Jp=J+1;
                
                if J==0
                    cJ1Jp1(indJ,indM,indK)=(-1)^(M-K)*sqrt(3/8*(2*J+1)*(2*Jp+1))/pi*Wigner3j([J 1 Jp],[-K 0 K])*approx_Wigner3j([J 1 Jp],[-M 0 M]);
                    cJ2Jp1(indJ,indM,indK)=(-1)^(M-K)*sqrt(5/8*(2*J+1)*(2*Jp+1))/pi*Wigner3j([J 2 Jp],[-K 0 K])*Wigner3j([J 2 Jp],[-M 0 M]);
                else
                    cJ1Jp1(indJ,indM,indK)=(-1)^(M-K)*sqrt(3/8*(2*J+1)*(2*Jp+1))/pi*approx_Wigner3j([J 1 Jp],[-K 0 K])*approx_Wigner3j([J 1 Jp],[-M 0 M]);
                    cJ2Jp1(indJ,indM,indK)=(-1)^(M-K)*sqrt(5/8*(2*J+1)*(2*Jp+1))/pi*approx_Wigner3j([J 2 Jp],[-K 0 K])*approx_Wigner3j([J 2 Jp],[-M 0 M]);
                end
                Jp=J+2;
                cJ2Jp2(indJ,indM,indK)=(-1)^(M-K)*sqrt(5/8*(2*J+1)*(2*Jp+1))/pi*approx_Wigner3j([J 2 Jp],[-K 0 K])*approx_Wigner3j([J 2 Jp],[-M 0 M]);
            end
        end
    end
end
varargout{1}=c0;
end
function [cJ2J,cJ1Jp1,cJ2Jp1,cJ2Jp2] = generate_D_coupling_old(maxJ)

cJ2J=zeros([maxJ+1 2*maxJ+1 2*maxJ+1]); % <J,M,K|2,0,0|J,M,K>
cJ1Jp1=zeros([maxJ 2*maxJ+1 2*maxJ+1]); % <J,M,K|1,0,0|J+1,M,K>
cJ2Jp1=zeros([maxJ 2*maxJ+1 2*maxJ+1]); % <J,M,K|2,0,0|J+1,M,K>
cJ2Jp2=zeros([maxJ-1 2*maxJ+1 2*maxJ+1]); % <J,M,K|2,0,0|J+2,M,K>
% Jvals=0:maxJ;
Mvals=-maxJ:maxJ;
Kvals=Mvals;
for indK=1:2*maxJ+1
    K=Kvals(indK);
    for indM=1:2*maxJ+1
        M=Mvals(indM);
        Jvals=max([abs(M) abs(K)]):maxJ;
        ind_offset=max([abs(M) abs(K)]);
%         ind_offset=0;
        for indJ=1:length(Jvals)
            J=Jvals(indJ);
            Jp=J;
            if J==0
    %             cJ1J(indJ+indM-1,indM)=(-1)^m*sqrt(5/4/pi*(2*J+1)*(2*Jp+1))*Wigner3j([J 1 Jp],[0 0 0])*Wigner3j([J 1 Jp],[-m 0 m]);
                cJ2J(indJ+ind_offset,indM,indK)=(-1)^(M-K)*sqrt(5/8*(2*J+1)*(2*Jp+1))/pi*Wigner3j([J 2 Jp],[-K 0 K])*Wigner3j([J 2 Jp],[-M 0 M]);
            else
    %             cJ1J(indJ+indM-1,indM)=(-1)^m*sqrt(5/4/pi*(2*J+1)*(2*Jp+1))*approx_Wigner3j([J 1 Jp],[0 0 0])*approx_Wigner3j([J 1 Jp],[-m 0 m]);
                cJ2J(indJ+ind_offset,indM,indK)=(-1)^(M-K)*sqrt(5/8*(2*J+1)*(2*Jp+1))/pi*approx_Wigner3j([J 2 Jp],[-K 0 K])*approx_Wigner3j([J 2 Jp],[-M 0 M]);
            end
        end
        Jvals=max([abs(M) abs(K)]):maxJ-1;
        for indJ=1:length(Jvals)
            J=Jvals(indJ);
            Jp=J+1;
            cJ1Jp1(indJ+ind_offset,indM,indK)=(-1)^(M-K)*sqrt(3/8*(2*J+1)*(2*Jp+1))/pi*approx_Wigner3j([J 1 Jp],[-K 0 K])*approx_Wigner3j([J 1 Jp],[-M 0 M]);
            cJ2Jp1(indJ+ind_offset,indM,indK)=(-1)^(M-K)*sqrt(3/8*(2*J+1)*(2*Jp+1))/pi*approx_Wigner3j([J 2 Jp],[-K 0 K])*approx_Wigner3j([J 2 Jp],[-M 0 M]);
        end
        Jvals=max([abs(M) abs(K)]):maxJ-2;
        for indJ=1:length(Jvals)
            J=Jvals(indJ);
            Jp=J+2;
            cJ2Jp2(indJ+ind_offset,indM,indK)=(-1)^(M-K)*sqrt(5/8*(2*J+1)*(2*Jp+1))/pi*approx_Wigner3j([J 2 Jp],[-K 0 K])*approx_Wigner3j([J 2 Jp],[-M 0 M]);
        end
    end
end
end
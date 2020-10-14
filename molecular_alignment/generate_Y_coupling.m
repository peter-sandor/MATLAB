function [cJ2J,cJJp1,cJJp2] = generate_Y_coupling(maxJ)

% cJ1J=zeros([maxJ+1 maxJ+1]);
cJ2J=zeros([maxJ+1 maxJ+1]);
cJJp1=zeros([maxJ maxJ]);
cJJp2=zeros([maxJ-1 maxJ-1]);
mvals=0:maxJ;
for indM=1:maxJ+1
    m=mvals(indM);
    Jvals=abs(m):maxJ;
    % coefficients calculated: <J,M|2,0|J,M>, <J,M|1,0|J+1,M> and <J,M|2,0|J+2,M>.
    for indJ=1:length(Jvals)
        J=Jvals(indJ);
        Jp=J;
        if J==0
%             cJ1J(indJ+indM-1,indM)=(-1)^m*sqrt(5/4/pi*(2*J+1)*(2*Jp+1))*Wigner3j([J 1 Jp],[0 0 0])*Wigner3j([J 1 Jp],[-m 0 m]);
            cJ2J(indJ+indM-1,indM)=(-1)^m*sqrt(5/4/pi*(2*J+1)*(2*Jp+1))*Wigner3j([J 2 Jp],[0 0 0])*Wigner3j([J 2 Jp],[-m 0 m]);
        else
%             cJ1J(indJ+indM-1,indM)=(-1)^m*sqrt(5/4/pi*(2*J+1)*(2*Jp+1))*approx_Wigner3j([J 1 Jp],[0 0 0])*approx_Wigner3j([J 1 Jp],[-m 0 m]);
            cJ2J(indJ+indM-1,indM)=(-1)^m*sqrt(5/4/pi*(2*J+1)*(2*Jp+1))*approx_Wigner3j([J 2 Jp],[0 0 0])*approx_Wigner3j([J 2 Jp],[-m 0 m]);
        end
    end
    Jvals=abs(m):maxJ-1;
    for indJ=1:length(Jvals)
        J=Jvals(indJ);
        Jp=J+1;
        cJJp1(indJ+indM-1,indM)=(-1)^m*sqrt(3/4/pi*(2*J+1)*(2*Jp+1))*approx_Wigner3j([J 1 Jp],[0 0 0])*approx_Wigner3j([J 1 Jp],[-m 0 m]);
    end
    Jvals=abs(m):maxJ-2;
    for indJ=1:length(Jvals)
        J=Jvals(indJ);
        Jp=J+2;
        cJJp2(indJ+indM-1,indM)=(-1)^m*sqrt(5/4/pi*(2*J+1)*(2*Jp+1))*approx_Wigner3j([J 2 Jp],[0 0 0])*approx_Wigner3j([J 2 Jp],[-m 0 m]);
    end
end
end
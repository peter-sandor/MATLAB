function [VJJ,VJ1,VJ2] = generate_Y_coupling_Bob(maxJ)

NJ=maxJ+1;
VJ1=zeros([NJ NJ]);
VJJ=zeros([NJ NJ]);
VJ2=zeros([NJ NJ]);

Jvec=0:maxJ;
for indJ=1:NJ
    J=Jvec(indJ);
    Mvec=0:J;
    for indM=1:length(Mvec)
        M=Mvec(indM);
        VJJ(indJ,M+1)=(2.0*J*(J+1)-2.0*M^2-1.0)/((2.0*J+3.0)*(2.0*J-1.0));
        J=Jvec(indJ)+1;
        VJ1(indJ,M+1)=sqrt((J.^2-M^2)/((2.0*J+1.0)*(2.0*J-1.0)));
        J=Jvec(indJ)+2;
        VJ2(indJ,M+1)=sqrt((J.^2-M^2)*((J-1).^2-M^2)/((2.0*J+1.0)*(2.0*J-1.0)*(2.0*J-1.0)*(2.0*J-3.0)));
    end
end
end
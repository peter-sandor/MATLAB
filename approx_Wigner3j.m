function W = approx_Wigner3j(J123,M123)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if J123(2)==2 && M123(2)==0 && M123(3)==-M123(1)
    J=J123(1);
    M=M123(3);
    if J123(3)==J123(1)
        W=(-1)^(J-M)*2*(3*M^2-J*(J+1))/sqrt((2*J+3)*(2*J+2)*(2*J+1)*(2*J)*(2*J-1));
    elseif J123(3)==J123(1)+1
        W=(-1)^(J-M+1)*2*M*sqrt(6*(J+M+1)*(J-M+1)/((2*J+4)*(2*J+3)*(2*J+2)*(2*J+1)*2*J));
    elseif J123(3)==J123(1)+2
        W=(-1)^(J-M)*sqrt((6*(J+M+2)*(J+M+1)*(J-M+2)*(J-M+1))/((2*J+5)*(2*J+4)*(2*J+3)*(2*J+2)*(2*J+1)));
    elseif J123(3)==J123(1)-2
        W=(-1)^(J+M)*sqrt((6*(J-M)*(J-M-1)*(J+M)*(J+M-1))/((2*J+1)*2*J*(2*J-1)*(2*J-2)*(2*J-3)));
%         J=J123(3);
%         M=M123(1);
%         W=(-1)^(J-M)*sqrt((6*(J+M+2)*(J+M+1)*(J-M+2)*(J-M+1))/((2*J+5)*(2*J+4)*(2*J+3)*(2*J+2)*(2*J+1)));
    else
        W=[];
        disp('Ran into unhandled case of Wigner 3j symbol approximations.');
    end
elseif J123(2)==1 && M123(2)==0 && M123(3)==-M123(1)
    J=J123(1);
    M=M123(3);
    if J123(3)==J123(1)
        W=(-1)^(J-M)*M/sqrt((2*J+1)*(J+1)*J);
    elseif J123(3)==J123(1)+1
        W=(-1)^(J-M-1)*sqrt((2*(J+M+1)*(J-M+1))/((2*J+3)*(2*J+2)*(2*J+1)));
    elseif J123(3)==J123(1)-1
        W=(-1)^(J-M-2)*sqrt((2*(J+M)*(J-M))/((2*J+1)*2*J*(2*J-1)));
%         J=J123(3);
%         M=M123(1);
%         W=(-1)^(J-M)*sqrt((6*(J+M+2)*(J+M+1)*(J-M+2)*(J-M+1))/((2*J+5)*(2*J+4)*(2*J+3)*(2*J+2)*(2*J+1)));
    else
        W=[];
        disp('Ran into unhandled case of Wigner 3j symbol approximations.');
    end
else
    disp('input error');
end
end


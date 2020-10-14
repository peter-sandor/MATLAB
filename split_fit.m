function FIT_parts = split_fit(FIT,Nparts)

% Nparts=4;
Neq=length(FIT.eq);
Neq_part=(Neq-Nparts+1)/Nparts;
Nparams=length(FIT.m)/Nparts;
for ind1=1:Nparts
    if ind1==1
        % the first part of the equation is not offset by '+' sign
        FIT_parts(ind1).eq=FIT.eq(1+(ind1-1)*Neq_part:Neq_part+(ind1-1)*Neq_part);
    else
        % the subsequent parts of the equation are offset by (ind1-1) '+' signs
        FIT_parts(ind1).eq=FIT.eq(1+ind1-1+(ind1-1)*Neq_part:Neq_part+ind1-1+(ind1-1)*Neq_part);
    end
    FIT_parts(ind1).param=FIT.param(ind1:Nparts:ind1+(Nparams-1)*Nparts);
    FIT_parts(ind1).m0=FIT.m0(ind1:Nparts:ind1+(Nparams-1)*Nparts);
    FIT_parts(ind1).m=FIT.m(ind1:Nparts:ind1+(Nparams-1)*Nparts);
    FIT_parts(ind1).x=FIT.x;
    FIT_parts(ind1).y=FIT.y;
    FIT_parts(ind1).r=FIT.r;
    FIT_parts(ind1).xvar=FIT.xvar;
    FIT_parts(ind1).yvar=FIT.yvar;
    FIT_parts(ind1).fitmode=FIT.fitmode;
    FIT_parts(ind1).name=FIT.name;
end
end
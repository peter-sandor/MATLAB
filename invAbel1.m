function f0 = invAbel(F0)

%increase resolution
samples = 1000;
rax = 0:1/(samples-1):1;
F = interp1(0 : 1/(length(F0)-1) : 1, F0, rax,'line');


dF = gradient(F);

for indr = 1 : samples-1
    f(indr) = -sum(dF(indr+1:end) ./ sqrt( rax(indr+1 : end).^2 - rax(indr)^2 )); 
end;
f(samples) = 0;

%downscale resolution again
f0 = interp1(rax,f,0:1/(length(F0)-1):1,'line');






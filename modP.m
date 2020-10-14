function out = modP(A,B)
out=A-floor(A/B)*B;
out(out==0)=B;
end
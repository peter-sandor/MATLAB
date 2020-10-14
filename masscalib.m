function mass = masscalib(delay,A,t0)
mass=(((1:8000)-t0)/A).^2;
end
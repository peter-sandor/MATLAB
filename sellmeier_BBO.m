function n = sellmeier_BBO (lambda, p)

% BBO kristaly Sellmeier-egyenlete

L1_o = 0;
L2_o = 0.1338;
L3_o = 10;
A1_o = 0.7126;
A2_o = 1.0279;
A3_o = 1.535;

L1_e = 0;
L2_e = 0.1249;
L3_e = 10;
A1_e = 0.5525;
A2_e = 0.8205;
A3_e = 0.423;

%lambda = 0.4 : 0.01 : 1.1;

if p == 0
    n = sqrt (1 + A1_o .* lambda.^2 ./(lambda.^2 - L1_o^2) + A2_o .* lambda.^2 ./(lambda.^2 - L2_o.^2) + A3_o .* lambda.^2 ./(lambda.^2 - L3_o^2));
else if p == 1    
    n = sqrt (1 + A1_e .* lambda.^2 ./(lambda.^2 - L1_e^2) + A2_e .* lambda.^2 ./(lambda.^2 - L2_e.^2) + A3_e .* lambda.^2 ./(lambda.^2 - L3_e^2));
    else 'nem megfelelo "p" ertek!'
    end

end
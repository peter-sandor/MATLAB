function n = sellmeier_FS (lambda)

% az omlesztett kvarc toresmutatoja a h.hossz fuggvenyeben

% az RP Photonics honlapja alapjan:
% http://www.rp-photonics.com/sellmeier_formula.html
% [lambda] = mikron

% az ordinarius koefficiensek
A = 0.6961663;
B = 0.0684043^2;
C = 0.4079426;
D = 0.1162414^2;
E = 0.8974794;
F = 9.896161^2;

%lambda = 0.400 : 0.01 : 2.000;
n = sqrt (1 + A * lambda.^2 ./(lambda.^2 - B) + C * lambda.^2 ./(lambda.^2 - D) + E * lambda.^2 ./(lambda.^2 - F));

end
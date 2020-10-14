function v_group = group_vel_calcite(lambda, p)

% Kalcit kristalyban az ordinarius es extraordinarius csoportsebessegek a
% h.hossz fuggvenyeben
% http://www.redoptronics.com/Calcite-crystal.html
% [lambda] = mikrometer
% lambda: oszlopvektor!

% fenysebesseg vakuumban [mikron/fs]:
c = 0.299792;

% k_max = 100;
% delta_lambda = 0.001;
% lambda = (lambda_c-delta_lambda*(k_max/2-1) : delta_lambda : lambda_c+delta_lambda*k_max/2)';
delta_lambda = lambda(2)-lambda(1);
k_max = length(lambda);

if p == 0
    % az ordinarius koefficiensek
    A_o = 2.69705;
    B_o = 0.0192064;
    C_o = 0.01820;
    D_o = 0.0151624;
    n_o = sqrt (A_o + B_o ./(lambda.^2-C_o) - D_o.*lambda.^2);
    d_n_o(1) = (n_o(2) - n_o(1)) / delta_lambda;
        for k = 2 : (k_max - 1);
            d_n_o(k) = (n_o(k+1) - n_o(k-1)) / (2*delta_lambda);
        end
    d_n_o(k_max) = d_n_o(k_max - 1);
    n_group_o = n_o - lambda .* d_n_o';
    v_group_o = c ./ n_group_o; % [mikron/fs]
    v_group = v_group_o;

elseif p == 1
    % az extraordinarius koefficiensek
    A_e = 2.18438;
    B_e = 0.0087309;
    C_e = 0.01018;
    D_e = 0.0024411;
    n_e = sqrt (A_e + B_e ./(lambda.^2-C_e) - D_e.*lambda.^2);
    d_n_e(1) = ((n_e(2) - n_e(1)) / delta_lambda)';
    for k = 2 : (k_max - 1);
        d_n_e(k) = (n_e(k+1) - n_e(k-1)) / (2*delta_lambda);
    end
    d_n_e(k_max) = d_n_e(k_max - 1);
    n_group_e = n_e - lambda .* d_n_e';
    v_group_e = c ./ n_group_e; % [mikron/fs]
    v_group = v_group_e;
end
end

function [T,SCH_qp_ode45] = solve_SCH_quadforcesine(delay,init_cond,param)

    function FUN = func_SCH_quadpot(delay,P)
%         mass = 9.1e-31;
%         h_bar = 6.626e-34;
%         omega = 2*pi*7.5e8;
        
        % FUN = [P(2)/param(1); -param(1)*param(3)^2*P(1);-i*2*param(2)/param(1)*P(3)^2+i*param(1)/(2*param(2))*param(3)^2;1/(2*param(1))*P(2)^2-param(2)^2/(2*param(1))*P(3)-param(1)*param(3)^2/2*P(1)^2;delay_step];
        FUN = [P(2)/param(1); -param(1)*param(3)^2*P(1)+param(6)*sin(param(4)*delay+param(5));-i*2*param(2)/param(1)*P(3)^2+i*param(1)/(2*param(2))*param(3)^2;1/(2*param(1))*P(2)^2-param(2)^2/(2*param(1))*P(3)-param(1)*param(3)^2/2*P(1)^2+param(6)*sin(param(4)*delay+param(5));];
    end

[T, SCH_qp_ode45] = ode45(@func_SCH_quadpot,delay,init_cond);
end
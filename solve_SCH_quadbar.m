function [T,SCH_qp_ode45] = solve_SCH_quadbar(delay,init_cond,param)

    function FUN = func_SCH_quadbar(delay,P)
    FUN = [P(2)/param(1); param(1)*param(3)^2*P(1);-i*2*param(2)/param(1)*P(3)^2-i*param(1)/(2*param(2))*param(3)^2;3/2*P(2)^2+param(2)^2/(2*param(1))*P(3)+param(1)*param(3)^2/2*P(1)^2;];
    end

[T, SCH_qp_ode45] = ode45(@func_SCH_quadbar,delay,init_cond);
end
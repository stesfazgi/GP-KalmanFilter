function [theta_second_deriv,M,n] = StateVariablesHigherDerivatives(x,tau,exo_system)
% Returns [q_second_deriv, q_third_deriv, q_fourth_deriv, theta_second_deriv, delta, M, M_dot, M_dotdot, n, n_dotdot].
%     workspace;
    %% System parameters and convenient variable naming
    l = exo_system.l; 
    m = exo_system.m;
    r = exo_system.r;
    k = exo_system.k;
    d = exo_system.d; 
    j = exo_system.j; 
    Nminv = exo_system.Nminv;
    
    g = 9.81;
    lc1 = l(1)/2; lc2 = l(2)/2;

    I1 = 1/4*m(1)*r(1)^2 + 1/12*m(1)*l(1)^2;
    I2 = 1/4*m(2)*r(2)^2 + 1/12*m(2)*l(2)^2;    

    %% Calculate second derivative of link angles 
    D = [d(1), 0; 0, d(2)];
    K = [k(1), 0; 0, k(2)];
    J = diag([j(1),j(2)]);
    
    % Inertia
    M = [m(1)*lc1^2 + I1 + m(2)*(l(1)^2+lc2^2+2*l(1)*lc2*cos(x(5))) + I2, m(2)*(lc2^2+l(1)*lc2*cos(x(5))) + I2; 
                m(2)*(lc2^2+l(1)*lc2*cos(x(5))) + I2,                       m(2)*lc2^2 + I2           ];
    
    % Gravitational forces
    G = [m(1)*lc1*g*cos(x(1)) + m(2)*g*(lc2*cos(x(1)+x(5))+l(1)*cos(x(1)));
                            m(2)*g*lc2*cos(x(1)+x(5))];  
    
    % Coriolis forces
    C = [(-m(2)*l(1)*lc2*sin(x(5)))*x(6)^2 - 2*m(2)*l(1)*lc2*sin(x(5))*x(2)*x(6); 
                        m(2)*l(1)*lc2*sin(x(5))*x(2)^2];                    
    
    n = C + G;                                
    theta_second_deriv = inv(J)*tau - inv(J)*( Nminv*(D*(Nminv*[x(4);x(8)]-[x(2);x(6)]) + K*(Nminv*[x(3);x(7)]-[x(1);x(5)])) );
end
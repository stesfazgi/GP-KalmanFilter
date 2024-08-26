function dxdt = twoDofPlanarRobotWithDamping(t,x,tau,taud,exo_sys)
%
    %% Hilfsbenennungen    
    m1 = exo_sys.m(1); m2 = exo_sys.m(2);
    l1 = exo_sys.l(1); l2 = exo_sys.l(2);
    r1 = exo_sys.r(1); r2 = exo_sys.r(2);
    k1 = exo_sys.k(1); k2 = exo_sys.k(2);
    d1 = exo_sys.d(1); d2 = exo_sys.d(2);
    Nm = diag(exo_sys.Nm);
    Nminv = inv(Nm);
    j = exo_sys.j(1);
    
    g = 9.81;
    lc1 = l1/2; lc2 = l2/2;

    I1 = 1/4*m1*r1^2 + 1/12*m1*l1^2;
    I2 = 1/4*m2*r2^2 + 1/12*m2*l2^2;     
    
    %% Trägheitsmatrizen mit Ableitungen
    M = [m1*lc1^2 + I1 + m2*(l1^2+lc2^2+2*l1*lc2*cos(x(5))) + I2, m2*(lc2^2+l1*lc2*cos(x(5))) + I2; 
                    m2*(lc2^2+l1*lc2*cos(x(5))) + I2,                       m2*lc2^2 + I2           ];
    
    %% Coriolis-, Zentrifugalkräfte mit Ableitungen
    C = [(-m2*l1*lc2*sin(x(5)))*x(6)^2 - 2*m2*l1*lc2*sin(x(5))*x(2)*x(6); 
                            m2*l1*lc2*sin(x(5))*x(2)^2];
    
    %% Gravitation mit Ableitungen
    G = [m1*lc1*g*cos(x(1)) + m2*g*(lc2*cos(x(1)+x(5))+l1*cos(x(1)));
                            m2*g*lc2*cos(x(1)+x(5))];                        
                        
    % C und G kombinieren analog Paper
    n = C + G;
    
    %% Federmatrix und Trägheiten              
    K = [k1, 0;
         0, k2];
     
    D = [d1, 0; 
         0, d2];
    
    J = j * eye(2);
    
    %% human passive torques
    [Mh, tau_h] = human_arm_model(exo_sys, x);

    %% Build system 
    q_second_derivatives = -inv(M+Mh)*( n + D*([x(2);x(6)]-Nminv*[x(4);x(8)]) + K*([x(1);x(5)]-Nminv*[x(3);x(7)]) + tau_h - taud');
    theta_second_derivatives = -inv(J)*( Nminv*(D*(Nminv*[x(4);x(8)]-[x(2);x(6)]) + K*(Nminv*[x(3);x(7)]-[x(1);x(5)])) - tau);
    
    f_x = [x(2);q_second_derivatives(1);x(4);theta_second_derivatives(1);
           x(6);q_second_derivatives(2);x(8);theta_second_derivatives(2)];
    
%     inv_motor_inertia = inv(J);
%     g_x = [0, 0; 0, 0; 0, 0; inv_motor_inertia(1,1), inv_motor_inertia(1,2);
%            0, 0; 0, 0; 0, 0; inv_motor_inertia(2,1), inv_motor_inertia(2,2)];
%     
    dxdt = f_x;
%     dxdt = f_x+g_x*tau;   
end
function [Mh, tau_h] = human_arm_model(exo_sys, x)

    %% Hilfsbenennungen    
    m1 = 2.4; m2 = 2.0;
    l1 = exo_sys.l(1); l2 = exo_sys.l(2)+0.1;
    r1 = 0.06; r2 = 0.04;
    d1 = 0.65; d2 = 0.65;
    
    g = 9.81;
    lc1 = l1/2; lc2 = l2/2;

    I1 = 1/4*m1*r1^2 + 1/12*m1*l1^2;
    I2 = 1/4*m2*r2^2 + 1/12*m2*l2^2;     
    
    %% Trägheitsmatrizen mit Ableitungen
    Mh = [m1*lc1^2 + I1 + m2*(l1^2+lc2^2+2*l1*lc2*cos(x(5))) + I2, m2*(lc2^2+l1*lc2*cos(x(5))) + I2; 
                    m2*(lc2^2+l1*lc2*cos(x(5))) + I2,                       m2*lc2^2 + I2           ];
    
    %% Coriolis-, Zentrifugalkräfte mit Ableitungen
    C = [(-m2*l1*lc2*sin(x(5)))*x(6)^2 - 2*m2*l1*lc2*sin(x(5))*x(2)*x(6); 
                            m2*l1*lc2*sin(x(5))*x(2)^2];
    
    %% Gravitation mit Ableitungen
    G = [m1*lc1*g*cos(x(1)) + m2*g*(lc2*cos(x(1)+x(5))+l1*cos(x(1)));
                            m2*g*lc2*cos(x(1)+x(5))];                        
                        
    % C und G kombinieren analog Paper
    n = C + G;
    
    % damping
    D = diag([d1,d2])*[x(2) x(6)]';
    
    tau_h = n + D;
end


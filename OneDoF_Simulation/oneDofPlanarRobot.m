function dxdt = oneDofPlanarRobot(t,x,tau,tau_d,sys_params, perturbatesystem, mp, Mp, dh)
% Angepasst auf das device im Labor, Schwerkraft zeigt also in die falsche Richtung  
   
  if perturbatesystem
   m1 = mp;
   l1 = sys_params.l1;
   r1 = sys_params.r1;
   d1 = sys_params.d1;
   k1 = sys_params.k1;
   J = sys_params.j;
   M = Mp;
  else
   m1 = sys_params.m1;
   l1 = sys_params.l1;
   r1 = sys_params.r1;
   d1 = sys_params.d1;
   k1 = sys_params.k1;
   J = sys_params.j;
   M = 1.32e-2;
  end
   g = 9.81;
   
%    M = m1*(l1/2)^2 + 1/4*m1*r1^2+1/12*m1*l1^2;

%    f_x = [x(2);
%           -1/M * ( m1*g*l1/2*sin(x(1)) + d1*(x(2)-x(4)) + k1*(x(1)-x(3)) );
%           x(4);
%           -1/J * ( d1*(x(4)-x(2)) + k1*(x(3)-x(1)) )];
   g_x = [0; 0; 0; 1/J];
   g_d_x = [0; 1/M; 0; 0];
   
   f_x = [x(2);
          -1/M * ( -m1*g*l1/2*cos(x(1)) + d1*(x(2)-x(4)) + k1*(x(1)-x(3)) + dh*x(2));
          x(4);
          -1/J * ( d1*(x(4)-x(2)) +  k1*(x(3)-x(1)) )];
%    dxdt = f_x;
   dxdt = f_x + g_x*tau + g_d_x * tau_d;
end
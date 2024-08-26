function exo_system = config()
% clear all;
% close all;
% pack all;
% clc;
% 
% addpath('aux_scripts')

%% Parameters for the system
exo_system = robot_system();
exo_system.dofs = 1;
exo_system.l = [0.3];
exo_system.m = [0.459];
exo_system.r = [1];
exo_system.k = [301.5];
exo_system.d = [0];
exo_system.j = [4.4e-4];
exo_system.Nm = [1];

exo_system.lh = [1];
exo_system.mh = [1];
exo_system.rh = [1];


exo_system.distorder = 4;
S = [zeros(exo_system.distorder-1,1) eye(exo_system.distorder-1);
    0 zeros(1,exo_system.distorder-1)];
H = [1, zeros(1,exo_system.distorder-1)];

% Output matrix C (both position and velocity measurements are used)
C = [eye(exo_system.dofs*4), zeros(exo_system.dofs*4, exo_system.distorder*exo_system.dofs)];
exo_system.C = C;

exo_system.S = S;
exo_system.H = H;
exo_system.initialize();
end

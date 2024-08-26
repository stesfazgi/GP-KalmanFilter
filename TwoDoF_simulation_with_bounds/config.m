clear all;
close all;
pack all;
clc;

addpath('aux_scripts')

%% Parameters for the system
exo_system = robot_system();
exo_system.dofs = 2;
exo_system.l = [0.3, 0.3];
exo_system.m = [0.45, 0.45];
exo_system.r = [0.03, 0.03];
exo_system.k = [10, 10];
exo_system.d = [0.12, 0.12];
exo_system.j = [0.001, 0.001];
exo_system.Nm = [1, 1];


exo_system.lh = [0.3, 0.4];
exo_system.mh = [2.4, 2.0];
exo_system.rh = [0.06, 0.04];


exo_system.distorder = 4;
S = [zeros(exo_system.distorder-1,1) eye(exo_system.distorder-1);
    0 zeros(1,exo_system.distorder-1)];
H = [1, zeros(1,exo_system.distorder-1)];

exo_system.S = S;
exo_system.H = H;
exo_system.initialize();

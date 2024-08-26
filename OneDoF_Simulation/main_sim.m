% Markus Kessler, Tim Gerstewitz
% 22.02.2023
% % 
%% Load system config
exo_system = config;

%% Initialize Augmented Kalman Filter (AKF)
AKF = GP_AEKF();
% set nominal model for AKF
AKF.nom_model = exo_system;
AKF.sample_time = 0.001;
AKF.enable_learning = 0;
% initialize LoG-GP model in AKF
AKF.init_log_gp()

%% Augmented Kalman Filter parameters
% uncertainty (variance) of assumed disturbance dynamics
Q_KF = zeros(exo_system.dofs*4);
Q_KF = 1.0e-9*eye(exo_system.dofs*4);
Qdist = 1*1.0e-1;
% Append disturbance state covariances
Q_KF = blkdiag(Q_KF, Qdist * eye(exo_system.distorder*exo_system.dofs));
R_KF = 1.0e-9*eye(exo_system.dofs*4);
AKF.Q = Q_KF;
AKF.R = R_KF;

% set simulation phases
sim_phases = [1,2]; % 1: training % 2: estimation with learning (no exo torque) 

%% Run simulation loop
for phase = 1:length(sim_phases)
  % AKF must be returned again to retain trained GP-model
  AKF = run_sim(AKF, sim_phases(phase));
end

function [AKF] = run_sim(AKF,phase)
%% run_sim: Main simulation function playing one phase

%% Parameters for the system;
sys_params.l1 = AKF.nom_model.l(1);     
sys_params.m1 = AKF.nom_model.m(1);  
sys_params.r1 = AKF.nom_model.r(1);
sys_params.k1 = AKF.nom_model.k(1);   
sys_params.d1 = AKF.nom_model.d(1);
sys_params.j = AKF.nom_model.j(1);

% load side inertia
M = 1.32e-2;
g = 9.81; 
% simulated human added to true simulated system
% as human damping
sys_params.dh = 0.12;
% and added load side mass
mh = 0.1;
% and added load side inertia
Mh = 0.8e-2;

% measurement noise covariance
loadsensvar = 1.0e-9;
motorsensvar = 1.0e-9;

% set flags for external torque estimation
use_perturbated_system_model = true;
use_meas_noise = false;
% set system/control flags
use_dynamic_compensator = false;
animate_system = false;
plot_results = true;

sample_time = 0.001;
stop_time = 5;
t_span = (0:sample_time:stop_time)';
learning_start = 3; %2.5/sample_time;
learning_stop = 5/sample_time-1;
% samples skipped during learning:
learning_step = 10;

x_ref.t = t_span;

constraint_max = -20;     constraint_max = deg2rad(constraint_max);
constraint_min = -70;     constraint_min = deg2rad(constraint_min);
beta_deg = 0;             beta_rad = deg2rad(beta_deg); % Safety margin in degrees

b = pi/15;
positional_offset = -45;    positional_offset = deg2rad(positional_offset);
amplitude = 30;             amplitude = deg2rad(amplitude);

x_ref.q1 = pi/2*ones(length(t_span),1);
x_ref.q1_first_deriv = 0*ones(length(t_span),1);
x_ref.q1_second_deriv = 0*ones(length(t_span),1);
x_ref.q1_third_deriv = 0*ones(length(t_span),1);
x_ref.q1_fourth_deriv = 0*ones(length(t_span),1);

% perturbated system parameters
if use_perturbated_system_model
   mp = sys_params.m1+mh;
   Mp = 1.32e-2+Mh;
%    Mp = M;
   dh = sys_params.dh;
else
   mp = sys_params.m1;
   Mp = M;
%    dh = 0.12;
   dh = 0;
   discmomsys.dh = dh;
end

 %% simulated disturbance torque
tau_d = zeros(size(t_span));
sigma_d_sq = 0.1;
tau_d = real(0.5/sqrt(sigma_d_sq*2*pi)*exp(-(t_span-1.5).^2./(2*sigma_d_sq)))-real(0.35/sqrt(sigma_d_sq*2*pi)*exp(-(t_span-3.5).^2./(2*sigma_d_sq)));


 %% System with damping
%x0 = [x_ref.q1(1), x_ref.q1_first_deriv(1), x_ref.q1(1), 0];
x0 = [0, 0, 0, 0];


var_names = {'t', 'q1', 'q1_dot', 'theta1', 'theta1_dot', 'v1', 'u_dot1', 'u1', 'tau1', 'tau1_enforced'};
var_types = {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
system_states = table('Size', [length(t_span),length(var_names)], 'VariableNames', var_names, 'VariableTypes', var_types);

inequality.A = zeros(length(t_span),2);
inequality.b = zeros(length(t_span),2);
inequality.f4x = zeros(length(t_span),1);
inequality.exitflag = zeros(length(t_span),1);

derivatives = zeros(length(t_span),8);

% real state initial covariances = 0 (since initial values known)
P_kfhat_t0 = zeros(4);
% append disturbance state initial covariance = 1 (since initial values unknown)
P_kfhat_t0 = blkdiag(P_kfhat_t0, eye(AKF.nom_model.distorder));

 %% Initialize plot arrays
 
% augmented KF
tau_d_KF_hat_plot = zeros(size(t_span));
theta_KF_plot = zeros(size(t_span));
xs_KF_plot = zeros(size(t_span));
thetad_KF_plot = zeros(size(t_span));
xsd_KF_plot = zeros(size(t_span));

KF_ubs = zeros(size(t_span,1),5);
KF_lbs = zeros(size(t_span,1),5);
KF_cov_plot = zeros(size(t_span,1),5);
KF_K_plot = zeros(size(t_span,1),1);

trans_traindata = zeros(size(t_span,1),AKF.nom_model.dofs);
x_traindata = zeros(size(t_span,1),AKF.nom_model.dofs*3);
u_traindata = zeros(size(t_span,1),1);

restrueplot = zeros(size(t_span,1), 1);

 %% Simulation loop
for i=1:(length(t_span)-1)
    % Start conditions
    if i==1
        current_x = x0;
        tau = 0;
        tloggp = 0;

        K_KF = 0;

        % Augmented KF initial conditions
        P_KF = P_kfhat_t0;
        x_KF_hat_0 = [current_x(3); current_x(1)-current_x(3);
        current_x(4); current_x(2)-current_x(4); zeros(AKF.nom_model.distorder,1)];
        x_KF_hat = x_KF_hat_0;
        AKF.initialize(x_KF_hat,P_KF)
        
        % initialize first plot value
        theta_KF_plot(i) = x_KF_hat(1);
        xs_KF_plot(i) = x_KF_hat(2);
        thetad_KF_plot(i) = x_KF_hat(3);
        xsd_KF_plot(i) = x_KF_hat(4);
        tau_d_KF_hat_plot(i) = x_KF_hat(5);
    else
        current_x = x(end,:);
    end
  
    %% Calculate derivatives
    % Numerical derivatives
    if i>1
        u_traindata(i) = tau;
        
        q_first_deriv_num = (current_x(1)-system_states.q1(i-1))/sample_time;
        q_second_deriv = (q_first_deriv_num-derivatives(i-1,2))/sample_time;
        q_third_deriv = (q_second_deriv-derivatives(i-1,3))/sample_time;

        % measurement from real system
        if use_meas_noise
            % load side angle + vel + acc measurements
            nq_meas_old = nq_meas;
            nq_meas = current_x(1) + sqrt(loadsensvar)*randn(1);
            nqd_meas_old = nqd_meas;
            nqd_meas = (nq_meas - nq_meas_old)/sample_time;
            nqdd_meas = (nqd_meas - nqd_meas_old)/sample_time;
            % motor side angle + vel + acc measurements
            ntheta_meas_old = ntheta_meas;
            ntheta_meas = current_x(3) + sqrt(motorsensvar)*randn(1);
            nthetad_meas_old = nthetad_meas;
            nthetad_meas = (ntheta_meas - ntheta_meas_old)/sample_time;
            nthetadd_meas = (nthetad_meas - nthetad_meas_old)/sample_time;
% 
%               % Apply moving average filter to noisy numerical derivatives
%             window_size = 5;
%             if i>=window_size
%                 q_first_deriv_filtered = mean([derivatives(i-(window_size-1):i-1,5);nqd_meas]);
%                 theta_first_deriv_filtered = mean([derivatives(i-(window_size-1):i-1,7);nthetad_meas]);
%                 qd_f_plot(i) = q_first_deriv_filtered;
%                 thetad_f_plot(i) = theta_first_deriv_filtered;
%             end
% 
%             window_size = 20;
%             if i>=window_size
%                 q_second_deriv_filtered = mean([derivatives(i-(window_size-1):i-1,6);nqdd_meas]);
%                 theta_second_deriv_filtered = mean([derivatives(i-(window_size-1):i-1,8);nthetadd_meas]);
%                 qdd_f_plot(i) = q_second_deriv_filtered;
%                 thetadd_f_plot(i) = theta_second_deriv_filtered;
%             end
            
            % set measured values used by KFs
            theta_meas = ntheta_meas;
            thetad_meas = nthetad_meas;
            thetadd_meas = nthetadd_meas;
            q_meas = nq_meas;
            qd_meas = nqd_meas;
            qdd_meas = nqdd_meas;
  
        else
            % load side angle + vel + acc measurements
            nnqd_meas_old = nnqd_meas;
            nnqd_meas = current_x(2);
            % numerical acceleration
            nnqdd_meas = (nnqd_meas - nnqd_meas_old)/sample_time;
            nnq_meas = current_x(1);
            
            % motor side angle + vel + acc measurements
            nntheta_meas = current_x(3);
            nnthetad_meas_old = nnthetad_meas;
            nnthetad_meas = current_x(4);
            % numerical acceleration
            nnthetadd_meas = (nnthetad_meas - nnthetad_meas_old)/sample_time;
            
            % set measured values used by KFs
            theta_meas = nntheta_meas;
            thetad_meas = nnthetad_meas;
            thetadd_meas = nnthetadd_meas;
            q_meas = nnq_meas;
            qd_meas = nnqd_meas;
            qdd_meas = nnqdd_meas;
        end
        
    else        
        q_first_deriv_num = x0(2);
        q_second_deriv = 0;
        q_third_deriv = 0;

        nqdd_meas = 0;
        nthetadd_meas = 0;
        
        nq_meas = x0(1);
        nqd_meas = x0(2);
        ntheta_meas = x0(3);
        nthetad_meas = x0(4);

        nnq_meas = x0(1);
        nnqd_meas = x0(2);
        nntheta_meas = x0(3);
        nnthetad_meas = x0(4);
    end
    
    % Filter numerical derivatives with moving average that takes the last window_size-1 values plus the present value
%     window_size = 10;
%     if i>=window_size
%         tau_first_deriv = mean([derivatives(i-(window_size-1):i,5);tau_first_deriv]);
%         
%         q_first_deriv_num = mean([derivatives(i-(window_size-1):i,2);q_first_deriv_num]);
%         q_second_deriv = mean([derivatives(i-(window_size-1):i,3);q_second_deriv]);
%         q_third_deriv = mean([derivatives(i-(window_size-1):i,4);q_third_deriv]);
%     else
%         tau_first_deriv = mean([derivatives(1:i-1,5);tau_first_deriv]);
%         
%         q_first_deriv_num = mean([derivatives(1:i-1,2);q_first_deriv_num]);
%         q_second_der2iv = mean([derivatives(1:i-1,3);q_second_deriv]);
%         q_third_deriv = mean([derivatives(1:i-1,4);q_third_deriv]);
%     end
    
    q_first_deriv = q_first_deriv_num;
%     q_first_deriv = current_x(2);

    % Save derivatives for plotting
    derivatives(i,:) = [current_x(1), current_x(2), q_second_deriv, q_third_deriv, nqd_meas, nqdd_meas, nthetad_meas, nthetadd_meas];

    %% Calculate nominal control
    L = [1e10, 1e9, 6.7e7, 250]; % undamped control
    L = [1e7, 1e6, 6.7e4, 250]; % undamped control

    index = find(x_ref.t==t_span(i));
    fprintf("%i: e_0 = %.4f | e_1 = %.4f | e_2 = %.4f | e_3 = %.4f\n", ...
        index, x_ref.q1(index)-current_x(1), x_ref.q1_first_deriv(index)-q_first_deriv, ...
        x_ref.q1_second_deriv(index)-q_second_deriv, x_ref.q1_third_deriv(index)-q_third_deriv);
    
    v = x_ref.q1_fourth_deriv(index) + L*[x_ref.q1(index) - current_x(1);
                                          x_ref.q1_first_deriv(index) - q_first_deriv;
                                          x_ref.q1_second_deriv(index) - q_second_deriv;
                                          x_ref.q1_third_deriv(index) - q_third_deriv];
    
    gamma = 2000;
    beta = beta_rad;
    
    B_max = gamma^3*(constraint_max - current_x(1) - beta) - 3*gamma^2*q_first_deriv - 3*gamma*q_second_deriv - q_third_deriv;
    grad_B_max = [-gamma^3; -3*gamma^2; -3*gamma; -1];

    B_min = gamma^3*(-constraint_min + current_x(1) - beta) + 3*gamma^2*q_first_deriv + 3*gamma*q_second_deriv + q_third_deriv; 
    grad_B_min = [gamma^3; 3*gamma^2; 3*gamma; 1];
                                      
        %% Control for undamped model
        f4x = sys_params.k1/(sys_params.j*Mp) * (-Mp*q_second_deriv + mp*g*sys_params.l1/2*cos(current_x(1)) - dh*current_x(2))...
            - 1/Mp*mp*g*sys_params.l1/2*(cos(current_x(1))*q_first_deriv^2+sin(current_x(1))*q_second_deriv) - sys_params.k1/Mp*q_second_deriv;
        g4x = sys_params.k1/(sys_params.j*Mp);
 
        tau = 1/g4x * (v - f4x);
%         tau_enforced = 1/g4x * (v_enforced - f4x);
        
        f_x = [q_first_deriv; q_second_deriv; q_third_deriv; f4x];
        g_x = [0; 0; 0; g4x];
        
        A = -[dot(grad_B_max,g_x); dot(grad_B_min,g_x)];
        b = [dot(grad_B_max,f_x)+B_max; dot(grad_B_min,f_x)+B_min];

        inequality.A(i,:) = A';
        inequality.b(i,:) = b';
        inequality.f4x(i) = f4x;
        
        opts_quadprog = optimset('Display','none');
        [tau_enforced,~,exitflag] = quadprog(eye(1),-tau,A,b,[],[],[],[],[],opts_quadprog);
        inequality.exitflag(i) = exitflag;

        u_dot = NaN;
        u = NaN;
        
    %% Save system states and compute next time step
    % Save system states
    system_states = copySystemStatesToTable(system_states,i,t_span(i),current_x,v,u_dot,u,tau,tau_enforced);
    tau=0;   
    
    %% Estimate external torque (GP-AEKF)
    if i>1
        %% Parametric model prediction & LoG-GP training data
        if phase == 1
          if i == learning_start
            disp('start training')
          end
          % set motor torque to excite system similarily to the disturbance
          tau = tau_d(i);
          % learn without disturbance
          tau_d(i) = 0;
          
          if (i >= learning_start) && (i < learning_stop)
            % check if current point should be added to training set
            % (determined by learning_step)
            if mod(i-learning_start ,learning_step) == 0
              % parametric model prediction from previous state
              Ghat = sys_params.m1*g*sys_params.l1/2*cos(q_meas);
              uhat = M*qdd_meas - Ghat + sys_params.d1*qd_meas + sys_params.j*thetadd_meas;
              
              fh = tau-uhat;

              % current state x
              xnew = [current_x(1) current_x(2) qdd_meas];
              x_traindata(i,:)  = xnew;
              % set training target
              trans_traindata(i,:) = fh;
                
              Gh = mh*g*sys_params.l1/2*cos(q_meas);
              truetrans = Mh*qdd_meas + dh*qd_meas - Gh; 
              restrueplot(i) = truetrans;
              
              tic
              % update GP for each output dimension
                AKF.train_log_gp(tau,xnew,thetadd_meas)
              tloggp = tloggp + toc;
            end
          end
        else
        %% Kalman Filter augmented state estimation approach
        % spring deformation measurement
        x_s_meas = q_meas-theta_meas;
        % augmented KF measurement vector
        y_meas = [theta_meas, x_s_meas, thetad_meas, qd_meas-thetad_meas]';

        AKF.predict(tau);
        AKF.update(y_meas);
        % Get state and covariance estimate from augmented state KF
        x_KF_hat = AKF.x_hat;
        P_KF = AKF.P;

        % save estimates for plot
        theta_KF_plot(i) = x_KF_hat(1);
        xs_KF_plot(i) = x_KF_hat(2);
        thetad_KF_plot(i) = x_KF_hat(3);
        xsd_KF_plot(i) = x_KF_hat(4);
        tau_d_KF_hat_plot(i) = x_KF_hat(5);
        P_KF_bounds = 2*sqrt(P_KF);
        for k = 1:5
          KF_ubs(i,k) = x_KF_hat(k) + P_KF_bounds(k,k);
          KF_lbs(i,k) = x_KF_hat(k) - P_KF_bounds(k,k);
          KF_cov_plot(i,k) = P_KF_bounds(k,k);
        end
        
        end
    
    end
    
    %% Compute system at next time step
    tau = 0;
    xold = [current_x(3) current_x(1)-current_x(3) current_x(4) current_x(2)-current_x(4)];
    [t,x] = ode45(@(t,x) oneDofPlanarRobot(t,x,tau,tau_d(i),sys_params,use_perturbated_system_model, mp, Mp, dh), [t_span(i),t_span(i+1)], current_x);
end

system_states = system_states(1:end-1,:);
if phase == 1
    tloggp = tloggp/(learning_stop-learning_start);
else
    tloggp = tloggp/(stop_time/sample_time);
end
%% Plots and Animation
if animate_system
    [fig_simulation, video_handle] = drawAnimatedSystem(system_states.t,system_states{:,2:5},x_ref,sys_params);
end

%% Plots
length_t = length(system_states.t);
plot_step = 1;

% Colors for plotting
% Define primary TUM colors
tum_blue = [0 1 189]./255;
% Define secondary TUM colors
tum_dark_gray = [51, 51, 51]./255; tum_gray = [128, 128, 128]./255; tum_light_gray = [204, 204, 204]./255;
% Define accent TUM colors
tum_green = [162, 173, 0]./255; tum_orange = [227, 114, 34]./255; tum_light_blue = [152, 198, 234]./255;

if plot_results
    % Plot position and torque
    % fig = figure('Name', 'Position and Torque')
%     figure(1*phase);
%     plot(system_states.t(1:plot_step:length_t),system_states.q1(1:plot_step:length_t),'Color',tum_blue,'LineWidth',1.5); hold on;
%     plot(system_states.t(1:plot_step:length_t),x_ref.q1(1:plot_step:end-1),'Color',tum_orange,'LineStyle','--');
%     plot(system_states.t(1:plot_step:length_t),constraint_max*ones(length_t/plot_step),'k','LineStyle','-.','LineWidth',1);
%     plot(system_states.t(1:plot_step:length_t),constraint_min*ones(length_t/plot_step),'k','LineStyle','-.','LineWidth',1);
%     legend(["Simulated trajectory"; "Reference trajectory"; "Constraints"], 'Location', 'Northwest');
%     title('Link positions');
%     grid on;
%     set(findall(gcf,'-property','FontSize'),'FontSize',16)

    %% Plot disturbance torque estimates
    figure(2*phase);
    plot(system_states.t(1:plot_step:length_t),tau_d_KF_hat_plot(1:plot_step:length_t),'b','LineWidth',2.5);hold on;
    plot(system_states.t(1:plot_step:length_t),KF_ubs(1:plot_step:length_t,5),'-.','Color', [0, 1, 1, 1],'LineWidth',1.5); hold on;
    plot(system_states.t(1:plot_step:length_t),KF_lbs(1:plot_step:length_t,5),'-.','Color', [0, 1, 1, 1],'LineWidth',1.5); hold on;

    plot(system_states.t(1:plot_step:length_t),tau_d(1:plot_step:length_t),'--','Color',tum_orange,'LineWidth',3); hold on;
    legend('estimated external torque augmented state KF','2\sigma ucb augmented state KF','2\sigma lcb augmented state KF','disturbance torque');

    title('External torque estimation');
    grid on;
    set(findall(gcf,'-property','FontSize'),'FontSize',20)

    %% Plot augmented state KF state estimates
    figure(3*phase);
    title('State estimates');

    subplot(2,2,1)
    plot(system_states.t(1:plot_step:length_t),KF_ubs(1:plot_step:length_t,1),'c-.','LineWidth',1.2); hold on;
    plot(system_states.t(1:plot_step:length_t),KF_lbs(1:plot_step:length_t,1),'c-.','LineWidth',1.2); hold on;
    plot(system_states.t(1:plot_step:length_t),system_states.theta1(1:plot_step:length_t),'r','LineWidth',1.5); hold on;
    plot(system_states.t(1:plot_step:length_t),theta_KF_plot(1:plot_step:length_t),'b--','LineWidth',1.5); hold on;
    legend('2\sigma ucb','2\sigma lcb','$\theta_m$ true','\theta_m estimate');
    grid on;
    subplot(2,2,2)
    plot(system_states.t(1:plot_step:length_t),KF_ubs(1:plot_step:length_t,3),'c-.','LineWidth',1.2); hold on;
    plot(system_states.t(1:plot_step:length_t),KF_lbs(1:plot_step:length_t,3),'c-.','LineWidth',1.2); hold on;
    plot(system_states.t(1:plot_step:length_t),system_states.theta1_dot(1:plot_step:length_t),'r','LineWidth',1.5); hold on;
    plot(system_states.t(1:plot_step:length_t),thetad_KF_plot(1:plot_step:length_t),'b--','LineWidth',1.5); hold on;
    legend('2\sigma ucb','2\sigma lcb','\dot{\theta}_m true','\dot{\theta}_m estimate');
    grid on;
    subplot(2,2,3)
    plot(system_states.t(1:plot_step:length_t),KF_ubs(1:plot_step:length_t,2),'c-.','LineWidth',1.2); hold on;
    plot(system_states.t(1:plot_step:length_t),KF_lbs(1:plot_step:length_t,2),'c-.','LineWidth',1.2); hold on;
    plot(system_states.t(1:plot_step:length_t),system_states.q1(1:plot_step:length_t) - system_states.theta1(1:plot_step:length_t),'r','LineWidth',1.5); hold on;
    plot(system_states.t(1:plot_step:length_t),xs_KF_plot(1:plot_step:length_t),'b--','LineWidth',1.5); hold on;
    legend('2\sigma ucb','2\sigma lcb','x_s true','x_s estimate');
    grid on;
    subplot(2,2,4)
    plot(system_states.t(1:plot_step:length_t),KF_ubs(1:plot_step:length_t,4),'c-.','LineWidth',1.2); hold on;
    plot(system_states.t(1:plot_step:length_t),KF_lbs(1:plot_step:length_t,4),'c-.','LineWidth',1.2); hold on;
    plot(system_states.t(1:plot_step:length_t),system_states.q1_dot(1:plot_step:length_t) - system_states.theta1_dot(1:plot_step:length_t),'r','LineWidth',1.5); hold on;
    plot(system_states.t(1:plot_step:length_t),xsd_KF_plot(1:plot_step:length_t),'b--','LineWidth',1.5); hold on;
    legend('2\sigma ucb','2\sigma lcb','\dot{x}_s true','\dot{x}_s estimate');
    grid on;
    set(findall(gcf,'-property','FontSize'),'FontSize',16)

    %% Plot augmented state KF state estimate errors
    % figure(4*phase);
    % title('State estimate errors');
    %
    % subplot(2,2,1)
    % plot(system_states.t(1:plot_step:length_t),theta_KF_plot(1:plot_step:length_t)-system_states.theta1(1:plot_step:length_t),'r','LineWidth',1.5); hold on;
    % legend('\theta_m error');
    % grid on;
    % subplot(2,2,2)
    % plot(system_states.t(1:plot_step:length_t),thetad_KF_plot(1:plot_step:length_t)-system_states.theta1_dot(1:plot_step:length_t),'r','LineWidth',1.5); hold on;
    % legend('\dot{\theta}_m error');
    % grid on;
    % subplot(2,2,3)
    % plot(system_states.t(1:plot_step:length_t),xs_KF_plot(1:plot_step:length_t) - (system_states.q1(1:plot_step:length_t) - system_states.theta1(1:plot_step:length_t)),'r','LineWidth',1.5); hold on;
    % legend('x_s error');
    % grid on;
    % subplot(2,2,4)
    % plot(system_states.t(1:plot_step:length_t),xsd_KF_plot(1:plot_step:length_t) - (system_states.q1_dot(1:plot_step:length_t) - system_states.theta1_dot(1:plot_step:length_t)),'r','LineWidth',1.5); hold on;
    % legend('\dot{x}_s error');
    % grid on;
    % set(findall(gcf,'-property','FontSize'),'FontSize',16)

    %% Plot GP-AKF estimate variance
%     figure(5*phase);
%     plot(system_states.t(1:plot_step:length_t),KF_cov_plot(1:plot_step:length_t,:),'LineWidth',2.5); hold on;
%     title('Augmented KF covariance estimate');
%     legend;
%     grid on;
%     set(findall(gcf,'-property','FontSize'),'FontSize',20)

end

if phase == 1
    x_traindata = x_traindata(learning_start:learning_step:learning_stop,:);
    trans_traindata = trans_traindata(learning_start:learning_step:learning_stop,:);
    restrueplot = restrueplot(learning_start:learning_step:learning_stop,:);
    save loggptraindat x_traindata trans_traindata restrueplot
end

end
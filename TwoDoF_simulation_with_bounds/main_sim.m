% Markus Kessler, 22.02.2023

%% Load system config
config

%% Initialize Augmented Kalman Filter (AKF)
AKF = GP_AEKF_bounds();
% set nominal model for AKF
AKF.nom_model = exo_system;
AKF.sample_time = 0.01;
% initialize LoG-GP model in AKF

% set training data region X for LoG-GP error bounds based on the given
% trajectory
b12 = pi/15;
Xte = [-pi/2+65*pi/180,-pi/2;
      3/4*pi,0;
      65*pi/180*b12,-65*pi/180*b12
      3/4*pi*b12,-3/4*pi*b12;
      0.0883,-0.2721;
      0.5596,-0.2055];
      
% initialize LoG-GP models
AKF.init_log_gp(Xte)

%% Augmented Kalman Filter parameters2
% uncertainty (variance) of assumed disturbance dynamics
Q_KF = 1*1.0e-3*eye(exo_system.dofs*4);
Qdist = 1*1.0e-2;
% Append disturbance state covariances
Q_KF = blkdiag(Q_KF, Qdist * eye(exo_system.distorder*2));
R_KF = 1.0e-9*eye(exo_system.dofs*4);
AKF.Q = Q_KF;
AKF.R = R_KF;

% set simulation phases
sim_phases = [1,2]; % 1: training % 2: estimation with learning (no exo torque) 
% 3: estimation with learning (exo active) % 4: with GP uncertainty (exo
% counteracting)

%% Run simulation loop
for phase = 1:length(sim_phases)
  % AKF must be returned again to retain trained GP-model
  AKF = run_sim(AKF, sim_phases(phase));
end

function AKF = run_sim(AKF, sim_phase)
%% run_sim: Main simulation function playing one phase

% plot results
draw_system_states      = true;
% animate simulation
animate_system          = false;
% save results in .mat file
save_results            = false;
save_training_data      = true;

% simulation sample time and duration
sample_time = AKF.sample_time;
stop_time = 15;
t_span = (0:sample_time:stop_time)';

%% LoG-GP torque learning phase
% step where training begins
learning_start = 3;
% step where training stops
learning_stop = stop_time/sample_time-1;
% samples skipped during learning:
learning_step = 1;      % 1 = all samples used, 10 = every 10th sample used
remove_data_start = 0;
remove_data_end = 0;

%% Construct reference trajectory
x_ref = constructReferenceStates([],[],t_span);

% simulated artificial motor torque (phase 3)
tau_h = zeros(size(t_span,1), 2);
sigma_d_sq = 1;
tau_exo1 = -real(3/sqrt(sigma_d_sq*2*pi)*exp(-(t_span-3.5).^2./(2*sigma_d_sq)))-real(2/sqrt(sigma_d_sq*2*pi)*exp(-(t_span-11.5).^2./(2*sigma_d_sq)));
sigma_d_sq = 0.8;
tau_exo2 = -real(2/sqrt(sigma_d_sq*2*pi)*exp(-(t_span-2.5).^2./(2*sigma_d_sq)))-real(1/sqrt(sigma_d_sq*2*pi)*exp(-(t_span-12.5).^2./(2*sigma_d_sq)));

%% Initialize system for simulation
x0 = [x_ref.q1(1), x_ref.q1_first_deriv(1), 0, 0, ...
      x_ref.q2(1), x_ref.q2_first_deriv(1), 0, 0];
var_names = {'t', 'q1', 'q1_dot', 'theta1', 'theta1_dot', 'q2', 'q2_dot', 'theta2', 'theta2_dot', ...
             'v1', 'v2', 'u_dot1', 'u_dot2', 'u1', 'u2', 'tau1', 'tau2'};
var_types = {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', ...
            'double', 'double', 'double', 'double', 'double', 'double' 'double', 'double'}; 
system_states = table('Size', [length(t_span),length(var_names)], 'VariableNames', var_names, 'VariableTypes', var_types);

% Start conditions
current_x = x0;
tau_m = zeros(2,1);
% Save system states
system_states = copySystemStatesToTable(system_states,1,t_span(1),current_x,tau_m);

%% Augmented KF initializationexo_sys
% initial AKF state
% initial system state known (uncertainty = 0), external torque
% unknown (uncertainty = 1)
P_KF = blkdiag(0.01*eye(AKF.nom_model.dofs*4),0.05*eye(AKF.nom_model.distorder*AKF.nom_model.dofs));

x_KF_hat = [current_x(3); current_x(7); current_x(1)-current_x(3)/AKF.nom_model.Nm(1); current_x(5)-current_x(7)/AKF.nom_model.Nm(2);...
    current_x(4); current_x(8); current_x(2)-current_x(4)/AKF.nom_model.Nm(1); current_x(6)-current_x(8)/AKF.nom_model.Nm(2);...
    0; zeros(AKF.nom_model.distorder-1,1);0.3800;zeros(AKF.nom_model.distorder-1,1)];
% initial error set (zero)
X_0 = zeros(AKF.nom_model.dofs*(4+AKF.nom_model.distorder));
% set AKF initial values
AKF.initialize(x_KF_hat,P_KF,X_0)

%% setup plotting storage variables
plots = init_plot_data(t_span);
% initialize first plot value
plots.x_KF_plot(1,:) = x_KF_hat(1:8);
plots.tau_d_KF_hat_plot(1,:) = [x_KF_hat(9); x_KF_hat(9+AKF.nom_model.distorder)];
% initialize computing time storage
tloggp = 0;
t_akf_pred = 0;
t_akf_upd = 0;

%% start simulation loop
for i=1:(length(t_span)-1)
    %% Compute control for system
    [theta_second_deriv,M,n] = ...
            StateVariablesHigherDerivatives(current_x,tau_m,AKF.nom_model);
    tau_m = ctc_controller(AKF.nom_model, current_x, theta_second_deriv, M, n, x_ref, i, sim_phase);
    if sim_phase >= 2
        % human becomes active
        tau_h(i,:) = ctc_controller_human(AKF.nom_model, current_x, theta_second_deriv, M, n, x_ref, i, sim_phase);
        if sim_phase >= 3
            % active exo (additionally)
            tau_m = AKF.nom_model.Nminv*[tau_exo1(i); tau_exo2(i)];
        else
            % passive exo
            tau_m = [0; 0];
        end
    end

    %% Save system states
    system_states = copySystemStatesToTable(system_states,i,t_span(i),current_x,tau_m);
    
    % wait for second step (numerical differentiation)
    if i>1
        %% Get fresh measurements from system
        meas = get_measurements(system_states,current_x,i,sample_time);

        if sim_phase == 1
            %% LoG-GP torque residual learning
            if i == learning_start
                disp('start training')
            end
            if (i >= learning_start) && (i < learning_stop)
                % check if current point should be added to training set
                % (determined by learning_step):
                if mod(i-learning_start ,learning_step) == 0
                    if i <= remove_data_start || i >= remove_data_end
                    % current state x
                    x_train = [meas.q' meas.qd' meas.qdd'];
                    tic
                    % update LoG-GP with new sample
                    y_train = AKF.train_log_gp(tau_m,x_train,meas.thetadd_m);

                    tupdate = toc;
                    tloggp = tloggp + tupdate;
                    plots.gptime_plot(i:i+learning_step) = tupdate/2;
                    plots.q(i,:) = meas.q;
                    plots.qd(i,:) = meas.qd;
                    plots.qdd(i,:) = meas.qdd;
                    
                    plots.xtrain(i,:) = x_train;
                    plots.ytrain(i,:) = y_train;
                    end
                end
            end

        else
            %% Kalman Filter augmented state estimation approach
            tic

            %% AKF prediction (x_k|k-1)
            jacvec = AKF.predict(tau_m);

            tpredict = toc;
            t_akf_pred = t_akf_pred+tpredict;

            tic
            %% AKF update step (x_k|k-1)->(x_k|k) using measurement
            % augmented KF measurement vector
            y_meas = [meas.theta_m; meas.q-AKF.nom_model.Nminv*meas.theta_m;...
                meas.thetad_m; meas.qd-AKF.nom_model.Nminv*meas.thetad_m];
            
            % update estimate using fresh measurement (y_meas)
            AKF.update(y_meas);

            tupdate = toc;
            t_akf_upd = t_akf_upd+tupdate;

            % Get state and covariance estimate from augmented state KF
            x_KF_hat = AKF.x_hat;
            P_KF = AKF.P;

            % save plot values
            plots.x_KF_plot(i,:) = x_KF_hat(1:8);
            plots.tau_d_KF_hat_plot(i,1) = x_KF_hat(9);
            plots.tau_d_KF_hat_plot(i,2) = x_KF_hat(9+AKF.nom_model.distorder);
            
            plots.q(i,:) = meas.q;
            plots.qd(i,:) = meas.qd;
            plots.qdd(i,:) = meas.qdd;
            plots.jacvec(i,:) = jacvec;
                    
            % calculate confidence bounds for estimate plot
            % uncertainty sigma bounds
            % MISSING CORRECT s FACTOR FROM CHI-SQUARE DISTRIBUTION!
            P_KF_bounds = sqrt(AKF.delta.*P_KF);
            % interval bounds from ellipsoid shape matrix
            X_mean_bounds = sqrt(diag(AKF.X_hat_est));
            % ids of estimated SEA states + active human torque
            ids = [1:8, 9, 9+AKF.nom_model.distorder];
            id = 0;
            % calcilate upper and lower bounds for each state
            for k = ids
                id = id+1;
                plots.AKF_ubs(i,id) = x_KF_hat(k) + P_KF_bounds(k,k) + X_mean_bounds(k);
                plots.AKF_lbs(i,id) = x_KF_hat(k) - P_KF_bounds(k,k) - X_mean_bounds(k);
                plots.AKF_cov_plot(i,id) = P_KF_bounds(k,k);
            end
        end
    end
    
    %% Compute system at next time step i+1
    [t,x] = ode45(@(t,x) twoDofPlanarRobotWithDamping(t,x,tau_m, tau_h(i,:), AKF.nom_model), [t_span(i),t_span(i+1)], current_x);
    fprintf("%i | %.4f\n", i, t(end));
    current_x = x(end,:);
end

%% Get final system states and computation times
system_states = system_states(1:stop_time/sample_time,:);
x_ref = x_ref(1:stop_time/sample_time,:);

if sim_phase == 1
    tloggp = tloggp/(learning_stop-learning_start);
    disp(['average LoG-GP update step time: ', num2str(tloggp), 's'])
else
    t_akf_pred = t_akf_pred/(length(t_span)-1);
    disp(['average AKF prediction step time: ', num2str(t_akf_pred), 's'])
    t_akf_upd = t_akf_upd/(length(t_span)-1);
    disp(['average AKF update step time: ', num2str(t_akf_upd), 's'])
end

%% Plots and Animation
if draw_system_states
    fig_system_states = drawSystemStatesNew(system_states,x_ref, tau_h, plots);
end

if animate_system
  tau_m = [system_states.tau1, system_states.tau2];
  l = AKF.nom_model.l;
  [fig_simulation, video_handle] = drawAnimatedSystem(system_states.t,system_states{:,2:9},x_ref,l,tau_m, tau_h,plots);
end
    
% save simulation results in file
if save_results == true
    tau_m = [system_states.tau1, system_states.tau2];
    t = system_states.t;
    x = system_states{:,2:9};
    results_file = ['simdata_phase_' num2str(sim_phase) '.mat'];
    l = AKF.nom_model.l;
    save(results_file, 't','x','x_ref','l','tau_m', 'tau_h','plots');
end

if save_training_data == true && sim_phase == 1
  x_train_plot = plots.xtrain(learning_start:learning_stop-1,:);
  y_train_plot = plots.ytrain(learning_start:learning_stop-1,:);
  save('loggptraindat.mat','x_train_plot','y_train_plot')
end

end
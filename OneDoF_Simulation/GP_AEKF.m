classdef GP_AEKF <handle
% class for external torque estimation using augmented Kalman filter (AKF)
% Estimates and their covariances can be accessed from properties x_hat and P
% Last modified: Markus Kessler 21.02.2023

    properties
        %% Parameters for the system
        enable_learning = 1;    % enable GP based residual learning 
        sample_time;    
        Q;                      % process noise covariance matrix
        R;                      % measurement noise covariance matrix

        % augmented state vector x_hat=[theta_m, theta_s, \dot{theta}_m, \dot{theta}_s,\omega_{ext}]
        x_hat;                  % augmented state estimate
        P;                      % estimate covariance matrix
        
        % nominal system model (specified by exo_system and config)
        nom_model;

        % select discretization mode: 1 = 'Euler' (faster), 2 = 'Matrix_exponential' (more accurate, slower)
        discretization_mode = 2;
    end

    properties(Access = protected)
        x_hat_pred;             % estimate prediction
        P_pred;                 % prediction covariance 
        x_hat_km1;              % previous estimate (differentiation)
        state_eye;
        disc_sys;               % discretized, linearized nominal system
        
        loggp;                  % LoG-GP model          
        mu_gp;                  % GP mean prediction
        var_gp;                 % GP variance
        jac_gp;                 % GP mean gradient
        dofs;                   % degrees of freedom
        distorder;              % order of assumed external torque dynamics
    end

    methods
        function initialize(obj,x0,P0)
            %% initialize: Setup constant variables and initial states
            % inputs:   -x0:    initial augmented state
            %           -P0:    initial covariance
            
            disp('Initializing augmented Kalman filter...');
            obj.dofs = obj.nom_model.dofs;
            obj.distorder = obj.nom_model.distorder;
            % initial state
            obj.x_hat = x0;
            obj.x_hat_km1 = obj.x_hat;
            % initial covariance
            obj.P = P0;

            obj.state_eye = eye(obj.dofs*(4+obj.distorder));
            obj.mu_gp = zeros(obj.dofs,1);
            obj.var_gp = zeros(obj.dofs,1);
            obj.jac_gp = zeros(obj.dofs*4,obj.dofs);
        end
        
        function init_log_gp(obj)
            %% init_log_gp: initialize LoG-GP models for torque residual learning
            
            disp('Setup of MoE-LoG-GP ...');
            % one model for each DOF
            for i = 1:obj.nom_model.dofs
                % use mixture of GP experts LoG-GPs
                MoE_LoG_GP(i) = MOE_LOG_GP();
                
                % set GP hyperparameters
                paramfile = ['hyperparam_dim_' num2str(i) '.mat'];
                load(paramfile, 'ls', 'sf', 'sn')
                MoE_LoG_GP(i).sigL = min(10^6,ls);
                MoE_LoG_GP(i).sigF = sf;
                MoE_LoG_GP(i).sigN = sn;
                
                % LoG-GP parameters (see MOE_LOG_GP.m for explanation)
                MoE_LoG_GP(i).pts = 100;
                MoE_LoG_GP(i).loadHyp = false;
                MoE_LoG_GP(i).wo = 100;
                MoE_LoG_GP(i).divMethod = 1;
                MoE_LoG_GP(i).xSize = 3*obj.nom_model.dofs;

                try
                    MoE_LoG_GP(i).N = 10;
                    MoE_LoG_GP(i).setupData();
                catch
                    MoE_LoG_GP(i).N = 10;
                    MoE_LoG_GP(i).setupData();
                end

            end
            % return LoG-GP models
            obj.loggp = MoE_LoG_GP;
        end

        function f_h = train_log_gp(obj,tau_m, x_train, theta_m_dd_meas)
            %% train_log_gp: Add training sample to LoG-GP models
            % inputs:   -tau_m:             motor torques vector
            %           -x_train:           training input data (q,qd,qdd)
            %           -theta_m_dd_meas:   motor acceleration measurements

            % load side variables
            q = x_train(1:obj.dofs)';
            qd = x_train(obj.dofs+1:obj.dofs*2)';
            qdd = x_train(2*obj.dofs+1:obj.dofs*3)';

            % parametric model transition
            uhat = obj.nom_model.inverse_dynamics(q, qd, qdd);
            % motor side solved for spring torque (substituted into load
            % side afterwards)
            u_motor = diag(obj.nom_model.Nm)*(tau_m-obj.nom_model.J*theta_m_dd_meas);

            % calculate torque residual
            f_h = u_motor-uhat;

            % update GP for each output dimension
            for dim = 1:obj.dofs
                obj.loggp(dim).update(x_train',f_h(dim));
            end
        end

        function predict(obj,tau_m)
            %% predict: Perform EKF prediction step based on nominal model
            % inputs:   -tau_m: motor torques vector
            
            % set nominal model state
            obj.nom_model.x = obj.x_hat;
            % Prediction from nominal model (time-continuous fwd dynamics)
            ft = obj.nom_model.fwd_dynamics(tau_m);
            % get linearized system matrices from nominal system
            [linsys.A, linsys.B, linsys.C] = obj.nom_model.fwd_dyn_jacobian(tau_m);
            % process noise
            linsys.Q = obj.Q;

            if obj.enable_learning == 1
                % Prediction from learned GP model
                [add_gp_trans, add_gp_jac, add_gp_var] = obj.get_GP_prediction();

                % add GP learned forward transition
                ft(obj.dofs*3+1:obj.dofs*4) = ft(obj.dofs*3+1:obj.dofs*4) - add_gp_trans;

                % add GP variance to model uncertainty Q
                linsys.Q(obj.dofs*3+1:obj.dofs*4,obj.dofs*3+1:obj.dofs*4) =...
                    linsys.Q(obj.dofs*3+1:obj.dofs*4,obj.dofs*3+1:obj.dofs*4)+add_gp_var;

                % add GP Jacobian to AKF linearization
                linsys.A(3*obj.dofs+1:4*obj.dofs,1:4*obj.dofs) = linsys.A(3*obj.dofs+1:4*obj.dofs,1:4*obj.dofs)-add_gp_jac;
            end
            % discretize linearized forward dynamics
            obj.disc_sys = obj.get_discrete_linear_sys(linsys);
            % perform estimate prediction step with Euler integration
            obj.x_hat_pred = obj.x_hat + obj.sample_time*ft;
            % calculate prediction uncertainty
            obj.P_pred = obj.disc_sys.A * obj.P * obj.disc_sys.A' + obj.disc_sys.Q;
        end

        function [add_gp_trans, add_gp_jac, add_gp_var] = get_GP_prediction(obj)
            %% get_GP_prediction: Predict torque residual from learned GP model
            % outputs:   -add_gp_trans: Forward dynamics portion calculated
            %                           from learned residual torque times inverted inertia 
            %            -add_gp_jac: Jacobian of add_gp_trans
            %            -add_gp_var: Variance of add_gp_trans (from GP variance) 

            % calculate GP input = (q,qd,qdd)
            q = obj.x_hat(1:obj.dofs)+obj.x_hat(obj.dofs+1:obj.dofs*2);
            qd = obj.x_hat(2*obj.dofs+1:obj.dofs*3)+obj.x_hat(3*obj.dofs+1:obj.dofs*4);
            theta_m_dd_num = (obj.x_hat(2*obj.dofs+1:obj.dofs*3)-obj.x_hat_km1(2*obj.dofs+1:obj.dofs*3))./obj.sample_time;
            theta_s_dd_num = (obj.x_hat(3*obj.dofs+1:obj.dofs*4)-obj.x_hat_km1(3*obj.dofs+1:obj.dofs*4))./obj.sample_time;
            qdd = obj.nom_model.Nminv * theta_m_dd_num+theta_s_dd_num;

            % GP input
            x_star = [q; qd; qdd];
            
            for dim = 1:obj.dofs
                % approximate torque residual from LoG-GP
                [mu,sig,jacvec,~] = obj.loggp(dim).predict(x_star);
                % mean approximation
                obj.mu_gp(dim) = mu;
                % approximation variance
                obj.var_gp(dim) = sig;
                % jacobian vector of GP mean to be added to AKF Jacobian(repeat entries since Jacobian
                % is calculated with respect to q,q_dot in the GP, not
                % theta_m,theta_s,theta_m_dot,theta_s_dot. However,
                % gradients with respect to theta_s are identical to
                % gradients with respect to q. Gradients with respect
                % to theta_m are identical to Nm times the gradient with
                % respect to q due to q=Nm^-1 * theta_m + theta_s.
                jac_theta = repmat(jacvec(1:obj.dofs),obj.dofs,1);
                jac_theta_d = repmat(jacvec(obj.dofs+1:obj.dofs*2),obj.dofs,1);
                obj.jac_gp(:,dim) = [jac_theta(:); jac_theta_d(:)];
            end
            
            % get inertia matrix
            M = obj.nom_model.get_M_numerical();
            % invert inertia matrix (for adding GP to fwd dynamics)
            % Minv IS REQUIRED SINCE THE FORWARD DYNAMICS OF AKF ARE:
            % f_fwd = Minv*(torques_nominal + torques_GP)
            Minv = inv(M);
            % get derivative of inverted inertia matrix (for AKF Jacobian)
            dMinvdq2 = obj.nom_model.get_Minv_deriv();

            % calculate Jacobian for AKF from LoG-GP Jacobian and inverted Inertia matrix with product rule
            add_gp_jac = Minv*obj.jac_gp';
            add_gp_jac(:,2) = add_gp_jac(:,2) + dMinvdq2*obj.mu_gp;
            add_gp_jac(:,4) = add_gp_jac(:,4) + dMinvdq2*obj.mu_gp;
            % get indices of augmented state corresponding to theta_m and
            % theta_m_dot
            theta_m_ids = [1:obj.dofs,2*obj.dofs+1:3*obj.dofs];
            % multiply Jacobian with respect to motor theta_m and
            % theta_m_dot by gear ratio 
            add_gp_jac(:,theta_m_ids) = diag(obj.nom_model.Nm) * add_gp_jac(:,theta_m_ids);

            % Calculate forward dynamics portion to be added to AKF prediction
            add_gp_trans = Minv*obj.mu_gp;

            % GP uncertainty added to Q
            add_gp_var = Minv*diag(obj.var_gp)*Minv';
        end

        function update(obj, y_meas)
            %% update: Update augmented state and covariance estimates with EKF
            % New STATE ESTIMATE will be stored and available in obj.x_hat
            % ESTIMATE COVARIANCE will  be stored and available in obj.P
            % inputs:   -y_meas: measurement vector (motor/spring positions&velocities)
            %                    (theta_m_meas, theta_s_meas,
            %                    thetad_m_meas, thetad_s_meas)

            % save old estimate for derivative
            obj.x_hat_km1 = obj.x_hat;

            % compute optimal Kalman gain
            K_KF = obj.P_pred * obj.disc_sys.C' * inv(obj.disc_sys.C * obj.P_pred * obj.disc_sys.C' + obj.disc_sys.R);

            % A posteriori state estimation update with measurement y_meas
            innovation = y_meas - obj.disc_sys.C * obj.x_hat_pred;
            obj.x_hat = obj.x_hat_pred + K_KF * innovation;

            % A posteriori covariance estimation update
            obj.P = (obj.state_eye - K_KF * obj.disc_sys.C) * obj.P_pred;
        end

        function disc_sys = get_discrete_linear_sys(obj, linsys)
            %% get_discrete_linear_sys: Discretize linearized forward dynamics for AKF
            % inputs:   -linsys:    Struct containing all time-continuous
            %                       system matrices (A,B,C,Q,R)
            % outputs:  -disc_sys:  Struct containing all discretized
            %                       system matrices (A,B,C,Q,R)

            if obj.discretization_mode == 1
                % simple Euler (zero-order hold) discretization
                disc_sys.A = obj.state_eye + linsys.A .*obj.sample_time;
                disc_sys.B = linsys.B .* obj.sample_time;
                disc_sys.Q = linsys.Q .* obj.sample_time;

            elseif obj.discretization_mode == 2
                % more accurate matrix exponential discretization (slower!)
                % prepare matrices for discretization
                H1 = [linsys.A linsys.B;
                    zeros(size(linsys.B,2),size(linsys.A,2)) zeros(size(linsys.B,2),size(linsys.B,2))];
                H2 = [-linsys.A linsys.Q;
                    zeros(size(linsys.A,1),size(linsys.A,2)) linsys.A'];
                % matrix exponential
                E1 = expm(H1*obj.sample_time);
                E2 = expm(H2*obj.sample_time);

                % discretized system matrizes (matrix exponential based)
                disc_sys.A = E1(1:size(linsys.A,1), 1:size(linsys.A,2));
                disc_sys.B = E1(1:size(linsys.B,1), size(linsys.A,2)+1:end);

                % discretized process noise
                M_2DKF = E2(1:size(linsys.A,1), size(linsys.A,2)+1:end);
                M_3DKF = E2(size(linsys.A,1)+1:end, size(linsys.A,2)+1:end);
                disc_sys.Q = M_3DKF.' * M_2DKF;
            end
            % discrete measurement noise
            disc_sys.R = obj.R./obj.sample_time;
            disc_sys.C = linsys.C;
        end
    end
end

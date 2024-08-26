classdef GP_AEKF_bounds <handle
% class for external torque estimation using augmented Kalman filter (AKF)
% Estimates and their covariances can be accessed from properties x_hat and P
% Additionally, the guaranteed estimation error bounds of the GP-AEKF are
% implemented using Lipschitz arguments here.
% Last modified: Markus Kessler 21.02.2023
% CAREFUL: THIS IS AN EARLIER DEBUGGING VERSION. IT DOES NOT CONTAIN CORRECT
% LIPSCHITZ CONSTANTS GIVING THEORETICAL GUARANTEES!
    properties
        %% Parameters for the system
        enable_learning = 1;    % enable GP based residual learning 
        sample_time;    
        Q;                      % process noise covariance matrix
        R;                      % measurement noise covariance matrix

        % augmented state vector x_hat=[theta_m, theta_s, \dot{theta}_m, \dot{theta}_s,\omega_{ext}]
        x_hat;                  % augmented state estimate
        P;                      % estimate covariance matrix

        % ellipsoid bounds (shape matrices)
        % predicted
        X_hat_pred;
        % estimated
        X_hat_est;
        % PROBABILITY LEVEL (P=1-delta)
        delta = 0.1;
        
        % nominal system model (specified by exo_system and config)
        nom_model;

        % select discretization mode: 1 = 'Euler' (faster), 2 = 'Matrix_exponential' (more accurate, slower)
        discretization_mode = 1;
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
        function initialize(obj,x0,P0,X0)
            %% initialize: Setup constant variables and initial states
            % inputs:   -x0:    initial augmented state
            %           -P0:    initial covariance
            %           -X0:    initial set of means (for error bounds)

            disp('Initializing augmented Kalman filter...');
            obj.dofs = obj.nom_model.dofs;
            obj.distorder = obj.nom_model.distorder;
            % initial state
            obj.x_hat = x0;
            % initial covariance
            obj.P = P0;
            % initial set of means
			obj.X_hat_est = X0;

            obj.state_eye = eye(obj.dofs*(4+obj.distorder));
            obj.x_hat_km1 = obj.x_hat;
            obj.jac_gp = zeros(obj.dofs*4,obj.dofs);
            obj.mu_gp = zeros(obj.dofs,1);
            obj.var_gp = zeros(obj.dofs,1);
        end
        
        function init_log_gp(obj, Xte)
            %% init_log_gp: initialize LoG-GP models for torque residual learning
            disp('Setup of MoE-LoG-GP ...');
            
            % one model for each DOF
            for i = 1:obj.nom_model.dofs
                % use mixture of GP experts LoG-GPs with error bounds
                MoE_LoG_GP(i) = MOE_LOG_GP_bounds();
                
                % set GP hyperparameters
                paramfile = ['hyperparam_dim_' num2str(i) '.mat'];
                load(paramfile, 'ls', 'sf', 'sn')
                MoE_LoG_GP(i).sigL = min(10^6,ls);
                MoE_LoG_GP(i).sigF = sf;
                MoE_LoG_GP(i).sigN = sn;
                
                % LoG-GP parameters (see MOE_LOG_GP_bounds.m for explanation)
                MoE_LoG_GP(i).pts = 100;
                MoE_LoG_GP(i).loadHyp = false;
                MoE_LoG_GP(i).wo = 100;
                MoE_LoG_GP(i).divMethod = 1;
                MoE_LoG_GP(i).xSize = 3*obj.nom_model.dofs;
                MoE_LoG_GP(i).delta = obj.delta;

                try
                    MoE_LoG_GP(i).N = 10;
                    MoE_LoG_GP(i).setupData(Xte);
                catch
                    MoE_LoG_GP(i).N = 10;
                    MoE_LoG_GP(i).setupData(Xte);
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

        function jacvec = predict(obj,tau_m)
            %% predict: Perform EKF prediction step based on nominal model
            % inputs:   -tau_m: motor torques vector

            % set nominal model state
            obj.nom_model.x = obj.x_hat;
            % Prediction from nominal model (time-continuous fwd dynamics)
            ft = obj.nom_model.fwd_dynamics(tau_m);
            % get linearized system matrices from nominal system
            [linsys.A, linsys.B, linsys.C] = obj.nom_model.fwd_dyn_jacobian(tau_m);
            
            % Lipschitz constant (of nominal model gradient and GP gradient) for linearization error
            l_mu = ones(obj.dofs*(4+obj.distorder),1);
            % process noise
            Q_GP_KF = obj.Q;
            % zero center vector for addition of zero centered ellipsoid
            zero_center = zeros(obj.dofs*(4+obj.distorder),1);

            if obj.enable_learning == 1
                % Prediction from learned GP model
                [add_gp_trans, add_gp_jac, add_gp_var, eta, beta,l_sigma,Minv,jacvec] = obj.get_GP_prediction();

                % add GP learned forward transition
                ft(obj.dofs*3+1:obj.dofs*4) = ft(obj.dofs*3+1:obj.dofs*4) - add_gp_trans;

                % add GP variance to model uncertainty Q
                Q_GP_KF(obj.dofs*3+1:obj.dofs*4,obj.dofs*3+1:obj.dofs*4) =...
                    Q_GP_KF(obj.dofs*3+1:obj.dofs*4,obj.dofs*3+1:obj.dofs*4)+add_gp_var;

                % add GP Jacobian to AKF linearization
                linsys.A(7:8,1:8) = linsys.A(7:8,1:8)-add_gp_jac;
                
                % calculate vector including standard deviations of GP for
                % inclusion of GP error bounds
            end
            % discretize linearized forward dynamics
            obj.disc_sys = obj.get_discrete_linear_sys(linsys, Q_GP_KF);
            % perform estimate prediction step with Euler integration
            obj.x_hat_pred = obj.x_hat + obj.sample_time*ft;
			
            % calculate confidence set from set of means and sigma bounds
            % on P_k-1
            [~,conf_set] = calculate_confidence_set(obj.x_hat,obj.X_hat_est,zero_center,obj.P,obj.delta);
            % compute linearization error rectangle from Lipschitz
            % constants
            [ub_EKF_lin, ub_GP_sigma] = compute_remainder_overapproximations(conf_set, l_mu,l_sigma,obj.dofs);
            % get GP prediction error
            GP_err = eta;% + sqrt(beta) .* ub_GP_sigma;
            % multiply pure GP error with sample time and inverted inertia
            % matrix (error is considered in fwd dynamics, where GP mean mu is
            % added as obj.sample_time.*Minv*mu
            GP_err = obj.sample_time.*Minv*GP_err;
            b_GP_err = zero_center;
            % fill GP error to augmented state dimension by appending zeros
            b_GP_err(obj.dofs*3+1:obj.dofs*4) = b_GP_err(obj.dofs*3+1:obj.dofs*4) + GP_err;
            
            % calculate ellipsoid bounding GP errors
            X_GP_err = ellipsoid_from_rectangle(b_GP_err);
            % calculate ellipsoid bounding linearization errors
            X_lin_err = ellipsoid_from_rectangle(ub_EKF_lin);
            
            % add error ellipsoids by Minkowski sum. Ellispoid outer
            % approximation is applied here.
            [~, X_pred_err] = sum_two_ellipsoids(zero_center,X_GP_err,zero_center,X_lin_err);
            
            % Linearly propagated ellipsoid
            X_pred_lin = obj.disc_sys.A * obj.X_hat_est * obj.disc_sys.A';
            % Sum linearly propagated ellipsoid with error ellispoid
            [~,obj.X_hat_pred] = sum_two_ellipsoids(obj.x_hat_pred,X_pred_lin,zero_center,X_pred_err);
            % calculate prediction uncertainty
            obj.P_pred = obj.disc_sys.A * obj.P * obj.disc_sys.A' + obj.disc_sys.Q;
        end

		function [add_gp_trans, add_gp_jac, add_gp_var,etavec,betavec,l_sigma,Minv,jacvec] = get_GP_prediction(obj)
            %% get_GP_prediction: Predict torque residual from learned GP model
            % outputs:   -add_gp_trans: Forward dynamics portion calculated
            %                           from learned residual torque times inverted inertia 
            %            -add_gp_jac: Jacobian of add_gp_trans
            %            -add_gp_var: Variance of add_gp_trans (from GP variance) 
            %            -etavec: Array of GP prediction errors
            %            -betavec: Array of beta coefficients for GP errors
            %            -lsigma: Array of Lipschitz constants for GP
            %                       standard deviation
            %            -Minv: inverted inertia Matrix
            %            -jacvec: Array containing 1-D gradient arrays of
            %                       GP mean gradient

            % calculate GP input = (q,qd,qdd)
            q = obj.nom_model.Nminv * obj.x_hat(1:obj.dofs)+obj.x_hat(obj.dofs+1:obj.dofs*2);
            qd = obj.nom_model.Nminv * obj.x_hat(2*obj.dofs+1:obj.dofs*3)+obj.x_hat(3*obj.dofs+1:obj.dofs*4);
            theta_m_dd_num = (obj.x_hat(2*obj.dofs+1:obj.dofs*3)-obj.x_hat_km1(2*obj.dofs+1:obj.dofs*3))./obj.sample_time;
            theta_s_dd_num = (obj.x_hat(3*obj.dofs+1:obj.dofs*4)-obj.x_hat_km1(3*obj.dofs+1:obj.dofs*4))./obj.sample_time;
            qdd = obj.nom_model.Nminv * theta_m_dd_num+theta_s_dd_num;
            
            % GP input
            x_star = [q; qd; qdd];
            etavec = zeros(obj.dofs,1);
            betavec = zeros(obj.dofs,1);
            l_sigma = zeros(obj.dofs,1);
            for dim = 1:obj.dofs
                % approximate transition residual from LoG-GP
                [mu,sig,jacvec,eta, beta,Lsig] = obj.loggp(dim).predict(x_star);
                % mean approximation
                obj.mu_gp(dim) = mu;
                % approximation variance
                obj.var_gp(dim) = sig;
                % jacobian vector of GP mean to be added to AKF Jacobian (repeat entries since Jacobian
                % is calculated with respect to q,q_dot in the GP, not
                % theta_m,theta_s,theta_m_dot,theta_s_dot. However,
                % gradients with respect to theta_s are identical to
                % gradients with respect to q. Gradients with respect
                % to theta_m are identical to Nm times the gradient with
                % respect to q due to q=Nm^-1 * theta_m + theta_s.
                jac_theta = repmat(jacvec(1:obj.dofs),obj.dofs,1);
                jac_theta_d = repmat(jacvec(obj.dofs+1:obj.dofs*2),obj.dofs,1);
                obj.jac_gp(:,dim) = [jac_theta(:); jac_theta_d(:)];
                
                % GP error
                etavec(dim) = eta;
                % GP error coefficient
                betavec(dim) = beta;
                % Standard deviation Lipschitz constants
                l_sigma(dim) = Lsig;
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

            % added GP uncertainty to Q
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
            
            % estimated set of means (for error bounds)
            obj.X_hat_est = (obj.state_eye - K_KF * obj.disc_sys.C) * obj.X_hat_pred * (obj.state_eye - K_KF * obj.disc_sys.C)';
            disp(sqrt(obj.X_hat_est(9,9)))
            % A posteriori covariance estimation update
            obj.P = (obj.state_eye - K_KF * obj.disc_sys.C) * obj.P_pred;
        end

        function disc_sys = get_discrete_linear_sys(obj, linsys, Q_KF)
            %% get_discrete_linear_sys: Discretize linearized forward dynamics for AKF
            % inputs:   -linsys:    Struct containing all time-continuous
            %                       system matrices (A,B,C,Q,R)
            % outputs:  -disc_sys:  Struct containing all discretized
            %                       system matrices (A,B,C,Q,R)

            if obj.discretization_mode == 1
                % simple Euler (zero-order hold) discretization
                disc_sys.A = obj.state_eye + linsys.A .*obj.sample_time;
                disc_sys.B = linsys.B .* obj.sample_time;
                disc_sys.Q = Q_KF .* obj.sample_time;

            elseif obj.discretization_mode == 2
                % more accurate matrix exponential discretization (slower!)
                % prepare matrices for discretization
                H1 = [linsys.A linsys.B;
                    zeros(size(linsys.B,2),size(linsys.A,2)) zeros(size(linsys.B,2),size(linsys.B,2))];
                H2 = [-linsys.A Q_KF;
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
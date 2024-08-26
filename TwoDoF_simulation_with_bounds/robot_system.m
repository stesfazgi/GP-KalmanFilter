classdef robot_system <handle
% Class defining the nominal human-robot model used by the AKF prediction.
% Set dofs, system parameters and assumed external torque dynamics to
% automatically generate all required forward & inverse dynamics.
% Last modified: Markus Kessler 16.09.2022
    
    properties
        %% Parameters for the system
        dofs;           % degrees of freedom
            
        % robot system parameters (as dofs-dimensional vectors)
        l;              % link lengths
        m;              % link masses
        r;              % link radius (cylindrical inertia)
        k;              % spring stiffness
        d;              % spring damping
        j;              % motor inertias
        
        lh;             % human limb length
        mh;             % human limb mass
        rh;             % human limb radius
        
        Nm;             % gear ratio
        Nminv;          % inverse of gear ratio matrix
        
        g = 9.81;
        
        % external torque assumed dynamics
        distorder;      % order of assumed external torque dynamics
        S;              % linear external torque dynamic matrix
        H;              % linear external torque output matrix
        
        K;              % stiffness matrix
        D;              % damping matrix
        J;              % motor inertia matrix
        C;              % output matrix
        
        use_param_human_model = false;

        % augmented state vector x=[theta_m, theta_s, \dot{theta}_m, \dot{theta}_s,\omega_{ext}]
        x;      
    end

    properties(Access = protected)
        paramvec;       % vector of system parameters 
        Mrfunc;         % Function handle for load side inertia
        Mhfunc;         % Function handle for human inertia
        Nrfunc;         % Function handle for load side nonlinearities
        Nhfunc;         % Function handle for human nonlinearities
        fwd_dyn;        % Function handle for forward dynamics
        Afunc;          % Function handle for linearized system matrix
        Bfunc;          % Function handle for linearized input matrix
        Minv_grad_func; % Function handle for gradient of inverse load side inertia
        state_eye;      % preallocated reuseable constant
    end

    methods

        function initialize(obj)
            %% initialize: Initialize system dynamics after parameters are set
            % Automatically generates symbolic dynamics matrices
            % M,C,G,N,... for given dofs and converts them to accelerated
            % matlabFunctions

            % stiffness matrix
            obj.K = diag(obj.k);
            % damping matrix (coupling)
            obj.D = diag(obj.d);
            % motor inertia matrix
            obj.J = diag(obj.j);
            % inverse of gear ratio matrix
            obj.Nminv = inv(diag(obj.Nm));
            % system output matrix (constant)
            obj.C = [eye(obj.dofs*4), zeros(obj.dofs*4, obj.distorder*obj.dofs)];

            % get symbolic load side inertia matrix
            [Mrs, symbols] = obj.get_M_robot_sym();
            % convert to symbolic to function (much faster!)
            obj.Mrfunc = matlabFunction(Mrs,'Vars',symbols);

            % get symbolic human inertia matrix
            [Mhs, symbols] = obj.get_M_human_sym();
            % convert to symbolic to function (much faster!)
            obj.Mhfunc = matlabFunction(Mhs,'Vars',symbols);

            % get gradient of inverse load side inertia (required for
            % adding GP jacobian in AKF)
            Ms = Mrs+Mhs;
            dMinvdqs = obj.get_Minv_jac_sym(Ms);
            % convert to symbolic to function (much faster!)
            obj.Minv_grad_func = matlabFunction(dMinvdqs,'Vars',symbols);

            % get symbolic robot nonlinear terms N(G,C,...) 
            [Nrs, symbols] = obj.get_N_robot_sym();
            % convert to symbolic to function (much faster!)
            obj.Nrfunc = matlabFunction(Nrs,'Vars',symbols);

            % get symbolic human nonlinear terms N(G,C,...) 
            [Nhs, symbols] = obj.get_N_human_sym();
            % convert to symbolic to function (much faster!)
            obj.Nhfunc = matlabFunction(Nhs,'Vars',symbols);
            
            Ns = Nrs + Nhs;
            % get symbolic forward dynamics of AKF state
            [fts, sym_state, sym_input]  = obj.fw_dynamics_sym(Ms, Ns);
            % convert to symbolic to function (much faster!)
            obj.fwd_dyn = matlabFunction(fts,'Vars',{sym_state, sym_input});

            % get symbolic linearizeation of forward dynamics for AKF
            [As, Bs] = obj.calc_sys_jac_sym(fts, sym_state, sym_input);
            % convert to symbolic to function (much faster!)
            obj.Afunc = matlabFunction(As, 'Vars',{sym_state, sym_input});
            obj.Bfunc = matlabFunction(Bs, 'Vars',{sym_state, sym_input});
        end

        function ft = fwd_dynamics(obj,tau_m)
            %% fwd_dynamics: Calculate forward dynamics of augmented state
            % inputs:   -tau_m: motor torques vector
            ft = obj.fwd_dyn(obj.x,tau_m);
        end

        function [Anum, Bnum, Cnum] = fwd_dyn_jacobian(obj,tau_m)
            %% fwd_dyn_jacobian: Calculate linearization of forward dynamics
            % inputs:   -tau_m: motor torques vector
            % outputs:  -Anum: Linearized system matrix
            %           -Bnum: Linearized input matrix
            %           -Cnum: Output matrix
            Anum = obj.Afunc(obj.x,tau_m);
            Bnum = obj.Bfunc(obj.x,tau_m);
            Cnum = obj.C;
        end

        function u_hat = inverse_dynamics(obj, q, qd, q_dd_meas)
            %% fwd_dyn_jacobian: Calculate passive torques from inverse dynamics without inputs
            % inputs:   -q:                 load side angles
            %           -qd:                load side velocities
            %           -qdd:               load side accelerations
            %           -theta_m_dd_meas:   motor side accelerations
            % outputs:  -u_hat: passive torques from inverse dynamics

%             q = obj.aug_state(1:obj.dofs) + obj.aug_state(obj.dofs+1:2*obj.dofs); % q = theta_m+theta_s
%             qd = obj.aug_state(2*obj.dofs+1:3*obj.dofs) + obj.aug_state(3*obj.dofs+1:4*obj.dofs); % qd = theta_m_d+theta_s_d
            M = obj.Mrfunc(q)+obj.Mhfunc(q);
            N = obj.Nrfunc(q,qd)+obj.Nhfunc(q,qd);
            u_hat = M * q_dd_meas + N;
        end
        
        function Mnum = get_M_numerical(obj)
            %% get_M_numerical: Return numerical load side inertia
            % outputs:  -Mnum: Load side inertia at current state q
            q = obj.Nminv*obj.x(1:obj.dofs)+obj.x(obj.dofs+1:obj.dofs*2);
            Mnum = obj.Mrfunc(q)+obj.Mhfunc(q);
        end

        function Minv_dq = get_Minv_deriv(obj)
            %% get_Minv_deriv: Return numerical gradient of load side inertia
            % outputs:  -Mnum: Gradient of load side inertia at current state q
            q = obj.Nminv*obj.x(1:obj.dofs)+obj.x(obj.dofs+1:obj.dofs*2);
            Minv_dq = obj.Minv_grad_func(q);
        end

        function [ft, sym_state, tau] = fw_dynamics_sym(obj, M, N)
            %% fw_dynamics_sym: Calculate symbolic forward dynamics
            % inputs:   -M:     symbolic load side inertia
            %           -N:     symbolic load side nonlinearities
            % outputs:  -ft:        symbolic forward dynamics
            %           -sym_state: state symbols
            %           -tau:       input symbols (motor torques)

            % input symbols
            tau = sym('tau_m',[obj.dofs, 1]);
            % state symbols
            theta_m = sym('theta_m',[obj.dofs, 1]);
            theta_m_d = sym('theta_m_d',[obj.dofs, 1]);
            theta_s = sym('theta_s',[obj.dofs, 1]);
            theta_s_d = sym('theta_s_d',[obj.dofs, 1]);
            q = sym('q',[obj.dofs, 1]);
            qd = sym('qd',[obj.dofs, 1]);
            omega_ext = sym('omega_ext',[obj.distorder,obj.dofs]);

            % external torque dynamics
            SCell = repmat({obj.S}, 1, obj.dofs);
            BigS = blkdiag(SCell{:});
            distdyn = BigS*omega_ext(:);
            HCell = repmat({obj.H}, 1, obj.dofs);
            BigH = blkdiag(HCell{:});
            % external torques
            taud = nonzeros(BigH*omega_ext(:));

            % get forward dynamics in symbolic variables
            q_second_derivatives = -inv(M)*( N + obj.D*theta_s_d + obj.K*theta_s - taud);
            theta_second_derivatives = -inv(obj.J)*( -obj.Nminv*(obj.D*theta_s_d + obj.K*theta_s) - tau);
            
            % forward dynamics of augmented state
            ft =...
                [theta_m_d;
                theta_s_d;
                simplify(theta_second_derivatives);
                simplify(q_second_derivatives-obj.Nminv*theta_second_derivatives);
                distdyn];
            % replace load side symbols by motor and spring symbols
            ft = subs(ft,[q, qd],[theta_s+obj.Nminv*theta_m, theta_s_d+obj.Nminv*theta_m_d]);
            % state symbols
            sym_state = [theta_m; theta_s; theta_m_d; theta_s_d; omega_ext(:)];
        end
        
        function [A, B] = calc_sys_jac_sym(obj, ft, sym_state, sym_input)
            %% calc_sys_jac_sym: Calculate symbolic linearization of forward dynamics
            % inputs:   -ft:        symbolic forward dynamics
            %           -sym_state: state symbols
            %           -sym_input: input symbols
            % outputs:  -A:         symbolic system matrix
            %           -B:         symbolic input matrix
                A = simplify(jacobian(ft, sym_state));
                B = simplify(jacobian(ft, sym_input));
        end

        function [N, symbols] = get_N_robot_sym(obj)
            %% get_N_sym: Calculate symbolic robot nonlinearities
            % outputs:  -N:         symbolic nonlinear torques
            %           -symbols:   state symbols
            q = sym('q',[obj.dofs, 1]);
            qd = sym('qd',[obj.dofs, 1]);

            if obj.dofs == 1
                G = -obj.m(1)*obj.g*obj.l(1)/2*cos(q(1));
                C = sym(0);
            else
                %% Coriolis- and centrifugal torques
                lc(1) = obj.l(1)/2; lc(2) = obj.l(2)/2;
                C = [(-obj.m(2)*obj.l(1)*lc(2)*sin(q(2)))*qd(2)^2 - 2*obj.m(2)*obj.l(1)*lc(2)*sin(q(2))*qd(1)*qd(2);
                    obj.m(2)*obj.l(1)*lc(2)*sin(q(2))*qd(1)^2];

                %% Gravitational torques
                G = [obj.m(1)*lc(1)*obj.g*cos(q(1)) + obj.m(2)*obj.g*(lc(2)*cos(q(1)+q(2))+obj.l(1)*cos(q(1)));
                    obj.m(2)*obj.g*lc(2)*cos(q(1)+q(2))];
            end
            % Combine nonlinearities to N
            N = C + G;
            symbols = {q,qd};
        end

        function [M, symbols] = get_M_robot_sym(obj)
            %% get_M_sym: Calculate symbolic load side inertia
            % outputs:  -M:         symbolic load side inertia
            %           -symbols:   state symbols
            q = sym('q',[obj.dofs, 1]);

            if obj.dofs == 1
                M = m1*(l1/2)^2 + 1/4*m1*r1^2+1/12*m1*l1^2;
                M = sym(1.32e-2);
            else

            % centers of mass
            lc(1) = obj.l(1)/2; lc(2) = obj.l(2)/2;
            
            % cylindrical constant inertia
            I(1) = 1/4*obj.m(1)*obj.r(1)^2 + 1/12*obj.m(1)*obj.l(1)^2;
            I(2) = 1/4*obj.m(2)*obj.r(2)^2 + 1/12*obj.m(2)*obj.l(2)^2;
            %% Load side inertia
            M = [obj.m(1)*lc(1)^2 + I(1) + obj.m(2)*(obj.l(1)^2+lc(2)^2+2*obj.l(1)*lc(2)*cos(q(2))) + I(2), obj.m(2)*(lc(2)^2+obj.l(1)*lc(2)*cos(q(2))) + I(2);
                obj.m(2)*(lc(2)^2+obj.l(1)*lc(2)*cos(q(2))) + I(2),                       obj.m(2)*lc(2)^2 + I(2)           ];
            end
            symbols = {q};   
        end

        function [N, symbols] = get_N_human_sym(obj)
            %% get_N_sym: Calculate symbolic human nonlinearities
            % outputs:  -N:         symbolic nonlinear torques
            %           -symbols:   state symbols
            
            if obj.dofs == 1
                G = -obj.m(1)*obj.g*obj.l(1)/2*cos(q(1));
                C = sym(0);
            else
                %% Coriolis- and centrifugal torques
                lc(1) = obj.lh(1)/2; lc(2) = obj.lh(2)/2;
                q = sym('q',[obj.dofs, 1]);
                qd = sym('qd',[obj.dofs, 1]);
                C = [(-obj.mh(2)*obj.lh(1)*lc(2)*sin(q(2)))*qd(2)^2 - 2*obj.mh(2)*obj.lh(1)*lc(2)*sin(q(2))*qd(1)*qd(2);
                    obj.mh(2)*obj.lh(1)*lc(2)*sin(q(2))*qd(1)^2];

                %% Gravitational torques
                G = [obj.mh(1)*lc(1)*obj.g*cos(q(1)) + obj.mh(2)*obj.g*(lc(2)*cos(q(1)+q(2))+obj.lh(1)*cos(q(1)));
                    obj.mh(2)*obj.g*lc(2)*cos(q(1)+q(2))];
            end
            % Combine nonlinearities to N
            N = C + G;
            if obj.use_param_human_model == false
                N = N*0;
            end
            symbols = {q,qd};
        end
        
        function [M, symbols] = get_M_human_sym(obj)
            %% get_M_human_sym: Calculate symbolic human inertia
            % outputs:  -M:         symbolic human inertia
            %           -symbols:   state symbols
            q = sym('q',[obj.dofs, 1]);
            % centers of mass
            lc(1) = obj.lh(1)/2; lc(2) = obj.lh(2)/2;
            
            % cylindrical constant inertia
            I(1) = 1/4*obj.mh(1)*obj.rh(1)^2 + 1/12*obj.mh(1)*obj.lh(1)^2;
            I(2) = 1/4*obj.mh(2)*obj.rh(2)^2 + 1/12*obj.mh(2)*obj.lh(2)^2;
            %% Human inertia
            M = [obj.mh(1)*lc(1)^2 + I(1) + obj.mh(2)*(obj.lh(1)^2+lc(2)^2+2*obj.lh(1)*lc(2)*cos(q(2))) + I(2), obj.mh(2)*(lc(2)^2+obj.lh(1)*lc(2)*cos(q(2))) + I(2);
                obj.mh(2)*(lc(2)^2+obj.lh(1)*lc(2)*cos(q(2))) + I(2),                       obj.mh(2)*lc(2)^2 + I(2)           ];
            symbols = {q};
            if obj.use_param_human_model == false
                M = M*0;
            end
        end

        function dMinvdq = get_Minv_jac_sym(obj,M)
            %% get_Minv_jac_sym: Calculate symbolic load side inertia
            % inputs:   -M:         symbolic load side inertia
            % outputs:  -dMinvdq:   gradient of inverted load side inertia
            q = sym('q',[obj.dofs, 1]);     
            Minv = inv(M);
            % calculate gradient
            dMinvdq = diff(Minv,q(2));
        end
    end
end


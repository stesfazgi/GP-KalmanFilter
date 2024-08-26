function tau = ctc_controller(exo_sys, current_x, theta_second_deriv, M, n, x_ref, i, sim_phase)
    J = diag(exo_sys.j);
    [Mh, tau_h] = human_arm_model(exo_sys, current_x);
    u_CTC = n + diag(exo_sys.Nm)*(J*theta_second_deriv) + tau_h;

    if sim_phase==3
        K_p = [1000 0;0 1000];
        K_d = [1000 0;0 1000];
        K_p = [50 0;0 50];
        K_d = [50 0;0 50];
    else
        K_p = [50 0;0 50];
        K_d = [50 0;0 50];
%         K_p = [1 0;0 1];
%         K_d = [1 0;0 1];
    end

    u_PD = -K_p * ([x_ref.q1(i); x_ref.q2(i)]  - [current_x(1); current_x(5)]) - K_d * ([x_ref.q1_first_deriv(i); x_ref.q2_first_deriv(i)]  - [current_x(2); current_x(6)]);

    %         disp('tracking error e')
    %         [x_ref.q1(i); x_ref.q2(i)]  - [current_x(1); current_x(5)];
    %         disp('tracking error e_d')
    %         [x_ref.q1_first_deriv(i); x_ref.q2_first_deriv(i)]  - [current_x(2); current_x(6)];

    qddref = [x_ref.q1_second_deriv(i); x_ref.q1_second_deriv(i)];
    tau = exo_sys.Nminv*((M+Mh)*(qddref - u_PD) + u_CTC);
end


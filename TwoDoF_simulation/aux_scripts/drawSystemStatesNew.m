function fig_system_states = drawSystemStatesNew(system_states_table,x_ref, tau_d, plots)
%
    % Define primary TUM colors
    tum_blue = [0 1 189]./255;
    % Define secondary TUM colors
    tum_dark_gray = [51, 51, 51]./255; tum_gray = [128, 128, 128]./255; tum_light_gray = [204, 204, 204]./255;
    % Define accent TUM colors
    tum_green = [162, 173, 0]./255; tum_orange = [227, 114, 34]./255; tum_light_blue = [152, 198, 234]./255;

%     fig_system_states = figure('Name', 'System States');
%     fig_system_states.Color = [1 1 1];
    
    fig_system_states = 0;
    
    length_t = length(system_states_table.t);
    plot_step = 1;
%     data = load('NominalDampedControl_SystemStates.mat');
%     tau_nom1 = data.system_states.tau1;
%     tau_nom2 = data.system_states.tau2;
    
    t = system_states_table.t;
    labels_angles = [];
    labels_torques = [];
        
    % Plot first link angle
    figure("Name", "q1");
%     subplot(2,2,1);
    plot(t,system_states_table.q1,'Color',tum_blue,'LineWidth',1.7); labels_angles = [labels_angles; "Simulated trajectory"]; hold on;
    plot(t,x_ref.q1,'Color',tum_orange,'LineStyle','--','LineWidth',1.2); labels_angles = [labels_angles; "Reference trajectory"];

    title('Link Angle First Joint ($q_1$)','Interpreter','latex'); xlabel('t in s'); ylabel('$q_1$ in rad','Interpreter','latex');
    grid on;
    
    % Plot second link angle
    figure("Name", "q2");
%     subplot(2,2,2);
    plot(t,system_states_table.q2,'Color',tum_blue,'LineWidth',1.7); hold on;
    plot(t,x_ref.q2,'Color',tum_orange,'LineStyle','--','LineWidth',1.2);

    title('Link Angle Second Joint ($q_2$)','Interpreter','latex'); xlabel('t in s'); ylabel('$q_2$ in rad','Interpreter','latex');
    lgd_angles = legend(labels_angles);
%     lgd_angles.Position = [0.454427083333333,0.885290889132821,0.124739583333333,0.08397365532382];
    grid on;
    
    % Plot torque for first motor
    figure("Name","tau1");
%     subplot(2,2,3);
    plot(t,system_states_table.tau1,'Color',tum_green,'LineWidth',1.7); hold on;
%     plot(t,tau_nom1,'Color',tum_orange,'LineStyle','--','LineWidth',1.2);
    labels_torques = ["Torques with Safety Filter"];%; "Nominal Controller"];
    title('Exerted Torque First Joint ($\tau_1$)','Interpreter','latex'); xlabel('t in s'); ylabel('$\tau_1$ in Nm','Interpreter','latex');
    grid on;
    
    % Plot torque for second motor
    figure("Name","tau2");
%     subplot(2,2,4);
    plot(t,system_states_table.tau2,'Color',tum_green,'LineWidth',1.7); hold on;
%     plot(t,tau_nom2,'Color',tum_orange,'LineStyle','--','LineWidth',1.2);
    title('Exerted Torque Second Joint ($\tau_2$)','Interpreter','latex'); xlabel('t in s'); ylabel('$\tau_2$ in Nm','Interpreter','latex');
    lgd_torques = legend(labels_torques);
%     lgd_torques.Position = [0.44921875,0.127332601536773,0.137239583333333,0.031833150384193];
    grid on;
    
    %% plot AKF estimates
    xt = system_states_table.t(1:plot_step:length_t);
    % link 1 states
    figure("Name","Link 1 AKF state estimates");
    title('Link 1 state estimates');

    subplot(2,2,1)
    plot(xt,plots.AKF_ubs(1:plot_step:length_t,1),'c-.','LineWidth',1.2); hold on;
    plot(xt,plots.AKF_lbs(1:plot_step:length_t,1),'c-.','LineWidth',1.2); hold on;
    plot(xt,system_states_table.theta1(1:plot_step:length_t),'r','LineWidth',1.5); hold on;
    plot(xt,plots.x_KF_plot(1:plot_step:length_t, 1),'b--','LineWidth',1.5); hold on;
    legend('$2\sigma$ ucb','$2\sigma$ lcb','$\theta_m$ true','$\theta_m$ estimate','Interpreter', 'latex');
    grid on;
    subplot(2,2,2)
    plot(xt,plots.AKF_ubs(1:plot_step:length_t,5),'c-.','LineWidth',1.2); hold on;
    plot(xt,plots.AKF_lbs(1:plot_step:length_t,5),'c-.','LineWidth',1.2); hold on;
    plot(xt,system_states_table.theta1_dot(1:plot_step:length_t),'r','LineWidth',1.5); hold on;
    plot(xt,plots.x_KF_plot(1:plot_step:length_t, 5),'b--','LineWidth',1.5); hold on;
    legend('$2\sigma$ ucb','$2\sigma$ lcb','$\dot{\theta}_m$ true','$\dot{\theta}_m$ estimate','Interpreter', 'latex');
    grid on;
    subplot(2,2,3)
    plot(xt,plots.AKF_ubs(1:plot_step:length_t,3),'c-.','LineWidth',1.2); hold on;
    plot(xt,plots.AKF_lbs(1:plot_step:length_t,3),'c-.','LineWidth',1.2); hold on;
    plot(xt,system_states_table.q1(1:plot_step:length_t) - system_states_table.theta1(1:plot_step:length_t),'r','LineWidth',1.5); hold on;
    plot(xt,plots.x_KF_plot(1:plot_step:length_t, 3),'b--','LineWidth',1.5); hold on;
    legend('$2\sigma$ ucb','$2\sigma$ lcb','$x_s$ true','$x_s$ estimate','Interpreter', 'latex');
    grid on;
    subplot(2,2,4)
    plot(xt,plots.AKF_ubs(1:plot_step:length_t,7),'c-.','LineWidth',1.2); hold on;
    plot(xt,plots.AKF_lbs(1:plot_step:length_t,7),'c-.','LineWidth',1.2); hold on;
    plot(xt,system_states_table.q1_dot(1:plot_step:length_t) - system_states_table.theta1_dot(1:plot_step:length_t),'r','LineWidth',1.5); hold on;
    plot(xt,plots.x_KF_plot(1:plot_step:length_t, 7),'b--','LineWidth',1.5); hold on;
    legend('$2\sigma$ ucb','$2\sigma$ lcb','$\dot{x}_s$ true','$\dot{x}_s$ estimate','Interpreter', 'latex');
    grid on;
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    
    % link 2 states
    figure("Name","Link 2 AKF state estimates");
    title('Link 2 state estimates');

    subplot(2,2,1)
    plot(xt,plots.AKF_ubs(1:plot_step:length_t,2),'c-.','LineWidth',1.2); hold on;
    plot(xt,plots.AKF_lbs(1:plot_step:length_t,2),'c-.','LineWidth',1.2); hold on;
    plot(xt,system_states_table.theta2(1:plot_step:length_t),'r','LineWidth',1.5); hold on;
    plot(xt,plots.x_KF_plot(1:plot_step:length_t, 2),'b--','LineWidth',1.5); hold on;
    legend('$2\sigma$ ucb','$2\sigma$ lcb','$\theta_m$ true','$\theta_m estimate$','Interpreter', 'latex');
    grid on;
    subplot(2,2,2)
    plot(xt,plots.AKF_ubs(1:plot_step:length_t,6),'c-.','LineWidth',1.2); hold on;
    plot(xt,plots.AKF_lbs(1:plot_step:length_t,6),'c-.','LineWidth',1.2); hold on;
    plot(xt,system_states_table.theta2_dot(1:plot_step:length_t),'r','LineWidth',1.5); hold on;
    plot(xt,plots.x_KF_plot(1:plot_step:length_t, 6),'b--','LineWidth',1.5); hold on;
    legend('$2\sigma$ ucb','$2\sigma$ lcb','$\dot{\theta}_m$ true','$\dot{\theta}_m$ estimate','Interpreter', 'latex');
    grid on;
    subplot(2,2,3)
    plot(xt,plots.AKF_ubs(1:plot_step:length_t,4),'c-.','LineWidth',1.2); hold on;
    plot(xt,plots.AKF_lbs(1:plot_step:length_t,4),'c-.','LineWidth',1.2); hold on;
    plot(xt,system_states_table.q2(1:plot_step:length_t) - system_states_table.theta2(1:plot_step:length_t),'r','LineWidth',1.5); hold on;
    plot(xt,plots.x_KF_plot(1:plot_step:length_t, 4),'b--','LineWidth',1.5); hold on;
    legend('$2\sigma$ ucb','$2\sigma$ lcb','$x_s$ true','$x_s$ estimate','Interpreter', 'latex');
    grid on;
    subplot(2,2,4)
    plot(xt,plots.AKF_ubs(1:plot_step:length_t,8),'c-.','LineWidth',1.2); hold on;
    plot(xt,plots.AKF_lbs(1:plot_step:length_t,8),'c-.','LineWidth',1.2); hold on;
    plot(xt,system_states_table.q2_dot(1:plot_step:length_t) - system_states_table.theta2_dot(1:plot_step:length_t),'r','LineWidth',1.5); hold on;
    plot(xt,plots.x_KF_plot(1:plot_step:length_t, 8),'b--','LineWidth',1.5); hold on;
    legend('$2\sigma$ ucb','$2\sigma$ lcb','$\dot{x}_s$ true','$\dot{x}_s$ estimate','Interpreter', 'latex');
    grid on;
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    
    % external torque estimate
    figure("Name","External torque estimates");
    
    x2 = [xt', fliplr(xt')];
    conf_region1 = [plots.AKF_lbs(1:plot_step:length_t,9)', fliplr(plots.AKF_ubs(1:plot_step:length_t,9)')];
    subplot(2,1,1)
    title('Link1 external torque estimation');
    fill(x2, conf_region1,[0, 1, 1],'EdgeColor','none'); hold on;
    plot(xt,plots.tau_d_KF_hat_plot(1:plot_step:length_t,1),'b','LineWidth',2.5);hold on;
%     plot(xt,plots.AKF_ubs(1:plot_step:length_t,9),'-.','Color', [1, 1, 0, 1],'LineWidth',1.5); hold on;
%     plot(xt,plots.AKF_lbs(1:plot_step:length_t,9),'-.','Color', [1, 1, 0, 1],'LineWidth',1.5); hold on;
    plot(xt,tau_d(1:plot_step:length_t,1),'--','Color',tum_orange,'LineWidth',3); hold on;
    legend('$90\%$ error bound','AKF estimate $\hat{\tau}_{ext,1}$','true external torque $\tau_{ext,1}$','Interpreter', 'latex');
    xlabel('time [s]')
    ylabel('torque [Nm]');

    grid on;
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    
    conf_region2 = [plots.AKF_ubs(1:plot_step:length_t,10)', fliplr(plots.AKF_lbs(1:plot_step:length_t,10)')];
    subplot(2,1,2)
    title('Link2 external torque estimation');
    fill(x2, conf_region2,[0, 1, 1],'EdgeColor','none'); hold on;
    plot(xt,plots.tau_d_KF_hat_plot(1:plot_step:length_t,2),'b','LineWidth',2.5);hold on;
%     plot(xt,plots.AKF_ubs(1:plot_step:length_t,10),'-.','Color', [1, 1, 0, 1],'LineWidth',1.5); hold on;
%     plot(xt,plots.AKF_lbs(1:plot_step:length_t,10),'-.','Color', [1, 1, 0, 1],'LineWidth',1.5); hold on;

    plot(xt,tau_d(1:plot_step:length_t,2),'--','Color',tum_orange,'LineWidth',3); hold on;
    legend('$90\%$ error bound','AKF estimate $\hat{\tau}_{ext,2}$','true external torque $\tau_{ext,2}$','Interpreter', 'latex');
    xlabel('time [s]')
    ylabel('torque [Nm]');
    grid on;
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    
    % LoG GP outputs
%     figure("Name","LoG-GP mean predictions");
%     subplot(2,1,1)
%     title('Link1 LoG-GP mean predictions');
%     plot(xt,mugp_plot(1:plot_step:length_t,1),'LineWidth',2.5); hold on;
%     plot(xt,residual_true_plot(1:plot_step:length_t,1),'LineWidth',2.5); hold on;
% %     plot(system_states.t(1:plot_step:length_t),errrestrueplot(1:plot_step:length_t),'LineWidth',2.5); hold on;
% %     plot(system_states.t(1:plot_step:length_t),mugp_avg_plot(1:plot_step:length_t),'LineWidth',2.5); hold on;
% %     legend('mu4','residual true', 'residual wo \tau_{ext}','mu4 moving avg');
%     legend('mu1','true residual');
%     grid on;
%     set(findall(gcf,'-property','FontSize'),'FontSize',16)
%     
%     subplot(2,1,2)
%     title('Link2 LoG-GP mean predictions');
%     plot(xt,mugp_plot(1:plot_step:length_t,2),'LineWidth',2.5); hold on;
%     plot(xt,residual_true_plot(1:plot_step:length_t,2),'LineWidth',2.5); hold on;
% %     plot(system_states.t(1:plot_step:length_t),errrestrueplot(1:plot_step:length_t),'LineWidth',2.5); hold on;
% %     plot(system_states.t(1:plot_step:length_t),mugp_avg_plot(1:plot_step:length_t),'LineWidth',2.5); hold on;
% %     legend('mu4','residual true', 'residual wo \tau_{ext}','mu4 moving avg');
%     legend('mu2','true residual');
%     grid on;
%     set(findall(gcf,'-property','FontSize'),'FontSize',16)
    
%     figure("Name","LoG-GP prediction variance");
%     subplot(2,1,1)
%     title('Link1 LoG-GP prediction variance');
%     plot(xt,vargp_plot(1:plot_step:length_t,1),'LineWidth',2.5); hold on;
%     title('LoG-GP Variance');
%     grid on;
%     set(findall(gcf,'-property','FontSize'),'FontSize',16)
%     legend('var1');
%     
%     subplot(2,1,2)
%     title('Link2 LoG-GP prediction variance');
%     plot(xt,vargp_plot(1:plot_step:length_t,2),'LineWidth',2.5); hold on;
%     title('LoG-GP Variance');
%     grid on;
%     set(findall(gcf,'-property','FontSize'),'FontSize',16)
%     legend('var2');
   
    % count of active local GP models
    
%     figure("Name","Count of locally active LoG-GP models");
%     plot(xt,nractivegps_plot(1:plot_step:length_t),'LineWidth',2.5); hold on;
%     title('Number of active local GPs');
%     grid on;
%     set(findall(gcf,'-property','FontSize'),'FontSize',16)
%     legend;
    
%     AKF estimate variance
    figure("Name","AKF estimate variances");
    plot(xt,plots.AKF_cov_plot(1:plot_step:length_t,:),'LineWidth',2.5); hold on;
    title('Augmented KF covariance estimate');
    legend;
    grid on;
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    
    % link positions
%     figure("Name","link positions");
%     plot(xt,plots.q(1:plot_step:length_t,:),'LineWidth',2.5); hold on;
%     title('link positions');
%     legend;
%     grid on;
%     set(findall(gcf,'-property','FontSize'),'FontSize',16)
%     
%     % link velocities
%     figure("Name","link velocities");
%     plot(xt,plots.qd(1:plot_step:length_t,:),'LineWidth',2.5); hold on;
%     title('link velocities');
%     legend;
%     grid on;
%     set(findall(gcf,'-property','FontSize'),'FontSize',16)
%     
%     % link accelerations
%     figure("Name","link accelerations");
%     plot(xt,plots.qdd(1:plot_step:length_t,:),'LineWidth',2.5); hold on;
%     title('link accelerations');
%     legend;
%     grid on;
%     set(findall(gcf,'-property','FontSize'),'FontSize',16)

    % GP jacobian
%     figure("Name","gp jacobian");
%     plot(xt,plots.jacvec(1:plot_step:length_t,:),'LineWidth',2.5); hold on;
%     title('gp jacobian');
%     legend;
%     grid on;
%     set(findall(gcf,'-property','FontSize'),'FontSize',16)
end
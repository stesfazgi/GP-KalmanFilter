function [fig_simulation,video_handle] = drawAnimatedSystem(t,x,x_ref,l,tau, tau_d,plots)
%
    % Anonymous function to comfortely draw an arrow
    drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:} );

%     tau_d = zeros(size(x));
%     tau_d(500:end)=1;
%     
%     tau = zeros(size(x));
%     tau(1:700)=1;
    
    % Define primary TUM colors
    tum_blue = [0 1 189]./255;
    % Define secondary TUM colors
    tum_dark_gray = [51, 51, 51]./255; tum_gray = [128, 128, 128]./255; tum_light_gray = [190, 190, 190]./255;
    
    % Define accent TUM colors
    tum_green = [162, 173, 0]./255; tum_orange = [227, 114, 34]./255; tum_light_blue = [152, 198, 234]./255;
    
    human_color = [250, 128, 114]./255;
    human_joint_color = [200, 30, 10]./255;
    human_body_color = [200, 200, 200]./255;
    active_color= [0, 205, 0]./255;
    purple = [128, 0, 128]./255;
    % Load system parameters
    l1 = l(1); l2 = l(2);
    l = l1+l2;
    radius_joints = 0.12;
    radius_man = radius_joints/2;
    l1h = l1 - 2*radius_joints;
    l2h = l2;
    lhand = 0.35;
    radiushand = radius_man*1.3;
%     widthupper = 0.1;
    necksize = radius_joints*3;
    head_radius = radius_joints*3;
    
    y_body_pos = -l1-l2+radius_joints;
    x_body_pos = (-l1-l2)/4;
    y_neck_pos = radius_joints;
    x_neck_pos = -necksize*0.75/2;
    x_head_pos = -head_radius;
    y_head_pos = necksize*1.2;
    
    x1_pos_strap1 = (l1/4)*cos(x(:,1));
    x2_pos_strap1 = x1_pos_strap1+l1/10*cos(x(:,1));
    y1_pos_strap1 = (l1/4)*sin(x(:,1)); 
    y2_pos_strap1 = y1_pos_strap1+l1/10*sin(x(:,1)); 
    
    x1_pos_strap2 = (2*l1/3)*cos(x(:,1));
    x2_pos_strap2 = x1_pos_strap2+l1/10*cos(x(:,1));
    y1_pos_strap2 = (2*l1/3)*sin(x(:,1)); 
    y2_pos_strap2 = y1_pos_strap2+l1/10*sin(x(:,1)); 
    
    x1_pos_strap3 = l1*cos(x(:,1))+l2/4*cos(x(:,1)+x(:,5));
    x2_pos_strap3 = x1_pos_strap3+l2/10*cos(x(:,1)+x(:,5));
    y1_pos_strap3 = l1*sin(x(:,1))+l2/4*sin(x(:,1)+x(:,5)); 
    y2_pos_strap3 = y1_pos_strap3+l2/10*sin(x(:,1)+x(:,5));
    
    x1_pos_strap4 = l1*cos(x(:,1))+(l2-radius_man*2-l2/10)*cos(x(:,1)+x(:,5));
    x2_pos_strap4 = x1_pos_strap4+l2/10*cos(x(:,1)+x(:,5));
    y1_pos_strap4 = l1*sin(x(:,1))+(l2-radius_man*2-l2/10)*sin(x(:,1)+x(:,5)); 
    y2_pos_strap4 = y1_pos_strap4+l2/10*sin(x(:,1)+x(:,5));
    
    x_circ = radius_joints*cos(2*pi*(0:0.01:1)); 
    y_circ = radius_joints*sin(2*pi*(0:0.01:1));
    
    % Compute joint positions upfront
    x_pos_link0 = radius_joints*cos(x(:,1));
    y_pos_link0 = radius_joints*sin(x(:,1));
    
    hx_pos_link0 = -2*radius_joints*sin(x(:,1));
    hy_pos_link0 = 2*radius_joints*cos(x(:,1));
    
    x_pos_link1 = (l1)*cos(x(:,1)); 
    y_pos_link1 = (l1)*sin(x(:,1)); 
    
    hx_pos_link1 = (l1h)*cos(x(:,1))-2*radius_joints*sin(x(:,1)); 
    hy_pos_link1 = (l1h)*sin(x(:,1))+2*radius_joints*cos(x(:,1));
    
    x_pos_link2 = l1*cos(x(:,1)) + (l2-0.5*radius_joints)*cos(x(:,1)+x(:,5));
    y_pos_link2 = l1*sin(x(:,1)) + (l2-0.5*radius_joints)*sin(x(:,1)+x(:,5));
    
    hx_pos_link2 = hx_pos_link1 + (l2h)*cos(x(:,1)+x(:,5));%-radius_joints*sin(x(:,1)+x(:,5));
    hy_pos_link2 = hy_pos_link1 + (l2h)*sin(x(:,1)+x(:,5));%+radius_joints*cos(x(:,1)+x(:,5));
    
    hx_pos_hand = l1*cos(x(:,1))+l2h*cos(x(:,1)+x(:,5));
    hy_pos_hand = l1*sin(x(:,1))+l2h*sin(x(:,1)+x(:,5));
    
    hx_pos_finger = hx_pos_hand+lhand*0.6*cos(x(:,1)+x(:,5));
    hy_pos_finger = hy_pos_hand+lhand*0.6*sin(x(:,1)+x(:,5));
    
    hx_pos_thumb = hx_pos_hand+lhand*0.65*cos(x(:,1)+x(:,5)+pi/6);
    hy_pos_thumb = hy_pos_hand+lhand*0.65*sin(x(:,1)+x(:,5)+pi/6);
    
    hx_pos_tip = hx_pos_finger+lhand*0.5*cos(x(:,1)+x(:,5)+pi/6);
    hy_pos_tip = hy_pos_finger+lhand*0.5*sin(x(:,1)+x(:,5)+pi/6);

    manipulator_positions.x = l1*cos(x(:,1))+l2*cos(x(:,1)+x(:,5));
    manipulator_positions.y = l1*sin(x(:,1))+l2*sin(x(:,1)+x(:,5));
  
    % Video settings
    video_handle = VideoWriter('animation');
    video_handle.FrameRate = 15;
    open(video_handle);
    
    fig_simulation = figure('Name', 'Simulation', 'Position', [0 0 2560 1440]);
    fig_simulation.Color = [1 1 1];
    labels = [];
    
    % Initialize objects
    txt1 = text(0.9,0.9,"");  txt1.Units = "normalized";
    txt2 = text(0.9,0.8,"");  txt2.Units = "normalized";
    txt3 = text(0.9,0.7,"");  txt3.Units = "normalized";

    % Animate system
    tic;
    i = 0; increment = 10;
    while i + increment < length(t)   
        i = i + increment;
        
        % Draw robot at current position
        figure(fig_simulation); clf(fig_simulation);
        subplot(2, 5, [1 2 3 6 7 8])
        % plot human torso
        
        if tau_d(i) == 0
          body_color = human_body_color;
        else
          body_color = active_color;
        end
        
        if tau(i) == 0
          man_color = tum_light_blue;
        else
          man_color = active_color;
        end
        
        rectangle('Position',[x_neck_pos y_neck_pos necksize*0.75 necksize],'EdgeColor',body_color,'FaceColor',[1 1 1],'LineWidth',2); hold on;
        rectangle('Position',[x_head_pos y_head_pos 2*head_radius 2*head_radius],'Curvature', [1 1], 'EdgeColor',body_color,'FaceColor',[1 1 1],'LineWidth',2); hold on;
        rectangle('Position',[x_body_pos y_body_pos (l1+l2)/2 (l1+l2)*1.05],'Curvature',[0.5,1],'EdgeColor',body_color,'FaceColor',[1 1 1],'LineWidth',2); hold on;
        
        % Plot links and trajectory
        plot([0, x_pos_link1(i)], [0, y_pos_link1(i)], 'Color', human_color, 'LineWidth', 30, 'HandleVisibility', 'off'); hold on; % First Link
        plot([l1*cos(x(i,1)), hx_pos_hand(i)], [l1*sin(x(i,1)), hy_pos_hand(i)],...
          'Color', human_color, 'LineWidth', 25, 'HandleVisibility', 'off'); hold on; % Second Link
        plot([hx_pos_hand(i), hx_pos_finger(i)], ...
          [hy_pos_hand(i), hy_pos_finger(i)],'-o', 'MarkerSize',2, 'MarkerEdgeColor',human_color,...
          'MarkerFaceColor',human_color, 'Color', human_color, 'LineWidth', 14, 'HandleVisibility', 'off');
        plot([hx_pos_hand(i), hx_pos_thumb(i)], ...
          [hy_pos_hand(i), hy_pos_thumb(i)],'-o', 'MarkerSize',2, 'MarkerEdgeColor',human_color,...
          'MarkerFaceColor',human_color, 'Color', human_color, 'LineWidth', 11, 'HandleVisibility', 'off');
        plot([hx_pos_finger(i), hx_pos_tip(i)], ...
          [hy_pos_finger(i), hy_pos_tip(i)],'-o', 'MarkerSize',2, 'MarkerEdgeColor',human_color,...
          'MarkerFaceColor',human_color, 'Color', human_color, 'LineWidth', 11, 'HandleVisibility', 'off');
        
        
        plot([x_pos_link0(i), x_pos_link1(i)], [y_pos_link0(i), y_pos_link1(i)], 'Color', tum_light_gray, ...
          'LineWidth', 6, 'HandleVisibility', 'off'); hold on; % First Link
        plot([l1*cos(x(i,1))+radius_joints*cos(x(i,1)+x(i,5)), x_pos_link2(i)], ...
          [l1*sin(x(i,1))+radius_joints*sin(x(i,1)+x(i,5)), y_pos_link2(i)], 'Color', tum_light_gray, 'LineWidth', 6, 'HandleVisibility', 'off'); hold on; % Second Link
        
        plot([x1_pos_strap1(i), x2_pos_strap1(i)], [y1_pos_strap1(i), y2_pos_strap1(i)], 'Color', tum_gray, 'LineWidth', 30, 'HandleVisibility', 'off'); hold on; % First Link
        plot([x1_pos_strap2(i), x2_pos_strap2(i)], [y1_pos_strap2(i), y2_pos_strap2(i)], 'Color', tum_gray, 'LineWidth', 30, 'HandleVisibility', 'off'); hold on; % First Link
        plot([x1_pos_strap3(i), x2_pos_strap3(i)], [y1_pos_strap3(i), y2_pos_strap3(i)], 'Color', tum_gray, 'LineWidth', 25, 'HandleVisibility', 'off'); hold on; % First Link
        plot([x1_pos_strap4(i), x2_pos_strap4(i)], [y1_pos_strap4(i), y2_pos_strap4(i)], 'Color', tum_gray, 'LineWidth', 25, 'HandleVisibility', 'off'); hold on; % First Link


        
        % Plot joints and end effector
        rectangle('Position', [-radius_joints, -radius_joints, 2*radius_joints, 2*radius_joints], 'Curvature', [1 1], ...
            'FaceColor', man_color, 'HandleVisibility', 'Off'); hold on;
%         rectangle('Position', [hx_pos_link0(i)-radius_joints, hy_pos_link0(i)-radius_joints, 2*radius_joints, 2*radius_joints], 'Curvature', [1 1], ...
%             'FaceColor', human_joint_color, 'HandleVisibility', 'Off'); hold on;
        rectangle('Position', [ l1*cos(x(i,1))-radius_joints, l1*sin(x(i,1))-radius_joints, 2*radius_joints, 2*radius_joints], ...
            'Curvature', [1 1], 'FaceColor', man_color, 'HandleVisibility', 'Off');
%         rectangle('Position', [hx_pos_link1(i)-radius_joints, hy_pos_link1(i)-radius_joints, 2*radius_joints, 2*radius_joints], 'Curvature', [1 1], ...
%             'FaceColor', human_joint_color, 'HandleVisibility', 'Off'); hold on;
        rectangle('Position', [hx_pos_hand(i)-radiushand,hy_pos_hand(i)-radiushand, ...
            2*radiushand, 2*radiushand], 'Curvature', [1 1], 'FaceColor', human_color, 'HandleVisibility', 'Off');
%         rectangle('Position', [manipulator_positions.x(i)-radius_man, manipulator_positions.y(i)-radius_man, ...
%            2*radius_man, 2*radius_man], 'Curvature', [1 1], 'FaceColor', tum_orange, 'HandleVisibility', 'Off');
                



%         plot([hx_pos_link2(i), hx_pos_hand(i)], ...
%           [hy_pos_link2(i), hy_pos_hand(i)],'-o', 'MarkerSize',3, 'MarkerEdgeColor',human_color,...
%           'MarkerFaceColor',human_color, 'Color', human_color, 'LineWidth', 18, 'HandleVisibility', 'off');


        plot(manipulator_positions.x(1:i), manipulator_positions.y(1:i), 'Color', tum_blue, 'LineWidth', 0.8);
        
        % Resize the plot
        xlim([-4/3*l, 4/3*l]); ylim([-4/3*l, 4/3*l]);
        
        if i <= increment
            labels = [labels; "Simulated trajectory"];
        end
        
        % Reference trajectory
        plot(l1*cos(x_ref.q1) + l2*cos(x_ref.q1+x_ref.q2), l1*sin(x_ref.q1) + l2*sin(x_ref.q1+x_ref.q2), ...
            'Color', tum_orange, 'LineStyle', '--', 'LineWidth', 0.3);
        
        if i <= increment
            labels = [labels; "Reference trajectory"];
        end
        
        % Plot motor angle indicators
        plot([0, radius_joints*cos(x(i,3))], [0, radius_joints*sin(x(i,3))], 'Color', 'k');
        plot([l1*cos(x(i,1)), l1*cos(x(i,1))+radius_joints*cos(x(i,7))], ...
            [l1*sin(x(i,1)), l1*sin(x(i,1))+radius_joints*sin(x(i,7))], 'Color', 'k', 'HandleVisibility', 'Off');
        
        if i <= increment
            labels = [labels; "Motor angles"];
        end
        
        % Constraint plotting
%         x_pos_constraint_link1 = 1/3*l1*cos(3/4*pi); y_pos_constraint_link1 = 1/3*l1*sin(3/4*pi);
%         
%         x_pos_constraint_link2_upper = l1*cos(x(i,1)) + 1/3*l2*cos(1/2*pi+x(i,1)); 
%         x_pos_constraint_link2_lower = l1*cos(x(i,1)) - 1/3*l2*cos(1/2*pi+x(i,1));
%         y_pos_constraint_link2_upper = l1*sin(x(i,1)) + 1/3*l2*sin(1/2*pi+x(i,1));
%         y_pos_constraint_link2_lower = l1*sin(x(i,1)) - 1/3*l2*sin(1/2*pi+x(i,1));
%         
%         plot([radius_joints*cos(3/4*pi), x_pos_constraint_link1], [radius_joints*sin(3/4*pi), y_pos_constraint_link1], 'k:', 'LineWidth', 1.5); 
%         plot([radius_joints*cos(3/4*pi), x_pos_constraint_link1], [-radius_joints*sin(3/4*pi), -y_pos_constraint_link1], 'k:', 'LineWidth', 1.5);
%         plot([l1*cos(x(i,1))+radius_joints*cos(pi/2+x(i,1)), x_pos_constraint_link2_upper], ...
%             [l1*sin(x(i,1))+radius_joints*sin(pi/2+x(i,1)), y_pos_constraint_link2_upper], 'k:', 'LineWidth', 1.5); 
%         plot([l1*cos(x(i,1))-radius_joints*cos(pi/2+x(i,1)), x_pos_constraint_link2_lower], ...
%             [l1*sin(x(i,1))-radius_joints*sin(pi/2+x(i,1)), y_pos_constraint_link2_lower], 'k:', 'LineWidth', 1.5);
%         
%         if i <= increment
%             labels = [labels; "Constraints"];
%         end
        
%         title("Animated motion of a 2-dof robot");
%         title('Planar robot with 2 degrees of freedom');

        xlabel("[m]"); ylabel("[m]");
        legend(labels,'Location','northwest');
        grid on;

        % Draw variable values
        txt1 = text(0.8, 0.95, "t="+t(i)+"s", 'Units', 'normalized');
        txt2 = text(0.8, 0.85, "q_1="+x(i,1)/pi*180+"°", 'Units', 'normalized'); hold off;
        txt3 = text(0.8, 0.80, "q_2="+x(i,5)/pi*180+"°", 'Units', 'normalized');
        format shortE
        txt4 = text(0.35, 0.95, "LoG-GP update time: " + sprintf('%.1e',plots.gptime_plot(i)) + "s", 'Units', 'normalized');
        format
        
        set(findall(gcf,'-property','FontSize'),'FontSize',18)
        
        % external torque estimate
        xt = t(1:increment:i);
        x2 = [xt', fliplr(xt')];
        conf_region1 = [plots.AKF_lbs(1:increment:i,9)', fliplr(plots.AKF_ubs(1:increment:i,9)')];
        conf_region2 = [plots.AKF_ubs(1:increment:i,10)', fliplr(plots.AKF_lbs(1:increment:i,10)')];
        
        subplot(2,5,[4, 5])
        title('Link1 external torque estimation');
        plot(xt,tau(1:increment:i,1),'Color',human_joint_color,'LineWidth',2); hold on;
        fill(x2, conf_region1,[0, 1, 1],'EdgeColor','none'); hold on;
        plot(xt,plots.tau_d_KF_hat_plot(1:increment:i,1),'b','LineWidth',3);hold on;
        plot(xt,tau_d(1:increment:i,1),'--','Color',tum_orange,'LineWidth',3); hold on;
        
        legend('Motor torque $\tau_{m,1}$','$2\sigma$ confidence region','Estimate $\hat{\tau}_{ext,1}$','Ground truth $\tau_{ext,1}$','Interpreter', 'latex');
        grid on;
        xlim([t(1),t(end)]);
        ylabel("Torque [Nm]"); xlabel("Time [s]");
        set(findall(gcf,'-property','FontSize'),'FontSize',18)
        
        subplot(2,5,[9, 10])
        title('Link2 external torque estimation');
        plot(xt,tau(1:increment:i,2),'Color',human_joint_color,'LineWidth',2); hold on;
        fill(x2, conf_region2,[0, 1, 1],'EdgeColor','none'); hold on;
        plot(xt,plots.tau_d_KF_hat_plot(1:increment:i,2),'b','LineWidth',3);hold on;       
        plot(xt,tau_d(1:increment:i,2),'--','Color',tum_orange,'LineWidth',3); hold on;
        legend('Motor torque $\tau_{m,1}$','$2\sigma$ confidence region','Estimate $\hat{\tau}_{ext,2}$','Ground truth $\tau_{ext,2}$','Interpreter', 'latex');
%         legend('AKF estimate $\hat{\tau}_{ext,2}$','$2\sigma$ ucb augmented state KF','$2\sigma$ lcb augmented state KF','true external torque $\tau_{ext,2}$','Interpreter', 'latex');
        grid on;
        xlim([t(1),t(end)]);
        ylabel("Torque [Nm]"); xlabel("Time [s]");
        set(findall(gcf,'-property','FontSize'),'FontSize',18)
        
        
        frame = getframe(gcf);
        writeVideo(video_handle,frame);
    end
    toc;
    
    close(video_handle);
end

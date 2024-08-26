function [fig_simulation,video_handle] = drawAnimatedSystem(t, x, x_ref, sys_params)
%
    % Load system parameters
    l1 = sys_params.l1; 
%     l2 = sys_params.l2;

    % Compute joint positions upfront
    x_pos_link1 = l1*cos(x(:,1)); 
    y_pos_link1 = l1*sin(x(:,1)); 
    
%     x_pos_link2 = x_pos_link1 + l2*cos(x(:,1)+x(:,5));
%     y_pos_link2 = y_pos_link1 + l2*sin(x(:,1)+x(:,5));

    manipulator_positions.x = x_pos_link1;
    manipulator_positions.y = y_pos_link1;
    
    video_handle = VideoWriter('animation');
    video_handle.FrameRate = 15;
    open(video_handle);
    
    fig_simulation = figure('Name', 'Simulation');
    
    % Initialize objects
    txt1 = text(0.9,0.9,"");  txt1.Units = "normalized";
    txt2 = text(0.9,0.8,"");  txt2.Units = "normalized";
    txt3 = text(0.9,0.7,"");  txt3.Units = "normalized";

    % Animate system
    tic;
    i = 0; increment = 1;
    while i + increment < length(t)   
        i = i + increment;
        
        % Draw robot at current position
        figure(fig_simulation); clf(fig_simulation);

        scatter(0, 0, 50, 'c', 'Filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5); hold on; % Base Point / First Joint
        scatter(x_pos_link1(i), y_pos_link1(i), 50, 'c', 'Filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5); % Second Joint
%         scatter(x_pos_link2(i), y_pos_link2(i), 50, 'r', 'Filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5); % Manipulator

        xlim([-l1, l1]); ylim([-l1, l1]);
        
        plot([0, x_pos_link1(i)], [0, y_pos_link1(i)], 'k', 'LineWidth', 1.5); % First Link
%         plot([x_pos_link1(i), x_pos_link2(i)], [y_pos_link1(i), y_pos_link2(i)], 'k', 'LineWidth', 1.5); % Second Link
        plot(manipulator_positions.x(1:i), manipulator_positions.y(1:i), 'g', 'LineWidth', 2);
        plot(l1*cos(x_ref.q1), l1*sin(x_ref.q1), 'r--');
        grid on;

        % Draw variable values
        txt1 = text(0.9, 0.95, "t="+t(i)+"s", 'Units', 'normalized');
        txt2 = text(0.9, 0.85, "q_1="+x(i,1)/pi*180+"°", 'Units', 'normalized'); hold off;
        txt3 = text(0.9, 0.80, "\theta_1="+x(i,3)/pi*180+"°", 'Units', 'normalized');
        
        frame = getframe(gcf);
        writeVideo(video_handle,frame);
    end
    toc;
    
    close(video_handle);
end
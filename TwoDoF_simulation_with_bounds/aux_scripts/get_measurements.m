function meas = get_measurements(system_states,current_x,i,sample_time)
qdd1meas = (current_x(2) - system_states.q1_dot(i-1))/sample_time;
qdd2meas = (current_x(6) - system_states.q2_dot(i-1))/sample_time;
meas.q = [current_x(1); current_x(5)];
meas.qd = [current_x(2); current_x(6)];
meas.qdd = [qdd1meas; qdd2meas];
meas.theta_m = [current_x(3); current_x(7)];
meas.thetad_m = [current_x(4); current_x(8)];
thetadd1meas = (current_x(4) - system_states.theta1_dot(i-1))/sample_time;
thetadd2meas = (current_x(8) - system_states.theta2_dot(i-1))/sample_time;
meas.thetadd_m = [thetadd1meas; thetadd2meas];
end


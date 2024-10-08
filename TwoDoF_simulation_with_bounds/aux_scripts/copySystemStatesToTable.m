function system_states_table = copySystemStatesToTable(system_states_table,i,current_t,current_x,current_tau)
%
    system_states_table.t(i) = current_t;
    
    system_states_table.q1(i)           = current_x(1);
    system_states_table.q1_dot(i)       = current_x(2);
    system_states_table.theta1(i)       = current_x(3);
    system_states_table.theta1_dot(i)   = current_x(4);
    system_states_table.q2(i)           = current_x(5);
    system_states_table.q2_dot(i)       = current_x(6);
    system_states_table.theta2(i)       = current_x(7);
    system_states_table.theta2_dot(i)   = current_x(8);
    
    system_states_table.tau1(i)     = current_tau(1);
    system_states_table.tau2(i)     = current_tau(2);
end
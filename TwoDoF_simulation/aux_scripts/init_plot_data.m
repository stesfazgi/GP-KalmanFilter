function plots = init_plot_data(t)
plots.tau_d_KF_hat_plot = zeros(size(t,1), 2);
plots.x_KF_plot = zeros(size(t,1), 8);

plots.AKF_ubs = zeros(size(t,1),10);
plots.AKF_lbs = zeros(size(t,1),10);
plots.AKF_cov_plot = zeros(size(t,1),10);

plots.gptime_plot = zeros(size(t,1), 1);

plots.q = zeros(size(t,1),2);
plots.qd = zeros(size(t,1),2);
plots.qdd = zeros(size(t,1),2);

plots.xtrain = zeros(size(t,1),6);
plots.ytrain = zeros(size(t,1),2);

plots.jacvec = zeros(size(t,1),6);
end


function [u_mu,u_sigma] = compute_remainder_overapproximations(X, l_mu, l_sigma,dofs)
    %% Compute the (hyper-)rectangle over-approximating the lagrangians of mu and sigma
% 
%     Parameters
%     ----------
%     X: n_s x n_s ndarray[float]
%         The shape matrix of the current state ellipsoid
%     l_mu: n x 0 numpy 1darray[float]
%         The lipschitz constants for the gradients of the predictive mean
%     l_sigma n x 0 numpy 1darray[float]
%         The lipschitz constans on the predictive variance
% 
%     Returns
%     -------
%     u_mu: n_s x 0 numpy 1darray[float]
%         The upper bound of the over-approximation of the mean lagrangian remainder
%     u_sigma: n_s x 0 numpy 1darray[float]
%         The upper bound of the over-approximation of the variance lagrangian remainder
    X = X(1:4*dofs,1:4*dofs);
    n_s = size(X,1);
    s = eye(n_s);
    b = s* s';
    qb = X* b;
    evals = eig(qb); 
    r_sqr = max(evals);
    
%     # This is equivalent to:
%     # q_inv = sLA.inv(X)
%     # evals,_,_ = sLA.eig(b,q_inv)
%     # however we prefer to avoid the inversion
%     # and breaking the symmetry of b and X

    u_mu = l_mu/2 * r_sqr;
    u_sigma = l_sigma * sqrt(r_sqr);

end


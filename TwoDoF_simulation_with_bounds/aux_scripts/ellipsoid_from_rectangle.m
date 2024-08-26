function X_ell = ellipsoid_from_rectangle(u_b)
    %% Compute ellipsoid covering box
% 
%     Given a box defined by
% 
%         B = [l_b[0],u_b[0]] x ... x [l_b[-1],u_b[-1]],
%     where l_b = -u_b (element-wise),
%     we compute the minimum enclosing ellipsoid in closed-form
%     as the solution to a linear least squares problem.
%     This can be either done by a diagonal shape matrix (axis-aligned)
%     or a rotated/shifted ellipsoid
% 
%     Method is described in:
%         [1] :
% 
%     TODO:   Choice of point is terrible as of now, since it contains linearly dependent
%             points which are not properly handled.
% 
%     Parameters
%     ----------
%         u_b: array_like, 1d
%             list of length n containing upper bounds of intervals defining box (see above)
%     Returns
%     -------
%         X_ell: np.ndarray[float, n_dim = 2], array of size n x n
%             Shape matrix of covering ellipsoid

%     p = length(u_b);
%     d = p * u_b .^ 2;
%     X_ell = diag(d);
dims = length(u_b);
x = zeros(dims);
X = zeros(dims,dims,dims);
for i = 1:dims
	X(i,i,i) = u_b(i)^2;
end
[~,X_ell] = sum_m_ellipsoids(x,X);
end


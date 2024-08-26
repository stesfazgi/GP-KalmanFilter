function [c_new,q_new] = sum_two_ellipsoids(c_1, X_1, c_2, X_2)
    %%  Sum of two ellipsoids
%     Computes the ellipsoidal overapproximation of the sum of two n-dimensional
%     ellipsoids.
% 
%     Parameters
%     ----------
%     c_1,c_2: n x 1 array
%         The centers of the ellipsoids to sum
%     X_1,X_2: n x n array
%         The shape matrices of the two ellipsoids
%     Returns
%     -------
%     p_new: n x 1 array
%         The center of the resulting ellipsoid
%     q_new: n x n array
%         The shape matrix of the resulting ellipsoid

p = sqrt(trace(X_1) / trace(X_2));
c_new = c_1 + c_2;
if isempty(nonzeros(X_1))
  q_new = X_2;
elseif isempty(nonzeros(X_2))
  q_new = X_1;
else
  q_new = (1 + (1. / p)) * X_1 + (1 + p) * X_2;
end
end


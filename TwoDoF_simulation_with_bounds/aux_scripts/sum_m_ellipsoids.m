function [c_new,X_new] = sum_m_ellipsoids(cm,Xm)
    %% Ellipsoidal overapproximation of sum of multiple ellipsoids
% 
%     Compute an ellipsoidal overapproximation of the sum
%     of n individual m-dimensional ellipsoids.
% 
%     Parameters
%     ----------
%     cm: n x m array[float]
%         The centers of the input ellipsoids
%     Xm: n x m x m array[float]
%         The shape matrices of the input ellipsoids
%     Returns
%     -------
%     c_new: n x 1 array
%         The center of the resulting ellipsoid
%     X_new: n x n array
%         The shape matrix of the resulting ellipsoid

m = size(cm,2);

sum_pi = 0;
sum_piQ = zeros(size(Xm(:,:,1)));
for i = 1:m
  p_i = sqrt(trace(Xm(:,:,i)));
  if p_i > 0
    sum_pi = sum_pi + p_i;
    sum_piQ = sum_piQ + (1./p_i).*Xm(:,:,i);
  end
end

c_new = sum(cm,2);
X_new = sum_pi*sum_piQ;
end


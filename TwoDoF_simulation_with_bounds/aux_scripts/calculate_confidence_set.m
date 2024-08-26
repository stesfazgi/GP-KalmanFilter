function [center,X_shape] = calculate_confidence_set(c_1,X_1,c_2,P_2,prob)
%CALCULATE_CONFIDENCE_SET Summary of this function goes here
%   Detailed explanation goes here
s = chi2inv(1-prob,length(c_1));
P_2_sigma = s.*P_2;
if isempty(nonzeros(X_1))
  center = c_1;
  X_shape = P_2_sigma;
else
  [center,X_shape] = sum_two_ellipsoids(c_1,X_1,c_2,P_2_sigma);
end
end


function train_loggp()
load loggptraindat.mat
Y = y_train_plot;
X = x_train_plot;

% ids = randperm(size(X,1));
% X = X(ids,:);
% Y = Y(ids,:);

% hyperparameter optimization
Npretrain = length(X);
disp('Hyperparameter training ...');
for dim = 1:size(Y,2)
  
  disp(['output dimension ' num2str(dim)]);
  gp = fitrgp(X(1:Npretrain-1,:),Y(1:Npretrain-1,dim),'KernelFunction','ardsquaredexponential','Standardize',false);
  
  
  ls = gp.KernelInformation.KernelParameters(1:end-1);
  sf = gp.KernelInformation.KernelParameters(end);
  sn = gp.Sigma;
  
  paramfile = ['hyperparam_dim_' num2str(dim) '.mat'];
  save(paramfile, 'ls', 'sf', 'sn');

end
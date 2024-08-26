function train_loggp()
load loggptraindat.mat
Y = trans_traindata;
X = x_traindata;

% hyperparameter optimization
Npretrain = length(X);
disp('Hyperparameter training ...');
for dim = 1:size(Y,2)

disp(['output dimension ' num2str(dim)]);
gp = fitrgp(X(1:Npretrain,:),Y(1:Npretrain,dim),'KernelFunction','ardsquaredexponential','Standardize',false);

ls = gp.KernelInformation.KernelParameters(1:end-1);
sf = gp.KernelInformation.KernelParameters(end);
sn = gp.Sigma;

paramfile = ['hyperparam_dim_' num2str(dim) '.mat'];
save(paramfile, 'ls', 'sf', 'sn');

end
function [phi, t_stat, var_res, resids, aic, LogV, std_coeff] = AR_LogV(p)
    %Fonction permettant d'estimer un modèle AR(p) par la méthode du maximum de
    %vraisemblance

global Yt;
warning('off');

%variable à expliquer
y = Yt;

%nombre d'observations
nobs = size(y,1) - p;

%estimation du modèle 
model = arima('ARLags',1:p);
[model, var_covar, LogV] = estimate(model, y, 'Display', 'off');

%AIC du modèle
aic = aicbic(LogV, p+1, nobs);

%Inverse de la Log Vraissemblance car estimate retourne une LogV négative
LogV = - LogV;

%résidus, SCR et variance résiduelle
resids = infer(model, y);
SCR = resids'*resids;
var_res = SCR/(nobs - (p+1));

%coefficients optimaux estimés 
phi = cell2mat([model.Constant; model.AR']);

%variables explicatives lagés 
X = lagmatrix(y, (1:p));
X = X(p+1:end,:);
X = [ones(size(X,1),1) X]; 

%matrice de variance-covariance du modèle
var_cov = var_res * inv(X'*X);

%ecart type des coefficients estimés 
std_coeff = sqrt(diag(var_cov));

%t-stat des coefficients
t_stat = phi ./ std_coeff;

end
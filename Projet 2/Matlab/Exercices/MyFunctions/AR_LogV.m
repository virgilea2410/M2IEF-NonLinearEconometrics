function [phi, t_stat, var_res, resids, aic, LogV, std_coeff] = AR_LogV(p)
    %Fonction permettant d'estimer un mod�le AR(p) par la m�thode du maximum de
    %vraisemblance

global Yt;
warning('off');

%variable � expliquer
y = Yt;

%nombre d'observations
nobs = size(y,1) - p;

%estimation du mod�le 
model = arima('ARLags',1:p);
[model, var_covar, LogV] = estimate(model, y, 'Display', 'off');

%AIC du mod�le
aic = aicbic(LogV, p+1, nobs);

%Inverse de la Log Vraissemblance car estimate retourne une LogV n�gative
LogV = - LogV;

%r�sidus, SCR et variance r�siduelle
resids = infer(model, y);
SCR = resids'*resids;
var_res = SCR/(nobs - (p+1));

%coefficients optimaux estim�s 
phi = cell2mat([model.Constant; model.AR']);

%variables explicatives lag�s 
X = lagmatrix(y, (1:p));
X = X(p+1:end,:);
X = [ones(size(X,1),1) X]; 

%matrice de variance-covariance du mod�le
var_cov = var_res * inv(X'*X);

%ecart type des coefficients estim�s 
std_coeff = sqrt(diag(var_cov));

%t-stat des coefficients
t_stat = phi ./ std_coeff;

end
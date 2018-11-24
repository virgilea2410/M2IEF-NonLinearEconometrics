function [phi, resid, lbstat, jbstat] = ARMA(y,p,q)

%%                      Estimation du modèle
%Estimation du modèle par max de vraissemblance
opt_modl = arima(p, 0, q);
[opt_est_mdl, opt_est_vcov, ~ , info] = estimate(opt_modl, y, 'print', false);

%Coefficients estimés du modèle
phi = cell2mat([opt_est_mdl.Constant opt_est_mdl.AR opt_est_mdl.MA]);

%résidus du modèle
resid = infer(opt_est_mdl, y);

%SCR du modèle
SCR = sum(resid'*resid);

%variance résiduelle
var_res = SCR/(size(y,1)-(1+p+q));

%ecart type résiduel
etyperes = sqrt(var_res);

%stat de student des coefficients
t_student = phi' ./ sqrt(diag(opt_est_vcov(1:3,1:3)));

%%                      Spécification du modèle

%autocorrélation des erreurs 
[lb, plb, lbstat, critikValLB] = lbqtest(resid);

%normalité des erreurs
[jb, pjb, jbstat, critikValJB] = jbtest(resid);

end

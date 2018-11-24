function [phi, resid, lbstat, jbstat] = ARMA(y,p,q)

%%                      Estimation du mod�le
%Estimation du mod�le par max de vraissemblance
opt_modl = arima(p, 0, q);
[opt_est_mdl, opt_est_vcov, ~ , info] = estimate(opt_modl, y, 'print', false);

%Coefficients estim�s du mod�le
phi = cell2mat([opt_est_mdl.Constant opt_est_mdl.AR opt_est_mdl.MA]);

%r�sidus du mod�le
resid = infer(opt_est_mdl, y);

%SCR du mod�le
SCR = sum(resid'*resid);

%variance r�siduelle
var_res = SCR/(size(y,1)-(1+p+q));

%ecart type r�siduel
etyperes = sqrt(var_res);

%stat de student des coefficients
t_student = phi' ./ sqrt(diag(opt_est_vcov(1:3,1:3)));

%%                      Sp�cification du mod�le

%autocorr�lation des erreurs 
[lb, plb, lbstat, critikValLB] = lbqtest(resid);

%normalit� des erreurs
[jb, pjb, jbstat, critikValJB] = jbtest(resid);

end

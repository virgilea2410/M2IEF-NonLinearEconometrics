function [phi, resids, var_res, var_phi, t_stat, pval] = SETAR(y,p,d,c)

T = size(y,1);
nobs = T - max(p,d);

%préparation de la matrice des regresseurs
X = zeros(nobs, 2*(p+1));

%préparation de la variable indicatrice
Q = lagmatrix(y,d);
I = Q<=c;

%Création de la matrice des explicatives x(t-1) ... x(t-p)
Xs = lagmatrix(y,(1:p));
Xs = Xs((p+1):end,:);
I = I((p+1):end,:);
Xs1 = Xs.*I;
Xs2 = Xs.*(1-I);

%Décalage du vecteur des variables expliquées
Y = y((p+1):end,:);

const1 = ones(size(Y,1),1).*I;
const2 = ones(size(Y,1),1).*(1-I);

%Matrice des régresseurs
X = [const1 Xs1 const2 Xs2];

%Regression
[phi, ~, resids, ~, stats] = regress(Y,X);

%Variance résiduelle
SCR = resids'*resids;
var_res = SCR/(nobs-2*(p+1));

var_phi = var_res * inv(X'*X);

t_stat = phi ./ sqrt(diag(var_phi));

pval = (1-tcdf(t_stat,nobs));

end
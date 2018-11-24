function [ phi, SCR, var_res, resids, var_phi, t_student ] = AR(y,p)

T = size(y,1);

nobs = (T-p);

%Vecteur des variables à expliquer
Y = y((p+1):end,:);

%construction des 3 vecteurs de variables explicatives : x(t-1) = Y1 ; x(t-2) = Y2 ; x(t-3) = Y3
Xs = lagmatrix(y,(1:p));
Xs = Xs((p+1):end,:);

%On construit ensuite la matrice de la régression, composé de 4 colonnes : une colonne de 1 pour la constante et trois colonnes pour les 3 varibales
%explicatives
X = [ones(size(Xs,1),1) Xs];

%on calcul les coefficients de la régression
[beta, ~, resids, ~, stats] = regress(Y,X);
phi = beta;

%ensuite, SCR = somme des résidus au carrés
SCR = resids'*resids;

%et l'ecart type est égale à SCR/nombre d'observations
var_res = SCR/(nobs-(p+1));

%variance des coefficients
var_phi = var_res * inv(X'*X);

%stats de student
t_student = phi ./ sqrt(diag(var_phi));

end



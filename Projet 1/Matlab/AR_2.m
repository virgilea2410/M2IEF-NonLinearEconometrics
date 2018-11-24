function [ phi, SCR, sigma ] = AR_2(z,y,p)

T = size(y,1);

nobs = (T-p);

%Vecteur des variables à expliquer
%Y = y((p+1):end,:);
Y = z((p+1):end,:);

%construction des 3 vecteurs de variables explicatives : x(t-1) = X1 ;
%x(t-2) = X2 ; x(t-3) = X3

X1 = y(1:(end-(p)),:);

X2 = y(2:(end-(p-1)),:);

X3 = y(3:(end-(p-2)),:);

%On construit ensuite la matrice de la régression, composé de 4 colonnes :
%une colonne de 1 pour la constante et trois colonnes pour les 3 varibales
%explicatives
X = [ ones(nobs,1) X1 X2 X3 ];

%on calcul les coefficients de la régression
beta = regress(Y,X);
phi = beta;

%on calcule le vecteur des variables à expliquer estimés - les variables à
%expliquer constaté
Y_estime = X*beta;
Y_constat = Y;


%on calcul les résidus : différence entre les xt constatés et les xt estimé
%
residu = Y_constat - Y_estime;

%ensuite, SCR = somme des résidus au carrés
SCR = sum(residu'*residu);

%et l'ecart type est égale à SCR/nombre d'observations
sigma = SCR/(nobs-(p+1));

end



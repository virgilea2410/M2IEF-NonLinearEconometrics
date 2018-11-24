function [LM_stat, seuil_critique, proba_critique, F_stat, seuil_critique_F, proba_critique_F] = TestSTAR(y,p,d)

Y= y(max(p,d) + 1:end, :);
T = size(y,1) - max(p,d);

rdt = y;

%% TEST DU MULTIPLICATEUR DE LAGRANGE 

%ESIMATION DE LA REGRESSION AUXILIAIRE SOUS H0 + CALCUL DES RESIDUS
[phi_AR, SCR_AR, sigma_AR, residus] = AR(rdt,p);

%REGRESSION DES RESIDUS SUR ZT
multipl = lagmatrix(rdt,1);
multipl2 = lagmatrix(rdt.^2, 1);
multipl3 = lagmatrix(rdt.^3, 1);
X1 = lagmatrix(rdt,(1:p));

X1 = X1(max(p,d) + 1:end,:);
multipl = multipl(max(p,d) + 1:end,:);
multipl2 = multipl2(max(p,d) + 1:end,:);
multipl3 = multipl3(max(p,d) + 1:end,:);

%matrice des regresseurs : p+1 + 3*p colonnes
X = [ones(T,1) X1 X1.*multipl X1.*multipl2 X1.*multipl3];

%on définit le Y comme le vecteur des résidus, afin d'effectuer la
%régression des résidus sur les x(t-1), x(t-2)
Y = residus((max(p,d) - min(p,d) + 1):end,:);

[beta, ~, ~, ~, r2]  = regress(Y,X);

% CALCUL DE LA STAT DU TEST 

r2 = r2(1,1);

LM_stat = T*r2;

seuil_critique = chi2inv(0.95, 3*p);

proba_critique = (1-chi2cdf(LM_stat, 3*p));

%% PARTIE 2 : TEST DE FISHER

%ESTIMATION DE LA REGRESSION AUXILIAIRE CONTRAINTE (SOUS HO, C'EST UN AR) + RESIDUS + SOMME DES CARRES RESIDUEL 

%Regression auxiliaire contraire = AR
[phi_AR2, SCR_AR2, sigma_AR2, residus2] = AR(rdt,p);

% QUESTION 2 : ESTIMATION DE LA REGRESSION AUXILIAIRE NON CONTRAINTE (PAS SOUS H0 DONC LES VECTEURS BETA NE SONT PAS EGAUX A 0) + RESIDUS + SOMME DES CARRES RESIDUEL

%on redéfinit le Y comme le vecteur des xt, afin d'effectuer la
%régression des xt sur les x(t-1), x(t-2)
Y = rdt(max(p,d) + 1:end, :);

[beta2, ~, residus_F]  = regress(Y,X);

SCR_F = sum(residus_F'*residus_F);

% CALCUL DE LA F STAT

F_stat = ((SCR_AR2 - SCR_F)/(3*p))/((SCR_F)/(T - 4*p - 1));

% CALCUL DU SEUIL CRITIQUE ET DE LA PROBA CRITIQUE

seuil_critique_F = finv(0.95, 3*p, T-4*p-1);
proba_critique_F = (1- fcdf(F_stat, 3*p, T-4*p-1));

end


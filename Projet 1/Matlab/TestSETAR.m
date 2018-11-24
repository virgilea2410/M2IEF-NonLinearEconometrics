function [decision, pvalue, statDuTest, Fstat] = TestSETAR(y,p,d,c)
%Retourne 1 si H0 est rejeté, 0 sinon.
% H0 : Le moidèle est linéaire
% H1 : le modèle est non linéaire

%%              Etape 1 : Calcul de la F-stat

[~, ~, std_AR] = AR(y,p);

[~, ~, std_SETAR] = SETAR(y,p,d,c);

T = size(y,1) - max(p,d);

statDuTest = T * (std_AR - std_SETAR)/std_SETAR;

%%            Etape 2 : Calcul de la distribution de la F-stat par simulation

%CALCUL DE LA DISTRIBUTION DE F(c_chapo) PAR SIMULATION 
J = 20;

z = ones(size(y,1),J);

%Génération de J vecteurs normale, puis estimation SETAR et AR de chacun
%des vecteurs
for i=1:J
    z(:,i) = normrnd(0,1,size(y,1), 1);
    [ARphi(:,i), ARresid(:,i), var_res_AR(i,1)] = AR_2(z(:,i),y,p);
    c(i,1) = SETAR_FIT_2(z(:,i), y, p,d);
    [ SETARphi(:,i), SCR_SETAR(:,i), var_res_SETAR(i,1)] = SETAR_2(z(:,i), y, p, d, c(i,1));
    Fstat(i,1) = T * (var_res_AR(i,1) - var_res_SETAR(i,1))/var_res_SETAR(i,1);
end

I = Fstat > statDuTest;
pvalue = sum(I)/J;

if pvalue <= 0.05
     decision = 1;
else
     decision = 0;
end

end

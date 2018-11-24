function [lb, pVal, stat_lb, critikVal] = MyLbTest(resids, k)
    %Fonction de Ljung-Box permettant de tester l'autocorrélation entre les
    %résidus d'un modèle.
    %retourne 1 si les résidus sont autocorrélés (rejet de 0)
    %retourne 0 sinon (non rejet de H0)

warning('off');

%Le vecteur res_autoc contient les autocorrélation entre les résidus pour
%un décalage allant de 1 jusqu'à 20
res_autoc = autocorr(resids);

%on supprime la première observation qui contient l'autocorrélation
%résiduelle d'ordre 0, qui vaut evidemment 1
res_autoc = res_autoc(2:end,1);

%H = 20
H = size(res_autoc,1);

%initialisation de la statistique du test
stat_lb = 0;

%taille de l'échantillon de résidus
n = size(resids,1);

%la statistique du test est la somme pour i allant de 1 à 20 des
%autocorrélations d'ordre i au carré, divisé par (n-i)...
for i=1:H
    stat_lb = stat_lb + res_autoc(i,1)^2/(n-i);
end

%... qu'on multiplie par n * (n-2)
stat_lb = n * (n+2) * stat_lb;

%calcul de la valeur critique
critikVal = chi2inv(0.95, k);

%décision du test
lb = stat_lb > critikVal;

%probabilité critique du test
pVal = 1- chi2cdf(stat_lb, k);

end
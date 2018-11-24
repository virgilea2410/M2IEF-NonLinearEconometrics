function [lb, pVal, stat_lb, critikVal] = MyLbTest(resids, k)
    %Fonction de Ljung-Box permettant de tester l'autocorr�lation entre les
    %r�sidus d'un mod�le.
    %retourne 1 si les r�sidus sont autocorr�l�s (rejet de 0)
    %retourne 0 sinon (non rejet de H0)

warning('off');

%Le vecteur res_autoc contient les autocorr�lation entre les r�sidus pour
%un d�calage allant de 1 jusqu'� 20
res_autoc = autocorr(resids);

%on supprime la premi�re observation qui contient l'autocorr�lation
%r�siduelle d'ordre 0, qui vaut evidemment 1
res_autoc = res_autoc(2:end,1);

%H = 20
H = size(res_autoc,1);

%initialisation de la statistique du test
stat_lb = 0;

%taille de l'�chantillon de r�sidus
n = size(resids,1);

%la statistique du test est la somme pour i allant de 1 � 20 des
%autocorr�lations d'ordre i au carr�, divis� par (n-i)...
for i=1:H
    stat_lb = stat_lb + res_autoc(i,1)^2/(n-i);
end

%... qu'on multiplie par n * (n-2)
stat_lb = n * (n+2) * stat_lb;

%calcul de la valeur critique
critikVal = chi2inv(0.95, k);

%d�cision du test
lb = stat_lb > critikVal;

%probabilit� critique du test
pVal = 1- chi2cdf(stat_lb, k);

end
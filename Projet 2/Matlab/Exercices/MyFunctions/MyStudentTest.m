function [decision, pvalue, critikVal] = MyStudentTest(x, alpha)
    %fonction retournant 1 si toutes les t-stat d'un vecteur sont
    %significatives, 0 sinon

global Yt;

warning('off');

%nombre de params estim�s, taille de la s�rie, nombre d'observations
k = size(x,1);
n = size(Yt,1);
nobs = n - k;

%calcul de la probabilit� critique, de la valeur critique
pvalue = 1 - tcdf(abs(x), nobs);
critikVal = tinv(1-alpha, nobs);

%vecteur de bool�an Z, 1 si la valeur absolue de la t-stat est > � la
%valeur critique, z�ron sinon
z = abs(x) > critikVal;

%on fais la somme des �lements du vecteur z, si la somme est �gal au nombre
%de param�tres estim�s, alors toutes les t-stat sont sup�rieures � la
%valeur critiques, et donc significatives
if sum(z) == k
    decision = 1;
elseif not(sum(z) == k)
    decision = z;
elseif sum(z) == 0
    decision = 0;
end

end
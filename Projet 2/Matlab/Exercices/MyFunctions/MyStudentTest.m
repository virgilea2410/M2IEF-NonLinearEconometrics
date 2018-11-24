function [decision, pvalue, critikVal] = MyStudentTest(x, alpha)
    %fonction retournant 1 si toutes les t-stat d'un vecteur sont
    %significatives, 0 sinon

global Yt;

warning('off');

%nombre de params estimés, taille de la série, nombre d'observations
k = size(x,1);
n = size(Yt,1);
nobs = n - k;

%calcul de la probabilité critique, de la valeur critique
pvalue = 1 - tcdf(abs(x), nobs);
critikVal = tinv(1-alpha, nobs);

%vecteur de booléan Z, 1 si la valeur absolue de la t-stat est > à la
%valeur critique, zéron sinon
z = abs(x) > critikVal;

%on fais la somme des élements du vecteur z, si la somme est égal au nombre
%de paramètres estimés, alors toutes les t-stat sont supérieures à la
%valeur critiques, et donc significatives
if sum(z) == k
    decision = 1;
elseif not(sum(z) == k)
    decision = z;
elseif sum(z) == 0
    decision = 0;
end

end
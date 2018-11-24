function [ c_optim, phi_optim ] = SETAR_FIT_2(z, x, p, d)

% on d�finit notre vecteur c comme �tant �gal au vecteur x des variables � expliquer (compos�e des x(t)) c = x;
c = z;

% on trie le vecteur c par ordre croissant
c = sort(c);

%on calcule le nombre d'observations de c, qui est �gale pour le moment au nombre d'observations du vecteur x
nobs = size(c,1);

%on prendre la partie enti�re de 15% du nombre d'observations, soit 15% * 560 = 84,xxx donc upDown = 84
upDown = floor(0.15 * nobs);

%on enleve les 84 premieres valeurs du vecteur c tri�
c = c(upDown:end,:);

%puis on enleve les 84 derni�res valeurs du vecteur c tri�
c = c(1:end - upDown,:);

% on d�finit un vecteur de 1, de la taille du vecteur c, afin qu'il contienne les SCR des regressions
result = ones(size(c,1),1);

% on fais une boucle de 1 jusqu'au nombre d'observations du vecteur c en appelant notre fonction SETAR, dans le but d'essayer chacune des valeurs possibles de C.
for i = 1:size(c,1)
    [ ~, result(i,1), ~] = SETAR_2(z, x, p, d, c(i,1));
end

% on cherche le scr minimum, ainsi que son index (son num�ro de ligne) dans le vecteur c
[~, index] = min(result);

%on affiche ensuite la valeur de c qui minimise le SCR
c_optim = c(index,1);

%enfin, on sort la matrice des coefficient optimaux, calcul�s a partir du c optimal question 3 du TD
phi_optim = ones(8,1);
[ phi_optim, ~, ~] = SETAR_2(z, x, p, d, c_optim);

%on calcule le SCR optimal
[ ~, SCR_optim, ~ ] = SETAR_2(z, x, p, d, c_optim);

%afin d'avoir la variance des r�sidus optimale (var des r�sidus = SCR/T) question 4 du TD
varRes_optim = SCR_optim/nobs;

% enfin, on calcule la matrice des variances convariances optimale question 5 du TD 
varcovar = varRes_optim .* inv(x'*x);
    
end

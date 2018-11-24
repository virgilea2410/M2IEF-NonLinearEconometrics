function [x, std_coeff, t_stat, fval] = maxVraissemblance_AR(x0)
    %Fonction permettant de maximiser, à l'aide d'un algorithme d'optimisation
    %numérique, la fonction de vraissemblance d'un modèle AR

global choice;
warning('off');

% Options de l'Optimisation
options = optimset('MaxIter', 10000, 'TolFun', 10^(-10), 'TolX', 10^(-10), 'Display', 'off');

% Appel de l'algorithme d'optimisation numérique
[x, fval, code, info, g, H] = fminunc(@vraissemblanceAR,x0, options);

%calcul de la valeur de l'ecart type des coefficients ainsi que de la t-stat
std_coeff = sqrt(diag(inv(H)));
t_stat = x./std_coeff;

end

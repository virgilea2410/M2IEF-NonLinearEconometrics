function [p_optim] = OptimLagByAutoc()
    %Fonction permettant d'optimiser le nombre de retard p* optimal pour les
    %variables explicatives du modèle AR sur Yt, par la méthode de
    %l'autocorrélation partielle

global Yt;
warning('off');

partial_autocorr = parcorr(Yt);

%on initialise le p* à 0
p_optim = 0;

%on itère p* tant que l'autocorrélation partielle d'ordre p est supérieure
%à 0,1
%on sélectionne p* tel que les autocorrélation partielles d'ordre 1->p soit
%toutes supérieures à 0,1
i = 1;
while partial_autocorr(i,1) > 0.1
    p_optim = p_optim + 1;
    i = i+1;
end

end

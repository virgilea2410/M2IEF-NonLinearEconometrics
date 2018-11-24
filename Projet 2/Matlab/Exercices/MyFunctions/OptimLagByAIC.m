function [p_optim] = OptimLagByAIC()
    %Fonction permettant d'optimiser le nombre de retard p* optimal pour les
    %variables explicatives du modèle AR sur Yt, par la méthode des critères
    %d'informations AIC

warning('off');

global Yt;

AR_aic = zeros(5,1);

%estimation des modèles AR(p) pour p allant de 1 jusqu'à 5
for i=1:5
    [~, ~, ~, ~, AR_aic(i,1)] = AR_LogV(i);
end

%séléction du p* optimal qui donne le critère AIC minimum
[~, index] = min(AR_aic);
p_optim = index;
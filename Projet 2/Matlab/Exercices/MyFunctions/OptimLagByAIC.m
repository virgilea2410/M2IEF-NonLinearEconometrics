function [p_optim] = OptimLagByAIC()
    %Fonction permettant d'optimiser le nombre de retard p* optimal pour les
    %variables explicatives du mod�le AR sur Yt, par la m�thode des crit�res
    %d'informations AIC

warning('off');

global Yt;

AR_aic = zeros(5,1);

%estimation des mod�les AR(p) pour p allant de 1 jusqu'� 5
for i=1:5
    [~, ~, ~, ~, AR_aic(i,1)] = AR_LogV(i);
end

%s�l�ction du p* optimal qui donne le crit�re AIC minimum
[~, index] = min(AR_aic);
p_optim = index;
function [p_optim] = OptimLagByAutoc()
    %Fonction permettant d'optimiser le nombre de retard p* optimal pour les
    %variables explicatives du mod�le AR sur Yt, par la m�thode de
    %l'autocorr�lation partielle

global Yt;
warning('off');

partial_autocorr = parcorr(Yt);

%on initialise le p* � 0
p_optim = 0;

%on it�re p* tant que l'autocorr�lation partielle d'ordre p est sup�rieure
%� 0,1
%on s�lectionne p* tel que les autocorr�lation partielles d'ordre 1->p soit
%toutes sup�rieures � 0,1
i = 1;
while partial_autocorr(i,1) > 0.1
    p_optim = p_optim + 1;
    i = i+1;
end

end

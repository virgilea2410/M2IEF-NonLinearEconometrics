function [MSAR_theta, MSAR_std_coeff, MSAR_t_stat, MSAR_logvraiss] = maxVraissemblance_2()
    %Fonction permettant de maximiser, à l'aide d'un algorithme d'optimisation
    %numérique, la fonction de vraissemblance d'un modèle MS-AR, et par la même
    %d'en calculer les probabilités filtrées

warning('off');

%initialisation des paramètres de l'optimisation
done = "False";
x0 = [0 0 0 0 0 0 0 0]';

%tant que la fonction maxVraissemblance n'est pas défini pour le vecteur de
%paramètres intiaux ...
while done == "False"
    try
        done = "True";
        [MSAR_theta, MSAR_std_coeff, MSAR_t_stat, MSAR_logvraiss] = maxVraissemblance(x0);
    catch
        %...on fais varier les probabilités de transition, et on itère sur le
        %vecteur des paramètres intiaux
        done = "False";
        x0(1:2,1) = normcdf(x0(1:2,1) + normrnd(0,1,2,1)); 
        x0(3:end,1) = x0(3:end,1) + 0.1;
    end
end

%/!\ Pour vérifier que les coefficients que nous avons trouvé via max de
%vraissemblance ne sont pas du à un maximum local pour la vraissemblance, 
%nous allons calculer, si possible, toute les vraissemblances pour un
%certain interval du vecteur de départ
% j=1;
% for i=1:15
%     try
%         [~, ~, ~, MSAR_logvraiss_3(:,j)] = maxVraissemblance(x0);
%     catch
%         MSAR_logvraiss_3(:,j) = 0;
%     end
%     j = j + 1;
%     x0 = x0 + 0.1;
% end

end

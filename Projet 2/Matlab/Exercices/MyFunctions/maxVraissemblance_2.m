function [MSAR_theta, MSAR_std_coeff, MSAR_t_stat, MSAR_logvraiss] = maxVraissemblance_2()
    %Fonction permettant de maximiser, � l'aide d'un algorithme d'optimisation
    %num�rique, la fonction de vraissemblance d'un mod�le MS-AR, et par la m�me
    %d'en calculer les probabilit�s filtr�es

warning('off');

%initialisation des param�tres de l'optimisation
done = "False";
x0 = [0 0 0 0 0 0 0 0]';

%tant que la fonction maxVraissemblance n'est pas d�fini pour le vecteur de
%param�tres intiaux ...
while done == "False"
    try
        done = "True";
        [MSAR_theta, MSAR_std_coeff, MSAR_t_stat, MSAR_logvraiss] = maxVraissemblance(x0);
    catch
        %...on fais varier les probabilit�s de transition, et on it�re sur le
        %vecteur des param�tres intiaux
        done = "False";
        x0(1:2,1) = normcdf(x0(1:2,1) + normrnd(0,1,2,1)); 
        x0(3:end,1) = x0(3:end,1) + 0.1;
    end
end

%/!\ Pour v�rifier que les coefficients que nous avons trouv� via max de
%vraissemblance ne sont pas du � un maximum local pour la vraissemblance, 
%nous allons calculer, si possible, toute les vraissemblances pour un
%certain interval du vecteur de d�part
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

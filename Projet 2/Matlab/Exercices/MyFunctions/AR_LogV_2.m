function [AR_phi, AR_t_stat, AR_var_res, AR_SCR, AR_aic, AR_logvraiss, std_coeff] = AR_LogV_2(p_optim)
    %Fonction permettant d'estimer un modèle AR(p_optim) par la méthode du maximum de
    %vraisemblance

global choice
global Yt;
warning('off');

%en fonction du retard p, définition du choix d'optimisation pour maxVraiss
if p_optim == 1
    choice = 1;
    x0_AR_MV = [0;0;0];
elseif p_optim == 2
    choice = 2;
    x0_AR_MV = [0;0;0;0];
elseif choice == 4
    choice = 4;
else
    disp("le retard p du modèle AR(p) passé en paramètre n est pas géré par la foction" ...
         + "\nVeuillez reessayer");
end

done = "False";

%taille de la variable à expliquer et nombre d'observations
n = size(Yt,1);
nobs = n - p_optim;

%tant que maxVraissemblance n'arrive pas à maximiser la fonction vraiss...
while done == "False"
    try
        done = "True";
        [AR_theta, std_coeff, AR_t_stat , AR_logvraiss] = maxVraissemblance_AR(x0_AR_MV);
    catch
        ... on itère le vecteur des paramètres initiaux, jusqu'a ce que la fonction soit correctement défini 
        done = "False";
        x0_AR_MV = x0_AR_MV + 0.1;
    end
end

%vecteurs des coefficients estimés et variance résiduelle du modèle 
if choice == 1
    AR_phi = AR_theta(1:2,1);
    AR_var_res = AR_theta(3,1);
elseif choice == 2
    AR_phi = AR_theta(1:3,1);
    AR_var_res = AR_theta(4,1);
end

%SCR du modèle
AR_SCR = AR_var_res * (n - p_optim - 1);

%nombre de paramètres estimés
AR_k = p_optim + 1;

%AIC du modèle
AR_aic = -2 * (-AR_logvraiss/nobs) + (2*(AR_k))/nobs;

end
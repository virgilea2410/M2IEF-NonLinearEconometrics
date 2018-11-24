function [AR_phi, AR_t_stat, AR_var_res, AR_SCR, AR_aic, AR_logvraiss, std_coeff] = AR_LogV_2(p_optim)
    %Fonction permettant d'estimer un mod�le AR(p_optim) par la m�thode du maximum de
    %vraisemblance

global choice
global Yt;
warning('off');

%en fonction du retard p, d�finition du choix d'optimisation pour maxVraiss
if p_optim == 1
    choice = 1;
    x0_AR_MV = [0;0;0];
elseif p_optim == 2
    choice = 2;
    x0_AR_MV = [0;0;0;0];
elseif choice == 4
    choice = 4;
else
    disp("le retard p du mod�le AR(p) pass� en param�tre n est pas g�r� par la foction" ...
         + "\nVeuillez reessayer");
end

done = "False";

%taille de la variable � expliquer et nombre d'observations
n = size(Yt,1);
nobs = n - p_optim;

%tant que maxVraissemblance n'arrive pas � maximiser la fonction vraiss...
while done == "False"
    try
        done = "True";
        [AR_theta, std_coeff, AR_t_stat , AR_logvraiss] = maxVraissemblance_AR(x0_AR_MV);
    catch
        ... on it�re le vecteur des param�tres initiaux, jusqu'a ce que la fonction soit correctement d�fini 
        done = "False";
        x0_AR_MV = x0_AR_MV + 0.1;
    end
end

%vecteurs des coefficients estim�s et variance r�siduelle du mod�le 
if choice == 1
    AR_phi = AR_theta(1:2,1);
    AR_var_res = AR_theta(3,1);
elseif choice == 2
    AR_phi = AR_theta(1:3,1);
    AR_var_res = AR_theta(4,1);
end

%SCR du mod�le
AR_SCR = AR_var_res * (n - p_optim - 1);

%nombre de param�tres estim�s
AR_k = p_optim + 1;

%AIC du mod�le
AR_aic = -2 * (-AR_logvraiss/nobs) + (2*(AR_k))/nobs;

end
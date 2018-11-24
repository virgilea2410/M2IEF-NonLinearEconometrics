function [LogV]= vraissemblance(x0)
    %fonction de vraisemblance d'un modèle MS-AR, calculant, puisque
    %nécessaire au calcul de la vraisemblance, également les probabilités
    %filtrées du modèle
    %
    %choice = 1 --> Y(t) = C(St) + phi1(St) + sigma_res(St)
    %choice = 2 --> Y(t) = C(St) + phi1 + phi2 + sigma_res(St)
    %choice = 3 --> Y(t) = C(St) + phi1 + phi2 + phi3 + phi4 + sigma_res

global choice;
global Yt;
warning('off');

p11= normcdf(x0(1,1));
p22= normcdf(x0(2,1));

%définition des paramètres initiaux
if choice == 1
    c_r1= x0(3,1);
    c_r2 =x0(4,1);
    phi1_r1 = x0(5,1);
    phi1_r2 = x0(6,1);
    sigma_res_r1 = x0(7,1);
    sigma_res_r2 = x0(8,1);
    sigma_res = [sigma_res_r1 sigma_res_r2]';
elseif choice == 2
    c_r1= x0(3,1);
    c_r2 =x0(4,1);
    phi1 = x0(5,1);
    phi2 = x0(6,1);
    sigma_res_r1 = x0(7,1);
    sigma_res_r2 = x0(8,1);
    sigma_res = [sigma_res_r1 sigma_res_r2]';
elseif choice == 3
    c_r1= x0(3,1);
    c_r2 =x0(4,1);
    phi1 = x0(5,1);
    phi2 = x0(6,1);
    phi3 = x0(7,1);
    phi4 = x0(8,1);
    sigma_res_r1 = x0(9,1);
    sigma_res = sigma_res_r1;
end

%vecteur et matrice des proba de transitions
PR_TR= [p11 1-p22; 1-p11 p22];
PR_TRF = [p11; 1-p11; 1-p22; p22];

% calcul des probabilites ergodiques pour initialiser le filtre
pi1 = (1-p22)/(2-p11-p22);   % P(St=1)
pi2 = 1-pi1;   % P(St=2)
PROB__T = [pi1; pi2];
PROB__ = [pi1; pi1; pi2; pi2]; 

%initialisation des paramètres calculés dans la boucle 
T = size(Yt,1);
LogV=0.0; 
PR_STL1 = zeros(T);  % init P(St|It-1)
PR_STT1 = zeros(T);  % init P(St|It)

if choice == 1
    i = 2;
elseif choice == 2
    i = 3;
elseif choice == 3
    i = 5;
end

for J_Iter =i:T  

    PROB_DD = PR_TRF .* PROB__; % P[St=i,St-1=j|It-1] (étape 1 algorithme filtrage)

    PROB_DD = PROB_DD(1:2,1)+ PROB_DD(3:4,1);   % P[St=i|It-1] (étape 2 algorithme filtrage)

    PR_STL1(J_Iter,1)=PROB_DD(1,1);  % P(St=1|It-1) 

    %calcul des résidus pour calculer la fonction de densité
    if choice == 1
        F_CAST = [Yt(J_Iter,1) - c_r1 - phi1_r1*Yt(J_Iter -1,1) ; Yt(J_Iter,1) - c_r2 - phi1_r2*Yt(J_Iter -1,1)];
    elseif choice == 2
        F_CAST = [Yt(J_Iter,1) - c_r1 - phi1*Yt(J_Iter -1,1) - phi2*Yt(J_Iter -2,1) ; Yt(J_Iter,1) - c_r2 - phi1*Yt(J_Iter -1,1)- phi2*Yt(J_Iter -2,1)];
    elseif choice == 3
        F_CAST = [Yt(J_Iter,1) - c_r1 - phi1*Yt(J_Iter -1,1) - phi2*Yt(J_Iter -2,1) - phi3*Yt(J_Iter -3,1) - phi4*Yt(J_Iter -4,1) ; Yt(J_Iter,1) - c_r2 - phi1*Yt(J_Iter -1,1)- phi2*Yt(J_Iter -2,1) - phi3*Yt(J_Iter -3,1) - phi4*Yt(J_Iter -4,1)];
    end

    PR_VL= (1./sqrt(2*pi*sigma_res)) .* exp(-0.5*F_CAST.^2./sigma_res) .* PROB_DD;  % f(yt,St|It-1) (étape 3 algorithme filtrage)
    %PR_VL= (1./sqrt(2*pi*(var_res.^2))) .* exp(-0.5*F_CAST.^2./(var_res.^2)) .* PROB_DD;  % f(yt,St|It-1) (étape 3 algorithme filtrage)

    %sommes des 2 fonctions de densités conditionnelles aux 2 états St
    PR_VAL= sum(PR_VL); % f(yt|It-1) (étape 4 algorithme filtrage)

    %calcul de la log-densité pour l'itération sur Y en cours, afin de
    %calculer la log-vraisemblance : produit des densités, somme des
    %log-densités
    %on prend l'inverse de la log-densité car on vas minimiser cette fonction
    LIK=-1*log(PR_VAL); 

    PROB__T= PR_VL ./ PR_VAL;   % P[St|It] (étape 5 algorithme filtrage)

    PR_STT1(J_Iter,1)=PROB__T(1,1);   % P(St=1|It) 

    PROB__ = [PROB__T(1,1);PROB__T(1,1);PROB__T(2,1);PROB__T(2,1)];  
    % nouveau PROB__ pour l'itération suivante

    %somme des log-densités pour chaque valeur de Y
    LogV = LogV + LIK;
end
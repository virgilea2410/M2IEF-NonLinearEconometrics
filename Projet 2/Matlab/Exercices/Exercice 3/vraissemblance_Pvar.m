function LIKV = vraissemblance_Pvar(x_v)
    %fonction de vraisemblance d'un modèle MS-AR à probabilités de
    %transition variables
    %Calcule également, puisque nécessaire au calcul de la vraisemblance,
    %les probabilités filtrées du modèle

global endog;
global exog_trans;

% definition des parametres 
c=x_v(3:4);
phi=x_v(5);
%VAR_L=[x_v(7)^2;x_v(8)^2];
VAR_L=[x_v(6);x_v(7)];

% une matrice de probabilites de transition initiale 
p11=normcdf(x_v(1));
p22=normcdf(x_v(2));

% calcul des probabilites ergodiques pour initialiser le filtre
pi1=(1-p22)/(2-p11-p22);
pi2=(1-p11)/(2-p11-p22);
PROB__= [pi1;pi1;pi2;pi2];  

T = size(endog,1);
LIKV=0.0; PR_STL1 = zeros(T,1); PR_STT1 = zeros(T,1);

for J_Iter=2:T  

    p11= normcdf(x_v(1) + x_v(end-1) * exog_trans(J_Iter));       
    p22= normcdf(x_v(2) + x_v(end) * exog_trans(J_Iter)); 

    PR_TRF=[p11;(1-p11);(1-p22);p22];

    PROB_DD=PR_TRF .* PROB__;  % Pr[St=i,St-1=j|It-1] (étape 1)

    PROB_DD = PROB_DD(1:2,1) + PROB_DD(3:4,1) ;   % Pr[St=i|It-1] (étape 2)

    PR_STL1(J_Iter,1)=PROB_DD(1,1);  % Pr(St=1|It-1) 

    F_CAST=endog(J_Iter,1)- c - phi*endog(J_Iter-1,1) ;  % calcul des residus

    PR_VL=(1./sqrt(2.*pi.*VAR_L)).*exp(-0.5*F_CAST.*F_CAST./VAR_L).*PROB_DD;
    % f(yt,St|It-1) (étape 3)

    PR_VAL=sum(PR_VL); 
    % f(yt|It-1), densite de yt sachant l'info passee (étape 4 du filtre)

    LIK=-1*log(PR_VAL); 

    PROB__T=PR_VL/PR_VAL; 
    % Pr[St|It] (étape 5)

    PR_STT1(J_Iter,1)=PROB__T(1,1);   % Pr(St=1|It) 

    PROB__ = [PROB__T(1,1);PROB__T(1,1);PROB__T(2,1);PROB__T(2,1)];  

    LIKV = LIKV + LIK;

end

end
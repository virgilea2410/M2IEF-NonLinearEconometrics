function LIKV=vraissemblance_PFixe(x)
    %fonction de vraisemblance d'un modèle MS-AR à probabilités de transition fixes
    %Calcule également, puisque nécessaire au calcul de la vraisemblance,
    %les probabilités filtrées du modèle

global endog

% definition des parametres 
p11=normcdf(x(1));
p22=normcdf(x(2));
c=x(3:4);
phi=x(5);
VAR_L=[x(6)^2;x(7)^2];

% une matrice de probabilites de transition initiale 
PR_TR=[p11 (1-p22);(1-p11) p22];
PR_TRF = PR_TR(:);

% calcul des probabilites ergodiques pour initialiser le filtre
pi1=(1-p22)/(2-p11-p22); pi2=(1-p11)/(2-p11-p22);
PROB__= [pi1;pi1;pi2;pi2];  

T = size(endog,1);
LIKV=0.0; PR_STL1 = zeros(T,1); PR_STT1 = zeros(T,1);

for J_Iter=2:T  

    PROB_DD=PR_TRF .* PROB__;  % Pr[St=i,St-1=j|It-1] 4x1  (étape 1 du filtre)

    PROB_DD = PROB_DD(1:2,1) + PROB_DD(3:4,1) ;   % Pr[St=i|It-1] 2x1  (étape 2 du filtre)

    PR_STL1(J_Iter,1)=PROB_DD(1,1);  % Pr(St=1|It-1) 

    F_CAST=endog(J_Iter,1)- c - phi*endog(J_Iter-1,1) ;  % calcul des residus

    PR_VL=(1./sqrt(2.*pi.*VAR_L)).*exp(-0.5*F_CAST.*F_CAST./VAR_L).*PROB_DD;
    % f(yt,St|It-1)  2x1 (étape 3 du filtre)

    PR_VAL=sum(PR_VL); 
    % f(yt|It-1), densite de yt sachant l'info passee (étape 4 du filtre)

    LIK=-1*log(PR_VAL); 
    % calcul de l'oppose de la log-vrais car fminunc mini les fonctions

    PROB__T=PR_VL/PR_VAL; 
    % Pr[St|It] 2x1 (étape 5 du filtre) 

    PR_STT1(J_Iter,1)=PROB__T(1,1);   % Pr(St=1|It) 

    PROB__ = [PROB__T(1,1);PROB__T(1,1);PROB__T(2,1);PROB__T(2,1)];  
    % nv PROB__ utilise ds l'iteration suivante 4x1

    LIKV = LIKV + LIK;

end
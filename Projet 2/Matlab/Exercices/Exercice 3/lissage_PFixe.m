function [PL,PF] = lissage_PFixe(x)
    %Fonction permettant de calculer les probabilités filtrées et lissées du
    %modèle MS-AR sur endog, à probabilités de transition fixes
    %PL = Proba Lissée = P[St=1 | I(T)]
    %PL = Proba Filtrée = P[St=1 | I(t)]
    %PP = Etape 2 du filtre = P[St=1 | I(t-1)]

global endog

p11=normcdf(x(1));      
[PF,PP]=filtrage(x);

T = size(endog,1);

PL=zeros(T,1);
PL(T,1)=PF(T,1);
for t=T-1:-1:1
    PLi1=PF(t,1)*p11/PP(t+1,1)*PL(t+1,1);
    PLi2=PF(t,1)*(1-p11)/(1-PP(t+1,1))*(1-PL(t+1,1));
    PL(t,1)=PLi1+PLi2;
end
end

function [PR_STT1,PR_STL1]=filtrage(x)

global endog

% definition des parametres 
p11=normcdf(x(1));    p22=normcdf(x(2));
c=x(3:4);
phi=x(5);
VAR_L=[x(6)^2;x(7)^2];

% une matrice de probabilites de transition initiale 
PR_TR=[p11 (1-p22);(1-p11) p22];
PR_TRF = PR_TR(:);  % ou PR_TRF=[p11;(1-p11);(1-p22);p22];

% calcul des probabilites ergodiques pour initialiser le filtre
pi1=(1-p22)/(2-p11-p22); pi2=(1-p11)/(2-p11-p22);
PROB__= [pi1;pi1;pi2;pi2];  

T = size(endog,1);
PR_STL1 = zeros(T,1); PR_STT1 = zeros(T,1);

for J_Iter=2:T  

    PROB_DD=PR_TRF .* PROB__;  
    % Pr[St=i,St-1=j|It-1] 4x1  (étape 1 du filtre)

    PROB_DD = PROB_DD(1:2,1) + PROB_DD(3:4,1) ;
    % Pr[St=i|It-1] 2x1  (étape 2 du filtre)

    PR_STL1(J_Iter,1)=PROB_DD(1,1);  % Pr(St=0|It-1) 

    F_CAST=endog(J_Iter,1)- c - phi*endog(J_Iter-1,1) ;  % calcul des residus (2x1)

    PR_VL=(1./sqrt(2.*pi.*VAR_L)).*exp(-0.5*F_CAST.*F_CAST./VAR_L).*PROB_DD;
    % f(yt,St|It-1)  2x1 (étape 3 du filtre)

    PR_VAL=sum(PR_VL); 
    % f(yt|It-1), densite de yt sachant l'info passee (étape 4 du filtre)

    PROB__T=PR_VL/PR_VAL; % Pr[St,St-1|It] 4x1
    % Pr[St|It] 2x1 (étape 5 du filtre) 

    PR_STT1(J_Iter,1)=PROB__T(1,1);   % Pr(St=0|It) 

    PROB__ = [PROB__T(1,1);PROB__T(1,1);PROB__T(2,1);PROB__T(2,1)];  
    % nv PRO__ utilise ds l'iteration suivante 4x1

end
end


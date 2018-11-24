function [PL,PF] = lissage_Pvar (x_v)
    %Fonction permettant de calculer les probabilités filtrées et lissées du
    %modèle MS-AR sur endog, à probabilités de transition variables, qui
    %sont fonction de exog_trans
    %PL = Proba Lissée = P[St=1 | I(T)]
    %PL = Proba Filtrée = P[St=1 | I(t)]
    %PP = Etape 2 du filtre = P[St=1 | I(t-1)]

global endog;
global exog_trans;

%p11 = normcdf[alpha0 + alpha1 * Z(t)]
p11= normcdf(x_v(1) + x_v(end-1) * exog_trans); 

[PF,PP]=filtrage_v(x_v);

T = size(endog,1);

PL=zeros(T,1); PL(T,1)=PF(T,1);
for t=T-1:-1:1
    PLi1= (p11(t+1,1) * PF(t,1)) / (PP(t+1,1)) * PL(t+1,1);
    PLi2= ((1-p11(t+1,1)) * PF(t,1)) / (1 - PP(t+1,1)) * (1 - PL(t+1,1));

    PL(t,1)=PLi1+PLi2;
end

function [PR_STT1,PR_STL1]=filtrage_v(x_v)

global endog;
global exog_trans;

% definition des parametres 
c=x_v(3:4);
phi=x_v(5);
VAR_L=[x_v(6)^2;x_v(7)^2];

% une matrice de probabilites de transition initiale 
p11=normcdf(x_v(1));    p22=normcdf(x_v(2));
PR_TR=[p11 (1-p22);(1-p11) p22];

% calcul des probabilites ergodiques pour initialiser le filtre
A = [eye(2)-PR_TR;ones(1,2)];
EN=[0;0;1];
PROB__T = (A'*A)\A'*EN;  % proba ergodiques Pr(St=1),Pr(St=2) 

PROB__= [PROB__T(1,1);PROB__T(1,1);PROB__T(2,1);PROB__T(2,1)];  

T = size(endog,1);
PR_STL1 = zeros(T,1); PR_STT1 = zeros(T,1);

for J_Iter=2:T  
    
    %p11 = normcdf[alpha0 + alpha1 * Z(J_Iter)]
    p11=normcdf(x_v(1)+x_v(end-1)*exog_trans(J_Iter));
    
    %p22 = normcdf[beta0 + beta1 * Z(J_Iter)]
    p22=normcdf(x_v(2)+x_v(end)*exog_trans(J_Iter));
    PR_TRF=[p11;(1-p11);(1-p22);p22];

    PROB_DD=PR_TRF .* PROB__;  
    % Pr[St=i,St-1=j|It-1] (étape 1)

    PROB_DD = PROB_DD(1:2,1) + PROB_DD(3:4,1) ;
    % Pr[St=i|It-1] (étape 2)

    PR_STL1(J_Iter,1)=PROB_DD(1,1);  % Pr(St=0|It-1) 

    F_CAST=endog(J_Iter,1)- c - phi*endog(J_Iter-1,1) ;  % calcul des residus (2x1)

    PR_VL=(1./sqrt(2.*pi.*VAR_L)).*exp(-0.5*F_CAST.*F_CAST./VAR_L).*PROB_DD;
    % f(yt,St|It-1) (étape 3)

    PR_VAL=sum(PR_VL); 
    % f(yt|It-1), densite de yt (étape 4)

    PROB__T=PR_VL/PR_VAL; % Pr[St,St-1|It
    % Pr[St|It] (étape 5) 

    PR_STT1(J_Iter,1)=PROB__T(1,1);   % Pr(St=0|It) 

    PROB__ = [PROB__T(1,1);PROB__T(1,1);PROB__T(2,1);PROB__T(2,1)];  
    % nv PRO__ utilise ds l'iteration suivante

end

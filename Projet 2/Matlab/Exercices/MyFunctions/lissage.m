function [PL,PF, PP] = lissage(x)
    %Fonction permettant de calculé les probabilités filtrées et lissées du modèle 
    %PL = Proba Lissée = P[St=1 | I(T)]
    %PL = Proba Filtrée = P[St=1 | I(t)]
    %PP = Etape 2 du filtre = P[St=1 | I(t-1)]

global Yt;
global choice;
warning('off');

p11 = normcdf(x(1,1));
[PF,PP]=filtrage(x);  % PF = Pr[St=1|It], PP = Pr(St=1|It-1)

T = size(Yt,1);
PL = zeros(T,1);
PL(T,1)= PF(T,1);% la derniere proba lissee=la derniere proba filtree

%Algorithme de Lissage
for t=T-1:-1:1
    PLi1 = (p11*PF(t,1))/(PP(t+1,1)) * PL(t+1,1);
    PLi2 = ((1-p11)*PF(t,1))/(1-PP(t+1,1))*(1-PL(t+1,1));
    PL(t,1)=PLi1+PLi2;
end

function [PR_STT1,PR_STL1]=filtrage(x)
% C'est la même fonction que vraissemblance.m, avec des outputs différents
%choice = 1 --> Y(t) = C(St) + phi1(St) + sigma_res(St)
%choice = 2 --> Y(t) = C(St) + phi1 + phi2 + sigma_res(St)
%choice = 3 --> Y(t) = C(St) + phi1 + phi2 + phi3 + phi4 + sigma_res

global Yt;

global choice;

% definition des parametres 
p11= normcdf(x(1,1));
p22= normcdf(x(2,1));
c= [x(3,1) x(4,1)]';
phi1= [x(5,1) x(6,1)]';
VAR_L= [x(7,1) x(8,1)]';

if choice == 1
    c= [x(3,1) x(4,1)]';
    phi1= [x(5,1) x(5,1)]';
    var_res = [x(6,1) x(7,1)]';
elseif choice == 2
    c= [x(3,1) x(4,1)]';
    phi1 = x(5,1);
    phi2 = x(6,1);
    var_res = [x(7,1) x(8,1)]';
elseif choice == 3
    c= [x(3,1) x(4,1)]';
    phi1 = x(5,1);
    phi2 = x(6,1);
    phi3 = x(7,1);
    phi4 = x(8,1);
    var_res = [x(9,1) x(9,1)]';
end

% une matrice de probabilites de transition initiale 
PR_TR= [p11 1-p22; 1-p11 p22];
PR_TRF= [p11; 1-p11; 1-p22; p22];

% calcul des probabilites ergodiques pour initialiser le filtre
pi1 = (1-p22)/(2-p11-p22);   % probabilité inconditionnelle P(St=1)
pi2 = 1-pi1;   % probabilité inconditionnelle P(St=2)
PROB__T = [pi1; pi2];     % vecteur des proba inconditionnelles (2x1)
PROB__= [pi1; pi1; pi2; pi2]; 

T = size(Yt,1);  % Taille de l'échantillon
LIKV=0.0; 
PR_STL1 = zeros(T);  % vecteur de taille T qui contiendra P(St|It-1)
PR_STT1 = zeros(T);  % vecteur de taille T qui contiendra P(St|It)

if choice == 1
    i = 2;
elseif choice == 2
    i = 3;
elseif choice == 3
    i = 5;
end

for J_Iter=i:T  

    PROB_DD = PR_TRF .* PROB__; % Pr[St=i,St-1=j|It-1] 4x1  (étape 1 du filtre)

    PROB_DD = PROB_DD(1:2,1)+ PROB_DD(3:4,1);   % Pr[St=i|It-1] 2x1  (étape 2 du filtre)

    PR_STL1(J_Iter,1)=PROB_DD(1,1);  % Pr(St=1|It-1) 

    if choice == 1
        F_CAST = [Yt(J_Iter,1) - c(1,1) - phi1(1,1)*Yt(J_Iter -1,1) ; Yt(J_Iter,1) - c(2,1) - phi1(2,1)*Yt(J_Iter -1,1)]; % calcul des residus et|St (2x1)
    elseif choice == 2
        F_CAST = [Yt(J_Iter,1) - c(1,1) - phi1*Yt(J_Iter -1,1) - phi2*Yt(J_Iter -2,1) ; Yt(J_Iter,1) - c(2,1) - phi1*Yt(J_Iter -1,1)]- phi2*Yt(J_Iter -2,1); % calcul des residus et|St (2x1)
    elseif choice == 3
        F_CAST = [Yt(J_Iter,1) - c(1,1) - phi1*Yt(J_Iter -1,1) - phi2*Yt(J_Iter -2,1) - phi3*Yt(J_Iter -3,1) - phi4*Yt(J_Iter -4,1) ; Yt(J_Iter,1) - c(2,1) - phi1*Yt(J_Iter -1,1)]- phi2*Yt(J_Iter -2,1) - phi3*Yt(J_Iter -3,1) - phi4*Yt(J_Iter -4,1); % calcul des residus et|St (2x1)
    end

    PR_VL= (1./sqrt(2*pi*VAR_L)) .* exp(-0.5*F_CAST.^2./VAR_L) .* PROB_DD;  % f(yt,St|It-1)  2x1 (étape 3 du filtre)

    PR_VAL= sum(PR_VL);    % f(yt|It-1), densite de yt sachant l'info passee (étape 4 du filtre)

    LIK=-1*log(PR_VAL); 
    % calcul de l'oppose de la log-vrais car fminunc mini les fonctions

    PROB__T= PR_VL ./ PR_VAL;   % Pr[St|It] 2x1 (étape 5 du filtre) 

    PR_STT1(J_Iter,1)=PROB__T(1,1);   % Pr(St=1|It) 

    PROB__ = [PROB__T(1,1);PROB__T(1,1);PROB__T(2,1);PROB__T(2,1)];  
    % nv PROB__ utilise ds l'iteration suivante 4x1

    LIKV = LIKV + LIK;

end


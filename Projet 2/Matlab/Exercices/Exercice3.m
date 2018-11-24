%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXERCICE 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
close all;
clear all;
warning('off');

SetPathExercices();

global choice;
global Yt;

global endog;
global exog_trans;

%% Question 1 : Générer une variable Z(t) de 500 observations suivant un AR(1)

%taille désiré de l'échantillon de la variable simulé
T=500;

%initialisation d'un vecteur contenant la variable simulé, avec 100
%observations de plus que prévues pour limiter l'effet des condtions
%intiales
Z= ones(T+100,1); 

%générataion du vecteur des résidus normalement et identiquement distribués
res=randn(T+100,1);

%simulation d'un processus AR(1) de T + 100 observations avec phi=0,5
for t=2:T+100
    Z(t)=0.5*Z(t-1)+res(t); 
end 

%la variable simulée Z est la variable dont vont dépendre les probabilités
%de transition variables
exog_trans = Z; 

%% Question 2 : Générer un processus MS-AR(1,1) de 500 observations

%valeurs des paramètres pour les probabilités de transition variables
alpha0 = 1.5;
alpha1 = 0.5; 
beta0 = 1;
beta1 = -0.5;

%2 vecteurs de probabilités de transation variables (p11 et p22), indéxé par le temps, fonction de
%Z(t)
p11 = normcdf(alpha0 + alpha1*Z,0,1);
p12 = 1 - p11;
p22 = normcdf(beta0 + beta1*Z,0,1);
p21 = 1 - p22;

%coefficients de départ du modèle 
c_r1 = 0.5;
c_r2 = -0.5;
phi = 0.5;
sig_r1 = 0.3;
sig_r2 = 0.2;

%Simulation d'un processus MS-AR  de 500 observations, Y, en fonction de probabilités de transition variables
%et indéxeés par le temps p11 et p22, et le vecteur d'états simulés associé, s
[Y, s] = Simul_MSAR([c_r1 c_r2 phi sig_r1 sig_r2], [p11 p22] ,500);

%La variable simulée Y est la variable à expliquer de notre modèle à
%probabilités de transition variables
endog = Y;
Yt = endog;

%% Question 3

%définition des probabilités de transition de départ pour l'estimation du
%modèle MS-AR
P = [0.9 0.2;0.1 0.8];
p11 = P(1,1);
p22 = P(2,2);

%définition de vecteur des paramètres initiaux pour l'optimisation
%numérique des coefficients du modèle MS-AR
%x0 --> paramètres initiaux du modèle à probas de transition fixes
%x0_v = paramètres initiaux du modèle à probas de transition variables

%x0 = [p11 p22 c1 c2 phi sig1 sig2]
x0 = [norminv(p11) norminv(p22) c_r1 c_r2 phi sig_r1 sig_r2]';
x0 = x0 + randn(7,1)./100;

%x0_v = [alpha0 beta0 c1 c2 phi sig1 sig2 alpha1 beta1]
x0_v = [x0;0;0];


%% Estimation par maximum de vraisemblance des paramètres du modèle Pfixe

% Options de l'Optimisation
options = optimset('TolX',10^(-10),'TolFun',10^(-10),'MaxIter',1000,'Display','final');

%Estimation du modèle MS-AR à probabilité de transitions fixes

[x,fval,code,info,g,H] = fminunc('vraissemblance_PFixe', x0, options);

%matrice de variance covariance du modèle et ecart type des coefficients
var_cov=inv(H);
sd = sqrt(diag(var_cov));

%t-stats des coefficients optimaux du modèle 
t_stat = x./sd;

%% Estimation par maximum de vraisemblance des paramètres du modèle Pvariable

% Options de l'Optimisation
options_v = optimset('TolX',10^(-10),'TolFun',10^(-10),'MaxIter',1000,'Display','final');

%Estimation du modèle MS-AR à probabilité de transitions variables

[x_v,fval_v,code_v,info_v,g_v,H_v] = fminunc('vraissemblance_Pvar', x0_v, options_v);

%matrice de variance covariance du modèle et ecart type des coefficients
var_cov_v=inv(H_v); % inverse du Hessien 
sd_v = sqrt(diag(var_cov_v)); % ecart-types des parametres estimes

%t-stats des coefficients optimaux du modèle 
t_stat_v = x_v./sd_v;

%test de significativité de student sur les coefficient du modèle
[student_PV, student_pvalue_PV, student_critikVal_PV] = MyStudentTest(t_stat_v, 0.05);

%% Affichage des résultats d'estimation

clc ;
format compact ;
disp('RESULTATS D''ESTIMATION DU MODELE');
disp('la valeur de la log-vraisemblance maximisee est');
disp([-fval -fval_v]);
disp('Test FTP versus TVTP: stat LR (p_value)');

LRobs= 2*(fval_v - fval);
pval = (1- cdf('chi2', LRobs, 2));

disp(strcat(num2str(LRobs),' (',num2str(pval),')'));
disp('----------------------------------------------');
disp('---    Proba Fixes       Proba Variables   ---');
disp('----------------------------------------------');
disp('---  estim      t-stat  estim      t-stat  ---');
disp('----------------------------------------------');
disp([[normcdf(x(1:2));x(3:7);0;0] [x./sd;0;0] x_v x_v./sd_v]);
disp('----------------------------------------------');

%% Représentation des probabilités lissées

[proba_liss,proba_filt] = lissage_PFixe(x);
[proba_liss_v,proba_filt_v] = lissage_Pvar(x_v);

%s vaut 1 dans l'état 1 et 2 dans l'état 2 
%s-1 vaut 0 dans l'état 1 et 1 dans l'état 2
%proba_**** est très proche de 0 dans l'état 2 (c'est la proba St = 1
%sachant qu'on est dans l'état 2).
%proba_**** est très proche de 1 dans l'état 1 (c'est la proba St = 1
%sachant qu'on est dans l'état 1).
%1 - proba_**** est très proche de 0 dans l'état 1 et très proche de 1 dans
%l'état 2.
%--> On compare donc s-1 à 1 - proba_**** pour comparer les états estimés
%aux états simulés

figure;
subplot(2,1,1)
hold on;
bar(s-1);
plot(1-proba_filt,'-r');
title('Probas Fixes : Etats estimés (probabilités filtrées) VS Etats simulés');
hold off;
legend('Etats simulés', 'P[St=2 | I(t)]','Pvariables');

subplot(2,1,2)
hold on;
bar(s-1);
plot(1-proba_liss,'-r');
title('Probas Fixes : Etats estimés (probabilités lissées) VS Etats simulés');
hold off;
legend('Etats simulés', '[St=2 | I(T)]','Pvariables');

figure;
subplot(2,1,1)
hold on;
bar(s-1);
plot(1-proba_filt_v,'-r');
title('Probas Variables : Etats estimés (probabilités filtrées) VS Etats simulés');
hold off;
legend('Etats simulés', '[St=2 | I(t)]','Pvariables');

subplot(2,1,2)
hold on;
bar(s-1)
plot(1-proba_liss_v,'-r');
title('Probas Variables : Etats estimés (probabilités lissées) VS Etats simulés');
hold off;
legend('Etats simulés', '[St=2 | I(T)]','Pvariables');

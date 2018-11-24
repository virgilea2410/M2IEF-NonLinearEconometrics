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

%% Question 1 : G�n�rer une variable Z(t) de 500 observations suivant un AR(1)

%taille d�sir� de l'�chantillon de la variable simul�
T=500;

%initialisation d'un vecteur contenant la variable simul�, avec 100
%observations de plus que pr�vues pour limiter l'effet des condtions
%intiales
Z= ones(T+100,1); 

%g�n�rataion du vecteur des r�sidus normalement et identiquement distribu�s
res=randn(T+100,1);

%simulation d'un processus AR(1) de T + 100 observations avec phi=0,5
for t=2:T+100
    Z(t)=0.5*Z(t-1)+res(t); 
end 

%la variable simul�e Z est la variable dont vont d�pendre les probabilit�s
%de transition variables
exog_trans = Z; 

%% Question 2 : G�n�rer un processus MS-AR(1,1) de 500 observations

%valeurs des param�tres pour les probabilit�s de transition variables
alpha0 = 1.5;
alpha1 = 0.5; 
beta0 = 1;
beta1 = -0.5;

%2 vecteurs de probabilit�s de transation variables (p11 et p22), ind�x� par le temps, fonction de
%Z(t)
p11 = normcdf(alpha0 + alpha1*Z,0,1);
p12 = 1 - p11;
p22 = normcdf(beta0 + beta1*Z,0,1);
p21 = 1 - p22;

%coefficients de d�part du mod�le 
c_r1 = 0.5;
c_r2 = -0.5;
phi = 0.5;
sig_r1 = 0.3;
sig_r2 = 0.2;

%Simulation d'un processus MS-AR  de 500 observations, Y, en fonction de probabilit�s de transition variables
%et ind�xe�s par le temps p11 et p22, et le vecteur d'�tats simul�s associ�, s
[Y, s] = Simul_MSAR([c_r1 c_r2 phi sig_r1 sig_r2], [p11 p22] ,500);

%La variable simul�e Y est la variable � expliquer de notre mod�le �
%probabilit�s de transition variables
endog = Y;
Yt = endog;

%% Question 3

%d�finition des probabilit�s de transition de d�part pour l'estimation du
%mod�le MS-AR
P = [0.9 0.2;0.1 0.8];
p11 = P(1,1);
p22 = P(2,2);

%d�finition de vecteur des param�tres initiaux pour l'optimisation
%num�rique des coefficients du mod�le MS-AR
%x0 --> param�tres initiaux du mod�le � probas de transition fixes
%x0_v = param�tres initiaux du mod�le � probas de transition variables

%x0 = [p11 p22 c1 c2 phi sig1 sig2]
x0 = [norminv(p11) norminv(p22) c_r1 c_r2 phi sig_r1 sig_r2]';
x0 = x0 + randn(7,1)./100;

%x0_v = [alpha0 beta0 c1 c2 phi sig1 sig2 alpha1 beta1]
x0_v = [x0;0;0];


%% Estimation par maximum de vraisemblance des param�tres du mod�le Pfixe

% Options de l'Optimisation
options = optimset('TolX',10^(-10),'TolFun',10^(-10),'MaxIter',1000,'Display','final');

%Estimation du mod�le MS-AR � probabilit� de transitions fixes

[x,fval,code,info,g,H] = fminunc('vraissemblance_PFixe', x0, options);

%matrice de variance covariance du mod�le et ecart type des coefficients
var_cov=inv(H);
sd = sqrt(diag(var_cov));

%t-stats des coefficients optimaux du mod�le 
t_stat = x./sd;

%% Estimation par maximum de vraisemblance des param�tres du mod�le Pvariable

% Options de l'Optimisation
options_v = optimset('TolX',10^(-10),'TolFun',10^(-10),'MaxIter',1000,'Display','final');

%Estimation du mod�le MS-AR � probabilit� de transitions variables

[x_v,fval_v,code_v,info_v,g_v,H_v] = fminunc('vraissemblance_Pvar', x0_v, options_v);

%matrice de variance covariance du mod�le et ecart type des coefficients
var_cov_v=inv(H_v); % inverse du Hessien 
sd_v = sqrt(diag(var_cov_v)); % ecart-types des parametres estimes

%t-stats des coefficients optimaux du mod�le 
t_stat_v = x_v./sd_v;

%test de significativit� de student sur les coefficient du mod�le
[student_PV, student_pvalue_PV, student_critikVal_PV] = MyStudentTest(t_stat_v, 0.05);

%% Affichage des r�sultats d'estimation

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

%% Repr�sentation des probabilit�s liss�es

[proba_liss,proba_filt] = lissage_PFixe(x);
[proba_liss_v,proba_filt_v] = lissage_Pvar(x_v);

%s vaut 1 dans l'�tat 1 et 2 dans l'�tat 2 
%s-1 vaut 0 dans l'�tat 1 et 1 dans l'�tat 2
%proba_**** est tr�s proche de 0 dans l'�tat 2 (c'est la proba St = 1
%sachant qu'on est dans l'�tat 2).
%proba_**** est tr�s proche de 1 dans l'�tat 1 (c'est la proba St = 1
%sachant qu'on est dans l'�tat 1).
%1 - proba_**** est tr�s proche de 0 dans l'�tat 1 et tr�s proche de 1 dans
%l'�tat 2.
%--> On compare donc s-1 � 1 - proba_**** pour comparer les �tats estim�s
%aux �tats simul�s

figure;
subplot(2,1,1)
hold on;
bar(s-1);
plot(1-proba_filt,'-r');
title('Probas Fixes : Etats estim�s (probabilit�s filtr�es) VS Etats simul�s');
hold off;
legend('Etats simul�s', 'P[St=2 | I(t)]','Pvariables');

subplot(2,1,2)
hold on;
bar(s-1);
plot(1-proba_liss,'-r');
title('Probas Fixes : Etats estim�s (probabilit�s liss�es) VS Etats simul�s');
hold off;
legend('Etats simul�s', '[St=2 | I(T)]','Pvariables');

figure;
subplot(2,1,1)
hold on;
bar(s-1);
plot(1-proba_filt_v,'-r');
title('Probas Variables : Etats estim�s (probabilit�s filtr�es) VS Etats simul�s');
hold off;
legend('Etats simul�s', '[St=2 | I(t)]','Pvariables');

subplot(2,1,2)
hold on;
bar(s-1)
plot(1-proba_liss_v,'-r');
title('Probas Variables : Etats estim�s (probabilit�s liss�es) VS Etats simul�s');
hold off;
legend('Etats simul�s', '[St=2 | I(T)]','Pvariables');

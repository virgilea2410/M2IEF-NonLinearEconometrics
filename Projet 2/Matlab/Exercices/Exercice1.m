%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXERCICE 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
close all;
clear all;
warning('off');

SetPathExercices();

global choice;
global Yt;

%importation et traitement des données 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = readtable('VIX.xlsx');
dates = datetime(data.Security(7:end));
dates = datenum(dates);
vix = str2double(data.VIXIndex(7:end));

%définition de la variable à expliquer du modèle Y
Yt = vix;

%taille de la série 
n = size(Yt,1);

%% Question 1 : Estimation processus AR(P) linéaire

%optimisation du retard p des variables explicatives du modèle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_optim = OptimLagByAIC();
p_optim_bis = OptimLagByAutoc();

%nombre de paramètres estimés du modèle, nombres d'observations du modèle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k_AR = p_optim + 1;
nobs = n - p_optim;

%estimation AR par max de vraisemblance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%méthode 1 
[AR_phi, AR_t_stat, AR_var_res, AR_resids, AR_aic, AR_logvraiss, AR_std_coeff] = AR_LogV(p_optim);

%méthode 2
choice = p_optim;
[AR_phi_2, AR_t_stat_2, AR_var_res_2, AR_resids_2, AR_aic_2, AR_logvraiss_2, AR_std_coeff_2] = AR_LogV_2(p_optim);

%test de diagnostics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%test de significativité des coefficients

[test_student_AR, student_AR_Pval, student_AR_critikVal] = MyStudentTest(AR_t_stat, 0.05);

%test d'autocorrélation des résidus 

%méthode 1
[AR_test_lb, AR_pVal_lb, AR_stat_lb, AR_critikVal_lb] = lbqtest(AR_resids, 'lags', 20, 'alpha', 0.05);

%méthode 2
[AR_test_lb_2, AR_pVal_lb_2, AR_stat_lb_2, AR_critikVal_lb_2] = MyLbTest(AR_resids, k_AR);

%test d'effet ARCH sur les résidus 

%méthode 1 
[AR_test_ARCH, AR_pVal_ARCH, AR_stat_ARCH, AR_critikVal_ARCH] = lbqtest(AR_resids.^2, 'lags', 20, 'alpha', 0.05);

%méthode 2 
[AR_test_ARCH_2, AR_pVal_ARCH_2, AR_stat_ARCH_2, AR_critikVal_ARCH_2] = MyLbTest(AR_resids.^2, k_AR);

%test de normalité sur les résidus 

[AR_test_jb, AR_pVal_jb, AR_stat_jb, AR_critikVal_jb] = jbtest(AR_resids);

%Affichage des résultats de l'estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

names = [string('C');string('Phi 1'); string('Phi 2'); string('Etype Res')];
z = zeros(2+p_optim,1);
AR_theta = [AR_phi_2;AR_var_res_2];

format compact ;
disp('--------------------------------------------------------------');
disp('RESULTATS D''ESTIMATION DU MODELE');
disp('--------------------------------------------------------------');
disp('--- Valeurs     Depart      estim         std      t-stat ----');
disp('--------------------------------------------------------------');
disp([names z AR_theta AR_std_coeff_2 AR_t_stat_2]);
disp('--------------------------------------------------------------');
disp('la valeur de la log-vraisemblance maximisee est');
disp(-AR_logvraiss);
disp('Les résidus sont non auto-corrélés :');
if AR_test_lb == 0
    disp('Vrai !');
else 
    disp('Faux !');
end
disp ('Les coefficients sont significatifs : ');
if test_student_AR == 1
    disp('Vrai !');
else 
    disp('Faux !');
end
disp('--------------------------------------------------------------');

%% Question 2 : Estimer le modèle MS-AR associé, avec St la variable de changement de régime suivant une chaine de Markov d'ordre 1

%initialisation des paramètres de l'optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%probabilités de transition & coefficients de départ du modèle
P = [0.9 0.2;0.1 0.8];
p11 = P(1,1);
p22 = P(2,2);
c_r1 = AR_phi(1,1);
c_r2 = AR_phi(1,1);
phi1 = AR_phi(2,1);
phi2 = AR_phi(3,1);
sig_r1 = sqrt(AR_var_res);
sig_r2 = sqrt(AR_var_res);

k_MSAR = 8;

%vecteur de paramètres intitiaux pour l'optimisation
x0 = [norminv(p11) norminv(p22) c_r1 c_r2 phi1 phi2 sig_r1 sig_r2]';
x0 = x0 + randn(8,1)./100 *2;

%Estimation du modèle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
choice = 2;

%Méthode 1
[MSAR_theta, MSAR_std_coeff, MSAR_t_stat, MSAR_logvraiss] = maxVraissemblance(x0);

%Méthode 2
[MSAR_theta_2, MSAR_std_coeff_2, MSAR_t_stat_2, MSAR_logvraiss_2] = maxVraissemblance_2();

%définition des vecteurs de proba optimales et des coefficients optimaux
MSAR_prob_optim = normcdf(MSAR_theta(1:2,1));
MSAR_phi_optim = MSAR_theta(3:end,1);

%cacul des probabilités lissées
[PL, PF, PP] = lissage(MSAR_theta);

%caractéristiques des 2 régimes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DureMoy1 = 1/(1-MSAR_prob_optim(1,1));
DureMoy2 = 1/(1-MSAR_prob_optim(2,1));
[~,info_Y1, info_Y2] = Res_Gen_MSAR(MSAR_theta,PF);
    

%Etats estimés 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
subplot(2,1,1);
plot(dates, PF);
title('Probabilités Filtrées');
datetick('x', 'mm-yy', 'keepticks');
legend('P[St=1 | I(t)]');
subplot(2,1,2);
plot(dates, PL);
title('Probabilités Lissées');
datetick('x', 'mm-yyyy', 'keepticks');
legend('P[St=1 | I(T)]');

figure;
bar(dates, PL);
title('Probabilités Lissées');
datetick('x', 'mm-yyyy', 'keepticks');
legend('P[St=1 | I(T)]');


%Affichage des résultats d'estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

names = [string('C (R1)'); string('C (R2)'); string('Phi 1'); string('Phi 2'); string('Etype Res (R1)'); string('Etype Res (R2)')];

format compact ;
disp('----------------------------------------------------------------------');
disp('RESULTATS D''ESTIMATION DU MODELE');
disp('----------------------------------------------------------------------');
disp('--- Valeurs             Depart        estim          std      t-stat ----');
disp('----------------------------------------------------------------------');
disp([names x0(3:end,1) MSAR_theta(3:end,1) MSAR_std_coeff(3:end,1) MSAR_t_stat(3:end,1)]);
disp('----------------------------------------------------------------------');
disp('la valeur de la log-vraisemblance maximisee est');
disp(MSAR_logvraiss);
disp('----------------------------------------------------------------------');
disp('----------------------------------------------------------------------');

%% Question 3 : Vérifier à l'aide de tests l'absence d'autocorrélation des erreurs

%Calcul des résidus généralisés du modèle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%résidus généralisés
resids_MSAR = Res_Gen_MSAR(MSAR_theta, PF);

%SCR généralisé
SCR_MSAR = resids_MSAR' * resids_MSAR;

%variance résiduelle généralisé
var_res_gen_MSAR = SCR_MSAR/(nobs-6);

%Test de diagnostics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%test d'autocorrélation

%méthode 1 
[lb_MSAR, pVal_lb_MSAR, stat_lb_MSAR, critikVal_lb_MSAR] = lbqtest(resids_MSAR, 'lags', 20, 'alpha', 0.05);

%méthode 2
[lb_MSAR_2, pVal_lb_MSAR_2, stat_lb_MSAR_2, critikVal_lb_MSAR_2] = MyLbTest(resids_MSAR, k_MSAR);
%--> Ne donne pas le même résulat car la fonction lbqtest considère un zone
%de rejet pour la loi du chi2 dont le degré de libérté est de 20, tandis
%que dans MyLbTest, le degré de libérté est de (20-k) = 12

%test d'absence d'effets ARCH

%méthode 1 
[ARCH_MSAR, pVal_ARCH_MSAR,  stat_ARCH_MSAR, critikVal_ARCH_MSAR] = lbqtest(resids_MSAR.^2, 'lags', 20, 'alpha', 0.05);

%méthode 2
[ARCH_MSAR_2, pVal_ARCH_MSAR_2,  stat_ARCH_MSAR_2, critikVal_ARCH_MSAR_2] = MyLbTest(resids_MSAR.^2, k_MSAR);

%test de normalité

[jb_MSAR, pVal_jb_MSAR] = jbtest(resids_MSAR);

%test de significativité des coefficients 

[test_student_MSAR, student_MSAR_Pval, student_MSAR_critikVal] = MyStudentTest(MSAR_t_stat, 0.05);

%Affichage des résultats des tests de diagnostics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--------------------------------------------------------------');
disp('TESTS DE DIAGNOSTIC DU MODELE');
disp('Les résidus sont non auto-corrélés :');
if lb_MSAR == 0
    disp('Vrai !');
else 
    disp('Faux !');
end
disp('Les résidus ne présentent pas d effets ARCH :');
if ARCH_MSAR == 0
    disp('Vrai !');
else 
    disp('Faux !');
end
disp('Les résidus peuvent être considérés comme normaux :');
if jb_MSAR == 0
    disp('Vrai !');
else 
    disp('Faux !');
end
disp ('Les coefficients sont significatifs : ');
if test_student_MSAR == 1
    disp('Vrai !');
else 
    disp('Faux !');
end
disp('--------------------------------------------------------------');
disp('--------------------------------------------------------------');

%% Question 4 : Comparer les critères AIC et SIC des modèles AR et MSAR, et commentez

%Définition du nombre de paramètres estimés de chaque modèle 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AR_k = p_optim + 1;
MSAR_k = p_optim + 4;

%calcul des critères d'Akaike (AIC) et de Scwhartz (SIC) des modèle AR et MS-AR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AR_AIC = -2 * (-AR_logvraiss/nobs) + (2*(AR_k))/nobs;
MSAR_AIC = -2 * (-MSAR_logvraiss/nobs) + (2*(MSAR_k))/nobs;

AR_SIC = -2 * (-AR_logvraiss/nobs) + (AR_k*log(nobs))/nobs;
MSAR_SIC = -2 * (-MSAR_logvraiss/nobs) + (MSAR_k*log(nobs))/nobs;

%vérification des résultats à l'aide d'une fonction Matlab
verif_ar_aic = aicbic(AR_logvraiss,AR_k,nobs);
verif_msar_aic = aicbic(-MSAR_logvraiss,MSAR_k,nobs);

%Afiichage de la comparaisons des critères d'informations AIC et SIC des 2
%modèles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format compact ;
disp('----------------------------------------------------------------------');
disp('COMPARAISON DES CRITERES D AKAIKE DES MODELES MSAR ET AR');
disp('----------------------------------------------------------------------');
disp('--- Modèle MSAR : ----');
disp('--- AIC         SIC : ----');
disp('----------------------------------------------------------------------');
disp([MSAR_AIC MSAR_SIC]);
disp('----------------------------------------------------------------------');
disp('--- Modèle AR : ----');
disp('--- AIC         SIC : ----');
disp('----------------------------------------------------------------------');
disp([AR_AIC AR_SIC]);
disp('----------------------------------------------------------------------');
disp('----------------------------------------------------------------------');

%% Question 5 : Effectuez le test d'un modèle AR contre un modèle MS-AR 

[MSAR_test_linearite, MSAR_StatDuTest_linearite, MSAR_pcritique_linearite, MSAR_critikVal_linearite] = myMSARTest(P);
Yt = vix;
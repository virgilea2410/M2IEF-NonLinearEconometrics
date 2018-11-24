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

%importation et traitement des donn�es 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = readtable('VIX.xlsx');
dates = datetime(data.Security(7:end));
dates = datenum(dates);
vix = str2double(data.VIXIndex(7:end));

%d�finition de la variable � expliquer du mod�le Y
Yt = vix;

%taille de la s�rie 
n = size(Yt,1);

%% Question 1 : Estimation processus AR(P) lin�aire

%optimisation du retard p des variables explicatives du mod�le
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_optim = OptimLagByAIC();
p_optim_bis = OptimLagByAutoc();

%nombre de param�tres estim�s du mod�le, nombres d'observations du mod�le
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k_AR = p_optim + 1;
nobs = n - p_optim;

%estimation AR par max de vraisemblance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%m�thode 1 
[AR_phi, AR_t_stat, AR_var_res, AR_resids, AR_aic, AR_logvraiss, AR_std_coeff] = AR_LogV(p_optim);

%m�thode 2
choice = p_optim;
[AR_phi_2, AR_t_stat_2, AR_var_res_2, AR_resids_2, AR_aic_2, AR_logvraiss_2, AR_std_coeff_2] = AR_LogV_2(p_optim);

%test de diagnostics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%test de significativit� des coefficients

[test_student_AR, student_AR_Pval, student_AR_critikVal] = MyStudentTest(AR_t_stat, 0.05);

%test d'autocorr�lation des r�sidus 

%m�thode 1
[AR_test_lb, AR_pVal_lb, AR_stat_lb, AR_critikVal_lb] = lbqtest(AR_resids, 'lags', 20, 'alpha', 0.05);

%m�thode 2
[AR_test_lb_2, AR_pVal_lb_2, AR_stat_lb_2, AR_critikVal_lb_2] = MyLbTest(AR_resids, k_AR);

%test d'effet ARCH sur les r�sidus 

%m�thode 1 
[AR_test_ARCH, AR_pVal_ARCH, AR_stat_ARCH, AR_critikVal_ARCH] = lbqtest(AR_resids.^2, 'lags', 20, 'alpha', 0.05);

%m�thode 2 
[AR_test_ARCH_2, AR_pVal_ARCH_2, AR_stat_ARCH_2, AR_critikVal_ARCH_2] = MyLbTest(AR_resids.^2, k_AR);

%test de normalit� sur les r�sidus 

[AR_test_jb, AR_pVal_jb, AR_stat_jb, AR_critikVal_jb] = jbtest(AR_resids);

%Affichage des r�sultats de l'estimation
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
disp('Les r�sidus sont non auto-corr�l�s :');
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

%% Question 2 : Estimer le mod�le MS-AR associ�, avec St la variable de changement de r�gime suivant une chaine de Markov d'ordre 1

%initialisation des param�tres de l'optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%probabilit�s de transition & coefficients de d�part du mod�le
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

%vecteur de param�tres intitiaux pour l'optimisation
x0 = [norminv(p11) norminv(p22) c_r1 c_r2 phi1 phi2 sig_r1 sig_r2]';
x0 = x0 + randn(8,1)./100 *2;

%Estimation du mod�le
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
choice = 2;

%M�thode 1
[MSAR_theta, MSAR_std_coeff, MSAR_t_stat, MSAR_logvraiss] = maxVraissemblance(x0);

%M�thode 2
[MSAR_theta_2, MSAR_std_coeff_2, MSAR_t_stat_2, MSAR_logvraiss_2] = maxVraissemblance_2();

%d�finition des vecteurs de proba optimales et des coefficients optimaux
MSAR_prob_optim = normcdf(MSAR_theta(1:2,1));
MSAR_phi_optim = MSAR_theta(3:end,1);

%cacul des probabilit�s liss�es
[PL, PF, PP] = lissage(MSAR_theta);

%caract�ristiques des 2 r�gimes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DureMoy1 = 1/(1-MSAR_prob_optim(1,1));
DureMoy2 = 1/(1-MSAR_prob_optim(2,1));
[~,info_Y1, info_Y2] = Res_Gen_MSAR(MSAR_theta,PF);
    

%Etats estim�s 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
subplot(2,1,1);
plot(dates, PF);
title('Probabilit�s Filtr�es');
datetick('x', 'mm-yy', 'keepticks');
legend('P[St=1 | I(t)]');
subplot(2,1,2);
plot(dates, PL);
title('Probabilit�s Liss�es');
datetick('x', 'mm-yyyy', 'keepticks');
legend('P[St=1 | I(T)]');

figure;
bar(dates, PL);
title('Probabilit�s Liss�es');
datetick('x', 'mm-yyyy', 'keepticks');
legend('P[St=1 | I(T)]');


%Affichage des r�sultats d'estimation
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

%% Question 3 : V�rifier � l'aide de tests l'absence d'autocorr�lation des erreurs

%Calcul des r�sidus g�n�ralis�s du mod�le
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%r�sidus g�n�ralis�s
resids_MSAR = Res_Gen_MSAR(MSAR_theta, PF);

%SCR g�n�ralis�
SCR_MSAR = resids_MSAR' * resids_MSAR;

%variance r�siduelle g�n�ralis�
var_res_gen_MSAR = SCR_MSAR/(nobs-6);

%Test de diagnostics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%test d'autocorr�lation

%m�thode 1 
[lb_MSAR, pVal_lb_MSAR, stat_lb_MSAR, critikVal_lb_MSAR] = lbqtest(resids_MSAR, 'lags', 20, 'alpha', 0.05);

%m�thode 2
[lb_MSAR_2, pVal_lb_MSAR_2, stat_lb_MSAR_2, critikVal_lb_MSAR_2] = MyLbTest(resids_MSAR, k_MSAR);
%--> Ne donne pas le m�me r�sulat car la fonction lbqtest consid�re un zone
%de rejet pour la loi du chi2 dont le degr� de lib�rt� est de 20, tandis
%que dans MyLbTest, le degr� de lib�rt� est de (20-k) = 12

%test d'absence d'effets ARCH

%m�thode 1 
[ARCH_MSAR, pVal_ARCH_MSAR,  stat_ARCH_MSAR, critikVal_ARCH_MSAR] = lbqtest(resids_MSAR.^2, 'lags', 20, 'alpha', 0.05);

%m�thode 2
[ARCH_MSAR_2, pVal_ARCH_MSAR_2,  stat_ARCH_MSAR_2, critikVal_ARCH_MSAR_2] = MyLbTest(resids_MSAR.^2, k_MSAR);

%test de normalit�

[jb_MSAR, pVal_jb_MSAR] = jbtest(resids_MSAR);

%test de significativit� des coefficients 

[test_student_MSAR, student_MSAR_Pval, student_MSAR_critikVal] = MyStudentTest(MSAR_t_stat, 0.05);

%Affichage des r�sultats des tests de diagnostics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--------------------------------------------------------------');
disp('TESTS DE DIAGNOSTIC DU MODELE');
disp('Les r�sidus sont non auto-corr�l�s :');
if lb_MSAR == 0
    disp('Vrai !');
else 
    disp('Faux !');
end
disp('Les r�sidus ne pr�sentent pas d effets ARCH :');
if ARCH_MSAR == 0
    disp('Vrai !');
else 
    disp('Faux !');
end
disp('Les r�sidus peuvent �tre consid�r�s comme normaux :');
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

%% Question 4 : Comparer les crit�res AIC et SIC des mod�les AR et MSAR, et commentez

%D�finition du nombre de param�tres estim�s de chaque mod�le 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AR_k = p_optim + 1;
MSAR_k = p_optim + 4;

%calcul des crit�res d'Akaike (AIC) et de Scwhartz (SIC) des mod�le AR et MS-AR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AR_AIC = -2 * (-AR_logvraiss/nobs) + (2*(AR_k))/nobs;
MSAR_AIC = -2 * (-MSAR_logvraiss/nobs) + (2*(MSAR_k))/nobs;

AR_SIC = -2 * (-AR_logvraiss/nobs) + (AR_k*log(nobs))/nobs;
MSAR_SIC = -2 * (-MSAR_logvraiss/nobs) + (MSAR_k*log(nobs))/nobs;

%v�rification des r�sultats � l'aide d'une fonction Matlab
verif_ar_aic = aicbic(AR_logvraiss,AR_k,nobs);
verif_msar_aic = aicbic(-MSAR_logvraiss,MSAR_k,nobs);

%Afiichage de la comparaisons des crit�res d'informations AIC et SIC des 2
%mod�les
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format compact ;
disp('----------------------------------------------------------------------');
disp('COMPARAISON DES CRITERES D AKAIKE DES MODELES MSAR ET AR');
disp('----------------------------------------------------------------------');
disp('--- Mod�le MSAR : ----');
disp('--- AIC         SIC : ----');
disp('----------------------------------------------------------------------');
disp([MSAR_AIC MSAR_SIC]);
disp('----------------------------------------------------------------------');
disp('--- Mod�le AR : ----');
disp('--- AIC         SIC : ----');
disp('----------------------------------------------------------------------');
disp([AR_AIC AR_SIC]);
disp('----------------------------------------------------------------------');
disp('----------------------------------------------------------------------');

%% Question 5 : Effectuez le test d'un mod�le AR contre un mod�le MS-AR 

[MSAR_test_linearite, MSAR_StatDuTest_linearite, MSAR_pcritique_linearite, MSAR_critikVal_linearite] = myMSARTest(P);
Yt = vix;
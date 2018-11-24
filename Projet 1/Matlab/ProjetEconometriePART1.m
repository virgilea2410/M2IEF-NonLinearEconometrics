%%                           Mise en place

clear all;
close all;
clc;

%Importation des données
data = readtable('VIX.xlsx');

%Suppression des lignes inutiles et définitions des tableaux dates et VIXclose
dates = datetime(data.Security(7:end));

VIXclose = str2double(data.VIXIndex(7:end));
data = VIXclose;

disp('--------------------- MISE EN PLACE ---------------------');
disp('Voici les données du VIX que nous allons étudier : ');
plot(datenum(dates), VIXclose);
title('Evolutions du VIX entre 1997 et 2017');
legend('VIX');
xlim([datenum('24/10/1997', 'dd/mm/yyyy') datenum('06/10/2017','dd/mm/yyyy')]);
xlabel('temps');
ylabel('USD');
datetick('x', 'mm/yy', 'keepticks');

%%                          Question 1

%On choisit entre les 3 modèles du test de Dickey Fuller Augmenté
[ModeleNo, statsADF, statsPP] = TestStationnarite(data);

ADFPval3 = statsADF(1,1);
ADFPval2 = statsADF(1,2);
ADFPval1 = statsADF(1,3);

PPPval3 = statsPP(1,1);
PPPval2 = statsPP(1,2);
PPPval1 = statsPP(1,3);

tstat3 = statsADF(2,1);
tstat2 = statsADF(2,2);
tstat1 = statsADF(2,3);

ppstat3 = statsPP(2,1);
ppstat2 = statsPP(2,2);
ppstat1 = statsPP(2,3);

ADFstat3 = statsADF(3,1);
ADFstat2 = statsADF(3,2);
ADFstat1 = statsADF(3,3);

disp('--------- Question Numéro 1 : Test de stationnarité --------');
disp('------------------------------------------------------------');
disp('------------------------------------------------------------');
disp('------------------------ Test ADF --------------------------');
disp('------ N° Modele ------- P Critique -------- Racine Unitaire -- T-stat --------')
disp('----------- (' + string(3) + ') --------------- (' + string(ADFPval3) + ') ------------- (' + string(tstat3) + ') ---------- ('+ string(ADFstat3) +') ------');
disp('----------- (' + string(2) + ') --------------- (' + string(ADFPval2) + ') ------------- (' + string(tstat2) + ') ---------- ('+ string(ADFstat2) +') ------');
disp('----------- (' + string(1) + ') --------------- (' + string(ADFPval1) + ') ------------- (' + string(tstat1) + ') ---------- ('+ string(ADFstat1) +') ------');
disp('------------------------------------------------------------');
disp('------------------------------------------------------------');
disp('------------------------ Test PP ---------------------------');
disp('------ N° Modele ------- P Critique -------- Racine Unitaire --------')
disp('----------- (' + string(3) + ') --------------- (' + string(PPPval3) + ') ------------- (' + string(ppstat3) + ') ---------');
disp('----------- (' + string(2) + ') --------------- (' + string(PPPval2) + ') ------------- (' + string(ppstat2) + ') ---------');
disp('----------- (' + string(1) + ') --------------- (' + string(PPPval1) + ') ------------- (' + string(ppstat1) + ') ---------');
disp('------------------------------------------------------------');
disp('------------------------------------------------------------');
disp('--------------- Choix du modèle numéro : -------------------');
disp('----------------- (' + string(ModeleNo) + ') --------------------------');


%%                          Question 2 

%On étudie l'autocorrélation générale et partielle de la variable expliquée
autoc = autocorr(data);
pautoc = parcorr(data);

%On cherche le nombre de retards p et q du modèle ARMA(p,q) qui minimise le
%BIC
[p_optim, q_optim] = ARMA_FIT(data);

%On estime notre modèle ARMA optimal
[ARMAphi, ARMAresid, lbstat, jbstat] = ARMA(data, p_optim, q_optim);

%On estime le modèle AR avec un peu plus de retards pour p
p_AR = p_optim + q_optim;

[ARphi, ARSCR, var_res_AR, ~, var_coeff_AR, t_stat_AR] = AR(data,p_AR);

%verification de nos estimations

%Modèle AR
%opt_modl = arima(p_AR, 0, 0);
%[opt_est_mdlAR, opt_est_vcovAR] = estimate(opt_modl, data, 'print', false);
%phi_est_AR = [opt_est_mdlAR.Constant opt_est_mdlAR.AR ]';

%Modèle ARMA
%opt_modlARMA = arima(p_optim, 0, q_optim);
%[opt_est_mdlARMA, opt_est_vcovARMA] = estimate(opt_modlARMA, data, 'print', false);
%phi_est_ARMA = [opt_est_mdlARMA.Constant opt_est_mdlARMA.AR opt_est_mdl.MA ]';

disp('------- Question Numéro 2 : Estimation du modèle ARMA ------');
disp('------------------------------------------------------------');
disp('On étudie l autocorrélation générale et partielle de la variable expliquée');
disp('Voir graphique "Sample Autocorrelation Function" et "Sample Partial Autocorrelation Function"');
figure;
autocorr(data);
figure;
parcorr(data);
disp('------------------------------------------------------------');
disp('--------------------- p et q optimaux ----------------------');
disp('--------- p optimal -------------------- q optimal ---------')
disp('----------- ' + string(p_optim) + ' -------------------------- ' + string(q_optim) + ' -----------------');
disp('------------------------------------------------------------');
disp('------------------------------------------------------------');
disp('----------- Coefficients du modèle ARMA(p*,q*) -------------');
disp('--------- Coefficients ----- Ecart type ------ T-stat ------');
disp('Constante --' + string(ARphi(1,1)) + '---------------' + string(sqrt(var_coeff_AR(1,1))) + '------------' + string(t_stat_AR(1,1)) + '---------');
disp('X(t-1) -----' + string(ARphi(2,1)) + '-------------' + string(sqrt(var_coeff_AR(2,2))) + '-----------' + string(t_stat_AR(2,1)) + '---------');
disp('X(t-2) -----' + string(ARphi(3,1)) + '-------------' + string(sqrt(var_coeff_AR(3,3))) + '-----------' + string(t_stat_AR(3,1)) + '---------');
disp('------------------------------------------------------------');
disp('------------------------------------------------------------');
disp('--------------- Stats du test de Ljung Box -----------------');
disp('--------------------------' + string(lbstat) + '-------------------------------');
disp('-------------- Stats du test de Jarque Bera ----------------');
disp('--------------------------' + string(jbstat) + '-------------------------------');

%%                              Question 3

%Initialisation du nombre de retards pour la variable de changement de
%régime, d
d = 1;

%définition du seuil de changement de régime optimal 
c_optim = SETAR_FIT(data,p_AR, d);

%Régression 
[SETARphi, ~, SETARvar_res, SETARvar_phi, SETARt_stat, SETARpval] = SETAR(data,p_AR, d, c_optim);

%calcul des proportions d'occurences entre nos régimes 1 et 2.
regime1 = data < c_optim;
regime2 = data > c_optim;

proportionR1 = sum(regime1) / (sum(regime1) + sum(regime2));
proportionR2 = sum(regime2) / (sum(regime1) + sum(regime2));

figure;
bar(datenum(dates), regime2);
datetick('x', 'mm/yy', 'keepticks');
title('Occurences du régime 2 entre 1997 et 2017');

disp('------- Question Numéro 3 : Estimation du modèle SETAR ------');
disp('------------------------------------------------------------');
disp('------------------------------------------------------------');
disp('------------------------ c optimal -----------------------');
disp('----------------------------' + string(c_optim) + '---------------------')
disp('------------------------------------------------------------');
disp('----------- Coefficients du modèle SETAR(p*,c*) -------------');
disp('--------- Coefficients ----- Ecart type ------ T-stat ------');
disp('Régime 1 -----------------------------------------------------');
disp('Constante --' + string(SETARphi(1,1)) + '---------------' + string(SETARvar_phi(1,1)) + '------------' + string(SETARt_stat(1,1)) + '---------');
disp('X(t-1) -----' + string(SETARphi(2,1)) + '-------------' + string(SETARvar_phi(2,2)) + '-----------' + string(SETARt_stat(2,1)) + '---------');
disp('X(t-2) -----' + string(SETARphi(3,1)) + '-------------' + string(SETARvar_phi(3,3)) + '-----------' + string(SETARt_stat(3,1)) + '---------');
disp('Régime 2 -----------------------------------------------------');
disp('Constante --' + string(SETARphi(4,1)) + '---------------' + string(SETARvar_phi(4,4)) + '------------' + string(SETARt_stat(4,1)) + '---------');
disp('X(t-1) -----' + string(SETARphi(5,1)) + '-------------' + string(SETARvar_phi(5,5)) + '-----------' + string(SETARt_stat(5,1)) + '---------');
disp('X(t-2) -----' + string(SETARphi(6,1)) + '-------------' + string(SETARvar_phi(6,6)) + '-----------' + string(SETARt_stat(6,1)) + '---------');
disp('------------------------------------------------------------');
disp('------------------------------------------------------------');
disp('--------------- Proportions d  occurences du régime 1 et 2 -----------------');
disp('---------------- Régime 1 ------------------------- Régime 2 --------------');
disp('-------------------' + string(proportionR1) + '--------------------' +  string(proportionR2) + '-----------');

%%                              Question 4

[decisionSETAR, pvalSETAR, statDuTest, Fstat] = TestSETAR(data,p_AR,d,c_optim);

disp('--------------------- Question Numéroe 4 : Test d un modèle AR contre un modèle SETAR ---------------------');
disp('--------------------- Résultat du test ---------------------');
disp('Statistique du test : ' + string(statDuTest));
disp('Distribtiion de la F-stat : ');
disp(Fstat);
disp('--------------------- Résultat du test ---------------------');
if decisionSETAR ==1
    disp('l hypothèse nulle de linéarité du modèle est rejeté, nous sommes en présence d un modèle SETAR non linéaire');
else
    disp('l hypothèse nulle de linéarité du modèle n est pas rejeté, nous sommes en présence d un modèle AR linéaire');
end


%%                              Question 5 

[d_optim, phi_optim_d, var_res_optim_d, decisionD, pvalD, c_optim_d] = SETAR_FIT_D(data,p_AR,c_optim);
[final_phi, final_resids, final_var_res, final_var_phi, final_t_stat, final_pval] = SETAR(data, p_AR, d_optim, c_optim_d);

disp('--------------------- Question Numéro 5 : esimation du décalage temporel optimal pour la variable de changement de régime ---------------------');
disp('--------------------- le décalage temporel d optimal pour la variable de changement de régime est de : ');
disp(d_optim);
disp('------------------------------------------------------------------------------------------------------------');
disp('----------- Coefficients du modèle SETAR(p*,c*,d*) -------------');
disp('--------- Coefficients ----- Ecart type ------ T-stat ------');
disp('Régime 1 -----------------------------------------------------');
disp('Constante --' + string(final_phi(1,1)) + '---------------' + string(final_var_phi(1,1)) + '------------' + string(final_t_stat(1,1)) + '---------');
disp('X(t-1) -----' + string(final_phi(2,1)) + '-------------' + string(final_var_phi(2,2)) + '-----------' + string(final_t_stat(2,1)) + '---------');
disp('X(t-2) -----' + string(final_phi(3,1)) + '-------------' + string(final_var_phi(3,3)) + '-----------' + string(final_t_stat(3,1)) + '---------');
disp('Régime 2 -----------------------------------------------------');
disp('Constante --' + string(final_phi(4,1)) + '---------------' + string(final_var_phi(4,4)) + '------------' + string(final_t_stat(4,1)) + '---------');
disp('X(t-1) -----' + string(final_phi(5,1)) + '-------------' + string(final_var_phi(5,5)) + '-----------' + string(final_t_stat(5,1)) + '---------');
disp('X(t-2) -----' + string(final_phi(6,1)) + '-------------' + string(final_var_phi(6,6)) + '-----------' + string(final_t_stat(6,1)) + '---------');
disp('------------------------------------------------------------');
disp('la p-value du test AR contre SETAR avec le d optimal est de : ');
disp(pvalD(d_optim,1));
disp('la variance résiduelle associée au d optimal est la suivante : ');
disp(var_res_optim_d);
disp('le niveau du seuil de changement de régime optimal c* associée au d optimal est le suivant : ');
disp(c_optim_d);

%%                              Question 6

[LM_stat, seuil_critique, proba_critique, F_stat, seuil_critique_F, proba_critique_F] = TestSTAR(data,p_AR,d_optim);

disp('--------------------- Question Numéro 6 : test d un modèle AR contre un modèle STAR ---------------------');
disp('la LM du test du Multiplicateur de Lagrange (AR contre STAR) est de : ');
disp(LM_stat);
disp('Le seuil critique du test du multiplicateur de Lagrange est le suivant : ');
disp(seuil_critique);
disp('La probabilité critique du test du multiplicateur de Lagrange est la suivante : ');
disp(proba_critique);
disp('la F du test de Fisher (AR contre STAR) est de : ');
disp(F_stat);
disp('Le seuil critique du test de Fisher est le suivant : ');
disp(seuil_critique_F);
disp('La probabilité critique de Fisher est la suivante : ');
disp(proba_critique_F);

if (LM_stat > seuil_critique)
    disp('On rejette l hypothèse nulle d un modèle AR pour le test du multiplicateur de lagrange');
else
    disp('On ne rejette pas l hypothèse nulle d un modèle AR pour le test du multiplicateur de lagrange');
end

if (F_stat > seuil_critique_F)
    disp('On rejette l hypothèse nulle d un modèle AR pour le test de Fisher');
else
    disp('On ne rejette pas l hypothèse nulle d un modèle AR pour le test de Fisher');
end

   
    

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXERCICE 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
close all;
clear all;
warning('off');

SetPathExercices();

global choice;
global Yt;

%Importation des données
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hamilton = readtable("Hamilton.xlsx");
GDP = hamilton.GNP;
dates = hamilton.Var1;
dates = datenum(dates);
n = size(GDP,1);

%% Question 1 : Estimation du modèle AR(4), avec saut de la constante exclusivement 

choice = 3;
Yt = GDP;

%définition de probabilités optimales + paramètres initiaux de l'énoncé du projet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p11 = 0.9;
p22 = 0.75;

c_r1 = 1.16;
c_r2 = -0.36;
phi1 = 0;
phi2 = 0;
phi3 = 0;
phi4 = 0;
var_res = 0.77;

%Estimation du modèle MS-AR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%vecteur des paramètres intiaux
x0_hamilton = [norminv(p11) norminv(p22) c_r1 c_r2 phi1 phi2 phi3 phi4 var_res]';
x0_hamilton = x0_hamilton + randn(9,1)./100;

%estimation du modèle par max de vraisemblance
[hamilton_theta, hamilton_std_coeff, hamilton_t_stat, hamilton_logvraiss] = maxVraissemblance(x0_hamilton);

%paramètres optimaux du modèle
p11_optim = normcdf(hamilton_theta(1,1));
p22_optim = normcdf(hamilton_theta(2,1));
c_optim = [hamilton_theta(3,1) hamilton_theta(4,1)]';
phi1_optim = hamilton_theta(5,1);
phi2_optim = hamilton_theta(6,1);
phi3_optim = hamilton_theta(7,1);
phi4_optim = hamilton_theta(8,1);
std_res_optim = hamilton_theta(9,1);

%nombre de paramètres estimés du modèle
k = 9;

%t-stat des coefficients optimaux estimés 
t_stat = hamilton_theta(3:end-1,1) ./ hamilton_std_coeff(3:end-1, 1);

%probabilités lissées et filtrées
[PL_hamilton, PF_hamilton, PP_hamilton] = lissage(hamilton_theta);

%Affichage des etats estimés 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,1,1);
plot(dates, PF_hamilton);
title('probabilités filtrées');
datetick('x','mm-yyyy', 'keepticks');

subplot(2,1,2);
plot(dates, PL_hamilton);
title('probabilités lissées');
datetick('x','mm-yyyy', 'keepticks');

%Affichage des résultats de l'estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

names = [string('C (R1)');string('C (R2)');string('Phi 1'); string('Phi 2');string('Phi 3');string('Phi 4'); string('Etype Res')];
z = zeros(3+4,1);

format compact ;
disp('------------------------------------------------------------------------------');
disp('RESULTATS D''ESTIMATION DU MODELE');
disp('------------------------------------------------------------------------------');
disp('--- Valeurs         Depart          estim            std         t-stat ----');
disp('------------------------------------------------------------------------------');
disp([names x0_hamilton(3:end,1) hamilton_theta(3:end,1) hamilton_std_coeff(3:end,1) hamilton_t_stat(3:end,1)]);
disp('------------------------------------------------------------------------------');
disp('la valeur de la log-vraisemblance maximisee est');
disp(hamilton_logvraiss);
disp('------------------------------------------------------------------------------');
disp('------------------------------------------------------------------------------');

%% Question 2 : Décrire les caractéristiques des deux régimes obtenus (taux de croissance moyen du PIB dans les deux états, durées moyennes des deux régimes ...)

%durées moyennes des deux régimes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DureMoy1 = 1/(1-p11_optim);
DureMoy2 = 1/(1-p22_optim);

%construction des vecteurs des explicatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X1 = lagmatrix(Yt,1);
X2 = lagmatrix(Yt,2);
X3 = lagmatrix(Yt,3);
X4 = lagmatrix(Yt,4);

%estimation des Yt dans le régime 1 et 2 + calcul de leurs caractéristiques
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Régime 1

Y_est_1 = c_optim(1,1) + phi1_optim*X1 + phi2_optim*X2 + phi3_optim*X3 + phi4_optim*X4;
Y_est_1 = Y_est_1(5:end,1);

Y_moy_1 = mean(Y_est_1); %moyenne
Y_var_1 = var(Y_est_1); %variance
Y_max_1 = max(Y_est_1); %valeur maximum
Y_min_1 = min(Y_est_1); %valeur minimum

%Régime 2

Y_est_2 = c_optim(2,1) + phi1_optim*X1 + phi2_optim*X2 + phi3_optim*X3 + phi4_optim*X4;
Y_est_2 = Y_est_2(5:end,1);

Y_moy_2 = mean(Y_est_2); %moyenne
Y_var_2 = var(Y_est_2); %variance
Y_max_2 = max(Y_est_2); %maximum
Y_min_2 = min(Y_est_2); %minimum


%Comparaison des caractéristiques des deux modèles 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

names = [string('Moyenne');string('Variance');string('Maximum');string('Minimum');];

format compact ;
disp('--------------------------------------------------------------');
disp('CARACTERISTIQUES DES MODELES');
disp('--------------------------------------------------------------');
disp('--- Valeurs       Régime 1        Régime 2 -------------------');
disp('--------------------------------------------------------------');
disp([names(1,1) Y_moy_1 Y_moy_2]);
disp([names(2,1) Y_var_1 Y_var_2]);
disp([names(3,1) Y_max_1 Y_max_2]);
disp([names(4,1) Y_min_1 Y_min_2]);
disp('--------------------------------------------------------------');
disp('--------------------------------------------------------------');


%% Question 3 : Comparer graphiquement les dates de récession obtenues au dates fournies par le NBER

%importation des dates de récessions fournies par le NBER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RecessDate = xlsread('Hamilton.xlsx', 'NBER');
RecessDate = RecessDate(:,2);

%comparaisons des "vrais" dates de récessions fournies par le NBER (dates
%du régime 2), avec les dates estimées du régime 2 du modèle MS-AR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on;
plot(dates, PL_hamilton, '--r');
plot(dates, RecessDate,'b');
title('Comparaison des états estimés du modèle et des "vrais" états du NBER');
datetick('x','mm-yyyy', 'keepticks');
legend('P[St=2|I(T)]', 'vrais états');
hold off;


%% Question 4 : Vérfier que les résidus du modèle satisfont les conditions usuelles 

%Tests de diagnostics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resGen = Res_Gen_MSAR(hamilton_theta, PF_hamilton);

%test d'autocorrélation

%méthode 1
[lb_Hamilton, pVal_lb_ham, stat_lb_ham, critikVal_lb_ham] = lbqtest(resGen, 'lags', 20, 'alpha', 0.05);

%méthode 2
[lb_Hamilton_2, pVal_lb_ham_2, stat_lb_ham_2, critikVal_lb_ham_2] = MyLbTest(resGen, k);

%effet ARCH

%méthode 1
[ARCH_Hamilton, pVal_ARCH_ham, stat_ARCH_ham, critikVal_ARCH_ham] = lbqtest(resGen.^2, 'lags', 20, 'alpha', 0.05);

%méthode 2
[ARCH_Hamilton_2, pVal_ARCH_ham_2, stat_ARCH_ham_2, critikVal_ARCH_ham_2] = MyLbTest(resGen.^2, k);

%normalité

[testNorm_Hamilton,testNorm_Hamilton_Pval, testNorm_Hamilton_StatTest, testNorm_Hamilton_CritikVal]  = jbtest(resGen);

%Test de significativité des coefficients

[test_student_hamilton, student_hamilton_pval, student_hamilton_critikVal] = MyStudentTest(hamilton_t_stat,0.05);


%Affichage des résultats des tests de diagnostics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('--------------------------------------------------------------');
disp('TESTS DE DIAGNOSTIC DU MODELE');
disp('Les résidus sont non auto-corrélés :');
if lb_Hamilton == 0
    disp('Vrai !');
else 
    disp('Faux !');
end
disp('Les résidus ne présentent pas d effets ARCH :');
if ARCH_Hamilton == 0
    disp('Vrai !');
else 
    disp('Faux !');
end
disp('Les résidus peuvent être considérés comme normaux :');
if testNorm_Hamilton == 0
    disp('Vrai !');
else 
    disp('Faux !');
end
disp ('Les coefficients sont significatifs : ');
if test_student_hamilton == 1
    disp('Vrai !');
else 
    disp('Faux !');
end
disp('--------------------------------------------------------------');
disp('--------------------------------------------------------------');


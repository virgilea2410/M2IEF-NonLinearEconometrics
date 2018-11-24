function [modeleChoisi, statsADF, statsPP] = TestStationnarite(y)
%%                                  Mise en place

y_diff = diff(y);

n = size(y,1) - 1;

x1 = lagmatrix(y,1);
x1_diff = lagmatrix(y_diff,1); 
y = y(3:end,:);
x1 = x1(3:end,:);
y_diff = y_diff(2:end,:);
x1_diff = x1_diff(2:end,:);


%%                                     Modèle 3 

%Matrice des regresseurs
X = [ones(size(x1)) x1 [(1:size(x1))]' x1_diff ];

%Regression
[beta_3, ~, EResidus3, ~, EStats3] = regress(y_diff, X);

%Variance des résidus
residus_3 = y_diff - X*beta_3;
SCR_3 = sum(residus_3'*residus_3);
sigma_residus_3 = SCR_3/(n-4);

%Variance des estimateurs
sigma_beta_3 = sigma_residus_3 * inv(X'*X);

%stat de student de beta
t_student_3 =  beta_3(2)/sqrt(sigma_beta_3(2,2));

%zone de rejet
rejet_3 = tinv(0.975,n);
pvalue_3 = tcdf(t_student_3, n);

%stat de student de b
t_student_3_bis =  beta_3(3)/sqrt(sigma_beta_3(3,3));

%%                                  Modèle 2

%Matrice des regresseurs
X = [ones(size(x1)) x1 x1_diff];

%Regression
beta_2 = regress(y_diff,X);

%Variance residuelle
residus_2 = y_diff - X*beta_2;
SCR_2 = sum(residus_2'*residus_2);
sigma_residus_2 = SCR_2/(n-3);

%Variance des coefficients
sigma_beta_2 = sigma_residus_2 * inv(X'*X);

%stat de student pour x(t-1)
t_student_2 = beta_2(2)/sqrt(sigma_beta_2(2,2));

%zone de rejet
rejet_2 = tinv(0.975,n);
pvalue_2 = tcdf(t_student_2, n);

%stat de student pour la constante
t_student_2_bis = beta_2(1)/sqrt(sigma_beta_2(1,1));

%%                                  Modèle 1  :

%Matrice des regresseurs
X = [x1 x1_diff];

%Regression
beta_1 = regress(y_diff,X);

%Variance residuelle
residus_1 = y_diff - X*beta_1;
SCR_1 = sum(residus_1'*residus_1);
sigma_residus_1 = SCR_1/(n-2);

%Variance des estimateurs
sigma_beta_1 = sigma_residus_1 * inv(X'*X);

%Stat de student pour x(t-1)
t_student_1 = beta_1(1,1)/sqrt(sigma_beta_1(1,1));

%stat de student pour la constante
t_student_1_bis = beta_1(2,1)/sqrt(sigma_beta_1(2,2));

%Zone de rejet
rejet_1 = tinv(0.975,n);
pvalue_1 = tcdf(t_student_1, n);

%%                          Verification des resultats

%Modèle 3
[h3, ADFPval3, tstat3, critikVal3]  = adftest(y, 'model', 'TS', 'lags', 1);

%Modèle 2
[h2, ADFPval2, tstat2, critikVal2] = adftest(y, 'model', 'ARD', 'lags', 1);

%Modèle 1
[h1, ADFPval1, tstat1, critikVal1] = adftest(y, 'model', 'AR', 'lags', 1);

statsADF = [ADFPval3 ADFPval2 ADFPval1; tstat3 tstat2 tstat1; t_student_3_bis t_student_2_bis t_student_1_bis ];

%%                          Test complementaires

% Test de Philippe Perron

%Modèle 3
[p3, PPPval3, ppstat3, PPcritikVal3] = pptest(y, 'model', 'TS', 'lags', 1);

%Modèle 2
[p2, PPPval2, ppstat2, PPcritikVal2] = pptest(y, 'model', 'ARD', 'lags', 1);

%Modèle 1
[p1, PPPval1, ppstat1, PPcritikVal1] = pptest(y, 'model', 'AR', 'lags', 1);

statsPP = [PPPval3 PPPval2 PPPval1; ppstat3 ppstat2 ppstat1];

%%                      Renvoi du modèle choisi

if ((t_student_3 < - 3.45) & (abs(t_student_3_bis) > 1.96))
    modeleChoisi = 3;
elseif ((t_student_2 < - 2.80) & (abs(t_student_2_bis) > 1.96))
    modeleChoisi = 2;
elseif (t_student_1 < - 1.96)
    modeleChoisi = 1;
else 
    modeleChoisi = 0;
end 

end
function [MSAR_test, StatDuTest, pcritique, critikVal] = myMSARTest(P)
    %fonction permettant d'effectuer le test d'un modèle non linéaire MS-AR
    %contre le modèle AR linéaire correspondant, sur la série Yt.
    %Renvoi 1 si la linéarité du modèle ne peut pas être rejeté (non rejet
    %de H0) : Modèle AR
    %Renvoi 0 si le modèle est MS-AR (rejet de H0)

global choice;
global Yt;
warning('off');

n = size(Yt,1);

%estimation d'un modèle AR(1)

%méthode 1 

[Test_AR_phi, Test_AR_t_stat, Test_AR_var_res, Test_AR_resids, Test_AR_aic, Test_AR_logvraiss, Test_AR_std_coeff] = AR_LogV(1);

%méthode 2
choice = 1;
[Test_AR_phi_2, Test_AR_t_stat_2, Test_AR_var_res_2, Test_AR_SCR, Test_AR_aic_2, Test_AR_logvraiss_2] = AR_LogV_2(1);

%estimation d'un modèle MS-AR(1,1)
p11 = P(1,1);
p22 = P(2,2);
c1_test = Test_AR_phi(1,1);
c2_test = Test_AR_phi(1,1);
phi1_test = Test_AR_phi(2,1);
phi2_test = Test_AR_phi(2,1);
sig1_test = sqrt(Test_AR_var_res);
sig2_test = sqrt(Test_AR_var_res);
    
x0_test = [norminv(p11) norminv(p22) c1_test c2_test phi1_test phi2_test sig1_test sig2_test]';
x0_test = x0_test + randn(8,1)./100;

choice = 1;
[Test_MSAR_theta, Test_MSAR_std_coeff, Test_MSAR_t_stat, Test_MSAR_logvraiss] = maxVraissemblance(x0_test);

StatDuTest = 2*(Test_MSAR_logvraiss - Test_AR_logvraiss);

%simulation pour obtenir le seuil de rejet de la stat du test
moy = mean(Yt);
etype = std(Yt);

%initialisation des paramètres de l'estimation de la probabilité critique
y_simul = zeros(n,10);
Test_ARsimul_phi = zeros(2,10);
Test_ARsimul_t_stat = zeros(2,10);
Test_ARsimul_var_res = zeros(1,10);
Test_ARsimul_resids = zeros(n,10);
Test_ARsimul_aic = zeros(1,10);
Test_ARsimul_logvraiss = zeros(1,10);
c1_simul = zeros(1,10);
c2_simul = zeros(1,10);
phi1_simul = zeros(1,10);
phi2_simul = zeros(1,10);
sig1_simul = zeros(1,10);
sig2_simul = zeros(1,10);
x0_simul = zeros(8,10);
Test_MSARsimul_theta = zeros(8,10);
Test_MSARsimul_t_stat = zeros(8,10);
Test_MSARsimul_logvraiss = zeros(1,10);
StatDuTestSimul = zeros(1,10);

for i=1:10
    y_simul(:,i) = normrnd(moy,etype, n, 1);
    
    choice = 1;
    Yt = y_simul(:,i);
  
    [Test_ARsimul_phi(:,i), Test_ARsimul_t_stat(:,i), Test_ARsimul_var_res(1,i), Test_ARsimul_resids(:,i), Test_ARsimul_aic(1,i), Test_ARsimul_logvraiss(1,i)] = AR_LogV(1);
    
    c1_simul(1,i) = Test_ARsimul_phi(1,i);
    c2_simul(1,i) = Test_ARsimul_phi(1,i);
    phi1_simul(1,i) = Test_ARsimul_phi(2,i);
    phi2_simul(1,i) = Test_ARsimul_phi(2,i);
    sig1_simul(1,i) = sqrt(Test_ARsimul_var_res(1,i));
    sig2_simul(1,i) = sqrt(Test_ARsimul_var_res(1,i));
    
    x0_simul(:,i) = [norminv(p11) norminv(p22) c1_simul(1,i) c2_simul(1,i) phi1_simul(1,i) phi2_simul(1,i) sig1_simul(1,i) sig2_simul(1,i)]';
    %x0_simul(:,i) = x0_simul(1,i) + randn(8,1)./100;
    
    [Test_MSARsimul_theta(:,i), ~, Test_MSARsimul_t_stat(:,i), Test_MSARsimul_logvraiss(1,i)] = maxVraissemblance(x0_simul(:,i));
    
    StatDuTestSimul(1,i) = 2*(Test_MSARsimul_logvraiss(1,i) - Test_ARsimul_logvraiss(1,i));
end

pcritique = StatDuTestSimul > StatDuTest;
pcritique = sum(pcritique)/size(pcritique,2);

critikVal = mean(StatDuTestSimul);

MSAR_test = pcritique > 0.05;

end
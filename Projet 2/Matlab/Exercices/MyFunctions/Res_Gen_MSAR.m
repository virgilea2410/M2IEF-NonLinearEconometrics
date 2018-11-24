function [resids_gen_MSAR, Y_info_1, Y_info_2] = Res_Gen_MSAR(MSAR_theta, PP)
    % Fonction permettant de calculer les R�sidus G�n�ralis�s d'un mod�le, en
    % fonction d'un vecteur de coefficients estim�s optimaux, MSAR_theta, et d'un vecteur
    % de Probabilit�s filtr�es PP, ind�x� par le temps et tel que
    % PP = P[St=1 | I(t)]
    %
    %La Fonction calcul aussi les caract�ristiques propres au deux r�gimes :
    %moyennes, variance, maximum et minimum du Y dans chacun des r�gimes
    %
    %choice = 1 --> Y(t) = C(St) + phi1(St) + sigma_res(St)
    %choice = 2 --> Y(t) = C(St) + phi1 + phi2 + sigma_res(St)
    %choice = 3 --> Y(t) = C(St) + phi1 + phi2 + phi3 + phi4 + sigma_res

global Yt;
global choice;
warning('off');

if choice == 2
    beg = 3;
elseif choice == 3
    beg = 5;
end

%d�finition des coefficients optimaux du mod�le
%en fonction de la variable gloable choice
MSAR_cste = MSAR_theta(3:4,1);

if choice == 2
    MSAR_phi1 = MSAR_theta(5,1);
    MSAR_phi2 = MSAR_theta(6,1);
    MSAR_std_res = MSAR_theta(7:8,1);
elseif choice == 3
    MSAR_phi1 = MSAR_theta(5,1);
    MSAR_phi2 = MSAR_theta(6,1);
    MSAR_phi3 = MSAR_theta(7,1);
    MSAR_phi4 = MSAR_theta(8,1);
    MSAR_std_res = MSAR_theta(9,1);
end

%initialisation des variables explicatives
X1 = lagmatrix(Yt,1);
X2 = lagmatrix(Yt,2);
X3 = lagmatrix(Yt,3);
X4 = lagmatrix(Yt,4);

X1 = X1(beg:end,1);
X2 = X2(beg:end,1);
X3 = X3(beg:end,1);
X4 = X4(beg:end,1);
Y = Yt(beg:end,1);

%calcul des Y estim� conditionnellement � chaque r�gime, puis calcul des
%r�sidus conditionnels � chaque r�gime, divis� par leurs ecart types
%r�siduels respectifs
if choice == 2
    Yt_est_1 = ( MSAR_cste(1,1) + MSAR_phi1 * X1 + MSAR_phi2 * X2); 
    resids_MSAR_1 = (Y - Yt_est_1) ./ MSAR_std_res(1,1);

    Yt_est_2 = ( MSAR_cste(2,1) + MSAR_phi1 * X1 + MSAR_phi2 * X2 );
    resids_MSAR_2 = (Y - Yt_est_2) ./MSAR_std_res(2,1);
elseif choice == 3
    Yt_est_1 = ( MSAR_cste(1,1) + MSAR_phi1 * X1 + MSAR_phi2 * X2+ MSAR_phi3 * X3 + MSAR_phi4 * X4); 
    resids_MSAR_1 = (Y - Yt_est_1) ./ MSAR_std_res;

    Yt_est_2 = ( MSAR_cste(2,1) + MSAR_phi1 * X1 + MSAR_phi2 * X2 + MSAR_phi3 * X3 + MSAR_phi4 * X4);
    resids_MSAR_2 = (Y - Yt_est_2) ./MSAR_std_res;
end

%calcul des r�sidus g�n�ralis�s en faisant la somme pond�r�s par leurs
%probabilit�s filtr�es des deux s�ries de r�sidus des 2 r�gimes
resids_gen_MSAR = resids_MSAR_1 .* PP(beg:end,1) + (resids_MSAR_2) .* (1-PP(beg:end,1));

%caract�ristiques du r�gime 1 r�gime 1
Y_info_1(1,1) = mean(Yt_est_1);
Y_info_1(2,1) = var(Yt_est_1);
Y_info_1(3,1) = max(Yt_est_1);
Y_info_1(4,1) = min(Yt_est_1);

%caract�ristiques du r�gime 1 r�gime 1
Y_info_2(1,1) = mean(Yt_est_2);
Y_info_2(2,1) = var(Yt_est_2);
Y_info_2(3,1) = max(Yt_est_2);
Y_info_2(4,1) = min(Yt_est_2);

end

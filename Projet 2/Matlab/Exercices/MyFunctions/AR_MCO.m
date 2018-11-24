function [phi, t_stat, SCR, aic] = AR_MCO(p)
    %Fonction permettant d'estimer un mod�le AR(p) par la m�thode des
    %moindres carr�s ordinaires

global Yt;
warning('off');

Y = Yt;

Y = y(p:end,1);

X = (lagmatrix(y,(1:p)));
X = X(p:end,:);

AR_model = fitlm(X,Y);

phi = AR_model.Coefficients.Estimate;

t_stat = AR_model.Coefficients.tStat;

end
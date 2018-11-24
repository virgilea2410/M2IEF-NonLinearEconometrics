function [phi, t_stat, SCR, aic] = AR_MCO(p)
    %Fonction permettant d'estimer un modèle AR(p) par la méthode des
    %moindres carrés ordinaires

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
function [y,s] = Simul_MSAR(theta,P,T)
    %Fonction permettant de simuler un processus MS-AR à probabilités de
    %transition variables, à l'aide d'une chaine de markov à probabilités
    %de transition variables.
    %la simulation a propremement parler et ensuite la même que pour un
    %processus de markov simulé à probabilité de transition fixes

y = zeros(T,1); 

s = Simul_Markov_PVar(P,T);

for t = 2:T
    eps = randn(1,1);
    y(t) = (theta(1) + theta(3)*y(t-1) + theta(4)*eps)*(s(t)==1) + ...
           (theta(2) + theta(3)*y(t-1) + theta(5)*eps)*(s(t)==2);
end
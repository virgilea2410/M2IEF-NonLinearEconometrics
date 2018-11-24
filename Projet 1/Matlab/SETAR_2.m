function [ phi, SCR, sigma ] = SETAR_2(z,y,p,d,c )

T = size(y,1);

%Définition du vecteur des variables à expliquer Y
Y = z((max(p,d)+1):end,:);

nobs = T - max(p,d);

X = zeros(nobs, 2*(p+1));

%definition de la variable de changement de seuil Pour x(t), la variable de changement de seuil est x(t-1), on décale donc notre matrice de 1
Q = lagmatrix(Y,d);

% création de la variable indicatrice : c'est un vecteur contenant des 1 quand x(t-1) < C, 0 sinon. 
I = Q<c;

% création des variables lagés 1 jusqu'a p (vecteur des variables explicatives X)
X1 = lagmatrix(Y,(1:p));

%assemblage de la matrice des regresseurs
X = [ones(nobs,1).*I(:,1) X1.*repmat(I(:,1),1,p) ones(nobs,1).*(1-I(:,1)) X1.*(1-repmat(I(:,1),1,p))];

% estimation des coefficient par MCO
phi = regress(Y,X);

%calcul du carré des résidus 
U = Y - X*phi;
U = U(4:end,:);

SCR = U'*U;

sigma = SCR/(nobs-2*(p+1));

%Calcul de la stats de significativé de Student pour tester la
%significativité des coefficients 

var_covar = inv(Y'*Y) *sigma;
sigma_B = diag(var_covar);
t_stat = phi./sqrt(sigma_B);

end


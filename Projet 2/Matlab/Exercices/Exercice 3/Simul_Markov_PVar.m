function s = Simul_Markov_PVar(Ptrans,T)
    %Fonction permettant de simuler une chaine de Markov à probabilités
    %de transition variables : Ptrans n'est plus une matrice 2x2 mais 2
    %vecteurs 1xT : Ptrans est de taille 2xT

s = zeros(T,1);
s(1,1) = 1; % initialisation : état 1 en t=1

for t=2:T 
    
  u = rand(1,1);
  
  % si état 1 en t, état en t+1
  if s(t-1)==1 
      if u <= Ptrans(t,1)
        s(t) = 1;
      else
        s(t) = 2;
      end
      
  % si état 2 en t, état en t+1
  elseif s(t-1)==2
      if u <= Ptrans(t,2)
        s(t) = 2;
      else
        s(t) = 1;
      end
  end
    
end


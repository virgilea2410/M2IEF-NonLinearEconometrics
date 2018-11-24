function [LogV]= vraissemblanceAR(x0)
    %fonction de vraisemblance d'un modèle AR
    %
    %choice = 1 --> Y(t) = AR(1)
    %choice = 2 --> Y(t) = AR(2)
    %choice = 3 --> Y(t) = AR(4)

global choice;
global Yt;
warning('off');

if choice == 1
    c= x0(1,1);
    phi1 = x0(2,1);       
    std_res = x0(3,1);
elseif choice == 2
    c= x0(1,1);
    phi1 = x0(2,1);
    phi2 = x0(3,1);
    std_res = x0(4,1);
elseif choice == 3
    c= x0(1,1);
    phi1 = x0(2,1);
    phi2 = x0(3,1);
    phi3 = x0(4,1);
    phi4 = x0(5,1);
    std_res = x0(6,1);
end

T = size(Yt,1);
LogV=0.0; 

if choice == 1
   i = 2;
elseif choice == 2
   i = 3;
elseif choice == 3
   i = 5;
end

for J_Iter =i:T  
    if choice == 1
        F_CAST = [Yt(J_Iter,1) - c - phi1*Yt(J_Iter -1,1)]; %AR(1)
    elseif choice == 2
        F_CAST = [Yt(J_Iter,1) - c - phi1*Yt(J_Iter -1,1) - phi2*Yt(J_Iter -2,1)]; %AR(2)
    elseif choice == 3
        F_CAST = [Yt(J_Iter,1) - c - phi1*Yt(J_Iter -1,1) - phi2*Yt(J_Iter -2,1) - phi3*Yt(J_Iter -3,1) - phi4*Yt(J_Iter -4,1)]; %AR(4)
    end
 
    PR_VAL = (1./sqrt(2*pi*std_res)) .* exp(-0.5*F_CAST.^2./std_res);
    %PR_VAL= (1./sqrt(2*pi*(std_res.^2))) .* exp(-0.5*F_CAST.^2./(std_res.^2)) .* PROB_DD;  % f(yt,St|It-1) (étape 3)
    %PR_VAL= (1./sqrt(2*pi*(std_res.^2))) .* exp(-0.5*F_CAST.^2./(std_res.^4)) .* PROB_DD;  % f(yt,St|It-1) (étape 3)
    
    LIK=-1*log(PR_VAL); 

    LogV = LogV + LIK;
end

end
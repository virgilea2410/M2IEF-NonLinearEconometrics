function [d_optim, phi_optim, var_res_optim, decision, pvalue, c_optim] = SETAR_FIT_D(y,p,c)

for d=1:5
    c_optim(d,1) = SETAR_FIT(y, p, d);
    
    [phi(:,d), ~ , var_res(d,1)] = SETAR(y,p,d,c_optim(d,1));

    [decision(d,1), pvalue(d,1)] = TestSETAR(y,p,d,c_optim(d,1));
end

[~, index] = min(var_res);

d_optim = index;
c_optim = c_optim(index,1);

phi_optim = phi(:,index);
var_res_optim = var_res(index,1);

end





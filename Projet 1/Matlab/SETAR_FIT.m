function [ c_optim, phi_optim ] = SETAR_FIT(y, p, d)

c = y;
c = sort(c);

nobs = size(c,1);

upDown = floor(0.15 * nobs);

c = c(upDown:(end-upDown),:);

result = ones(size(c,1),1);

for i = 1:size(c,1)
    [ ~, ~, result(i,1)] = SETAR(y, p, d, c(i,1));
end

[minimum, index] = min(result);

c_optim = c(index,1);

phi_optim = ones(8,1);
[ phi_optim, SCR_optim, ~] = SETAR(y, p, d, c_optim);

varRes_optim = SCR_optim/nobs;

varcovar = varRes_optim .* inv(y'*y);

end

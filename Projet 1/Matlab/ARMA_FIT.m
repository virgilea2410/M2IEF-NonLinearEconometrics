function [p_optim, q_optim] = ARMA_FIT(y)

data = y;

p_max = 4;
q_max = 4;

for i=0:p_max
    for j=0:q_max
        model = arima(i,0,j);
        [~, ~, logV] = estimate(model,data, 'print', false);
        est_bic(i+1,j+1) = aicbic(logV,1+(i+j),size(data,1));
    end
end

[v_min, row_min] = min(est_bic);
[~, col_min] = min(v_min);

p_optim = row_min(1, col_min) - 1;
q_optim = col_min - 1;

end

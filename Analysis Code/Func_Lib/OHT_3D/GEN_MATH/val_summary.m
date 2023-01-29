function [resp_summary] = val_summary(values)

if size(values,2) > 1
    error('Values should be an n x 1 vector');
end
n = size(values,1);

nanct = sum(isnan(values));

resp_opt = unique(values);
if nanct > 0;
    resp_opt = [resp_opt; NaN];
end
m = size(resp_opt,1);

resp_summary = zeros(m,2);
for i = 1:1:m
    resp_summary(i,1) = resp_opt(i);
    resp_summary(i,2) = sum(values == resp_opt(i));
    if (i == m) && (nanct > 0)
        resp_summary(i,2) = nanct;
    end
end
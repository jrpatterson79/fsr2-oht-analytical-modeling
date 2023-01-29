function [s_update] = KF_Ensemble(s_ens,model,y,R)

num_unk = size(s_ens,1);
num_ens = size(s_ens,2);
num_data = size(y,1);

hs_ens = zeros(num_data,num_ens);
for i = 1:1:num_ens
    hs_ens(:,i) = feval(model,s_ens(:,i));
end

HQHt_approx = cov(hs_ens');

d_cov = HQHt_approx + R;

mu_hs_ens = mean(hs_ens,2);
mu_s = mean(s_ens,2);

s_resids = s_ens - repmat(mu_s,1,num_ens);
hs_resids = hs_ens - repmat(mu_hs_ens,1,num_ens);

QHt_approx = zeros(num_unk,num_data);
for i = 1:1:num_unk
    for j = 1:1:num_data
        QHt_approx(i,j) = 1./(num_ens-1)*(s_resids(i,:)*(hs_resids(j,:)'));
    end
end

dfit_resids = repmat(y,1,num_ens) - hs_ens;

s_update = zeros(num_unk,num_ens);
for i = 1:1:num_ens
    s_update(:,i) = s_ens(:,i) + QHt_approx*(d_cov\dfit_resids(:,i));
end
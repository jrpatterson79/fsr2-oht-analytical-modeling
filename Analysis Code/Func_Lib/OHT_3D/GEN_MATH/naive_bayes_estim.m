function [cat_estimate] = naive_bayes_estim(features, cat_distr_assumed, cat_distr_params)

num_cats = size(cat_distr_assumed,1);
num_features = size(cat_distr_assumed,2);

if size(features,2) ~= num_features
    error('Features array does not appear to have the correct number of values');
end
if size(cat_distr_params,2) ~= num_features
    error('cat_distr_params array does not appear to have the correct number of values');
end

num_datapoints = size(features,1);

% Go through each, given classifier figure out which category it should
% belong to
negloglike_mat = zeros(num_datapoints, num_cats);
for i = 1:1:num_datapoints
    for j = 1:1:num_cats
        for k = 1:1:num_features
            if(~isnan(features(i,k)))
                if (cat_distr_assumed(j,k) == 1)
                    mean = cat_distr_params(j,k,1); var = cat_distr_params(j,k,2);
                    negloglike_mat(i,j) = negloglike_mat(i,j) + .5*log(var) + .5*(features(i,k)-mean).^2./var;
                elseif (cat_distr_assumed(j,k) == 2)
                    meanlog = cat_distr_params(j,k,1); varlog = cat_distr_params(j,k,2);
                    negloglike_mat(i,j) = negloglike_mat(i,j) + .5*log(varlog) + log(features(i,k)) + .5*(log(features(i,k))-meanlog).^2./varlog;
                elseif (cat_distr_assumed(j,k) == 3)
                    lambda = cat_distr_params(j,k,1);
                    negloglike_mat(i,j) = negloglike_mat(i,j) - log(lambda) + (lambda*features(i,k));
                end
            end
        end
    end
end
        
best_negloglike = (min(negloglike_mat'))';
cat_estimate = zeros(num_datapoints,1);
for i = 1:1:num_datapoints
    index_best = find(negloglike_mat(i,:) == best_negloglike(i));
    cat_estimate(i) = index_best;
end

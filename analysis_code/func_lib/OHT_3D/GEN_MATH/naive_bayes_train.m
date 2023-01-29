function [cat_distr_params] = naive_bayes_train(features,categories,cat_distr_assumed)

%naive_bayes_train: Function which computes distributional parameters
%during training of a Naive Bayes Classifier. Once trained, the NBC can be
%used in order to guess the class of a set of un-classified datapoints,
%given their feature values. If any datapoint does not have a known value
%for a particular feature, assign NaN to that datapoint's feature in the
%features matrix.
%
%Syntax:
%   [cat_distr_params] = naive_bayes_train(features,categories,cat_distrs)
%where:
%   -cat_distr_params: a (numcategories x numfeatures x 2) tensor
%   where numcategories is the number of categories to which a datapoint
%   can be classified, numfeatures is the number of features used to
%   estimate the category / class, and the third dimension of the tensor
%   corresponds to the various parameters of the distributions to be
%   estimated
%   -features: a (numdatapoints x numfeatures) array giving the values for
%   each feature for each datapoint.
%   -categories: a (numdatapoints x 1) vector giving the category to which
%   each datapoint is known to be assigned.
%   -cat_distrs: a (numcategories x numfeatures) matrix giving the assumed
%   distribution of each feature value within each category. Possible
%   values:
%       1 = Gaussian distribution, cat_distr_params(i,j,[1 2]) will contain
%       the mean and variance of the distribution
%       2 = log-normal distribution, cat_distr_params(i,j,[1 2]) will
%       contain the mean and variance of the distribution
%       3 = exponential distribution, cat_distr_params(i,j,[1]) will
%       contain the rate parameter (lambda)

% NOTE: More checking / debugging is needed on this. Need to check for
% cases where no datapoints belong to a particular class, among others.

num_cats = size(cat_distr_assumed,1);
num_datapoints = size(features,1);
num_features = size(features,2);

%Cases for distributions (in cat_distr_assumed variable)
%1) Guassian - cat_distr_params = mean,var
%2) log-normal - cat_distr_params = mean,var
%3) exponential - cat_distr_params = rate
cat_distr_params = zeros(num_cats,num_features,2);

for i = 1:1:num_cats
    index_list = find(categories == i);
    for j = 1:1:num_features
        if (cat_distr_assumed(i,j) == 1)
            cat_distr_params(i,j,1) = nanmean(features(index_list,j));
            cat_distr_params(i,j,2) = nanvar(features(index_list,j));
            if cat_distr_params(i,j,2) == 0
                cat_distr_params(i,j,2) = 1e6;
            end
        elseif (cat_distr_assumed(i,j) == 2)
            cat_distr_params(i,j,1) = nanmean(log(features(index_list,j)));
            cat_distr_params(i,j,2) = nanvar(log(features(index_list,j)));
            if cat_distr_params(i,j,2) == 0
                cat_distr_params(i,j,2) = 1e6;
            end
        elseif (cat_distr_assumed(i,j) == 3)
            cat_distr_params(i,j,1) = 1/nanmean(features(index_list,j));
        end
    end
end

function [relz, varargout] = mh_sampling(nlogP_func,init_params,covar_mat,varargin)

%mh_sampling: Simplified Metropolis-Hastings sampling which uses a given
%covariance matrix to generate Gaussian-based steps between candidates.
%Includes burn-in period. Takes advantage of fact that probability of steps
%is the same forwards and backwards, so does not include second term in
%general Metropolis-Hastings algorithm
%[relz, [nl_probs, num_accept]] = mh_sampling(nlogP_func,init_params,[num_relz, num_toss, visual])
%where:
%   -relz is an m x num_relz output matrix where m is the number of
%   unknowns and num_relz is the number of generated realizations
%   -nl_probs (optional) is the negative-log-likelihood of each realization
%   (1 x num_relz)
%   -num_accept (optional) is the number of NEW markov chain steps that
%   were accepted during generation of the chain
%   -nlogP_func is the function giving the negative log-probability of the
%   parameters given the data, within a constant (created using anonymous @
%   function syntax). nlogP_func should have one input (the parameter
%   vector), and one output (the negative log-probability value, within a
%   constant).
%   -init_params is an initial guess at the parameter vector. init_params
%   should be a column vector of size m x 1
%   -covar_mat is the covariance matrix used to generate new candidate
%   realizations. (m x m matrix)
%   -num_relz (optional) is the number of realizations to generate
%   -num_toss (optional) is the number of steps used to initialize the
%   Markov Chain.
%   -visual (optional) is a 1/0 variable indicating whether visualizations
%   of the chain's mean and variance should be shown to the user. (Off by
%   default)

%Initial error checking
m = numel(init_params);
if size(covar_mat,1)~= size(covar_mat,2) || numel(covar_mat) ~= (m^2)
    error(['The covariance matrix must be square with each dimension equal ', ...
        'to the number of parameters']);
end

%Calculate the cholesky decomposition of the covariance matrix - used for
%generating new steps.
sd_mat = (chol(covar_mat))';

%Default values for function parameters
num_relz = 10000*m;
num_toss = 1000*m;
visual = 0;

%Assign additional arguments, if supplied by user
[num_relz, num_toss, visual] = process_extra_args(varargin,num_relz,num_toss, visual);

%Start the Markov chain & perform burn-in. Do nothing with the accepted
%realizations.
curr_cand = init_params;
curr_neglogprob = nlogP_func(curr_cand);
for i = 1:1:num_toss
    new_cand = curr_cand + sd_mat*randn(m,1);
    neglogalpha = -log(rand);
    new_neglogprob = nlogP_func(new_cand);
    if (new_neglogprob < curr_neglogprob) || (neglogalpha > (new_neglogprob - curr_neglogprob))
        curr_cand = new_cand;
        curr_neglogprob = new_neglogprob;
    end
    disp(['Burn in, # ', num2str(i)]);
end

%Begin the sampling. Initialize the realizations matrix, then use MCMC to
%accept or reject new steps. Store all steps.
num_accept = 1;
relz = zeros(m,num_relz);
nl_probs = zeros(1,num_relz);
for i = 1:1:num_relz
    new_cand = curr_cand + sd_mat*randn(m,1);
    neglogalpha = -log(rand);
    new_neglogprob = nlogP_func(new_cand);
    if (new_neglogprob < curr_neglogprob) || (neglogalpha > (new_neglogprob - curr_neglogprob))
        curr_cand = new_cand;
        curr_neglogprob = new_neglogprob;
        num_accept = num_accept + 1;
    end
    relz(:,i) = curr_cand;
    nl_probs(i) = curr_neglogprob;
    disp(['Run # ', num2str(i),', Accepted ', num2str((num_accept./(i+1))*100), '%']);
    %Display visualizations if requested
    if visual == true
        hold on
        figcols = floor(m^.5);
        figrows = ceil(m/figcols);
        figure(1)
        for j = 1:1:m
            subplot(figrows,figcols,j)
            hold on
            plot(i,mean(relz(j,1:i)),'o')
            hold off
        end
        figure(2)
        for j = 1:1:m
            subplot(figrows,figcols,j)
            hold on
            plot(i,var(relz(j,1:i)),'o')
            hold off
        end
    end
end

if nargout > 1
    varargout{1} = nl_probs;
end

if nargout > 2
    varargout{2} = num_accept;
end
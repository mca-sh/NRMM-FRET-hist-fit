function [prm,Nmax] = get_inference_param(inference_method)

switch inference_method
    case {'EM','EM-GMM'}
        Nmax = Inf;  % max. nb of observation (for speed)
        prm.niter = 600; % nb. of EM iterations
        prm.ninit = 500; % nb. of EM initializations
        prm.Kmax = 10;   % max. nb. of clusters
        prm.nSpl = 100;  % nb. of bootstrap samples
        prm.nsig = 3;    % nb. of sigmas for confidence intervals
        prm.tol = 1E-6;  % convergence criterion of log-likelihood ratio

    case 'DPMM'
        Nmax = 5000;                     % max. nb of observation (for speed)
        prm.niter = 4000;                % nb of iterations in Gibbs sampling
        prm.burnin = round(prm.niter/2); % nb. of burn in iterations
        prm.nruns = 5;                   % nb. of DPMM runs
        prm.K_init = 10;                 % initial nb. of clusters
        prm.clvl = 0.95;                 % credibility level for interval calculations
        prm.prior_flx = 'uniform';       % prior p(flx), 'uniform' or 'normal'
        prm.prior_E = 'uniform';         % prior p(E), 'uniform' or 'beta'
        prm.hyp.alpha_beta = 0.5;        % hyperparameter alpha for E's beta prior
        prm.hyp.beta_beta = 0.5;         % hyperparameter beta for E's beta prior
        
        prm.alpha0 = 0.003;                % alpha value of DP (fixed or initial value if sampled)
        % UNCOMMENT TO SAMPLE ALPHA
        % prm.hyp.alpha_a = 1;             % hyperparameter a for alpha's gamma prior 
        % prm.hyp.alpha_b = 0.1;           % hyperparameter b for alpha's gamma prior

    otherwise
        error('Input 1 must be ''EM'', ''EM-GMM'', or ''DPMM''.');
end
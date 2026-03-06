function [E_post,flx_post,w_post,E_spl,flx_spl,w_spl,plt,bestrun] = ...
    DPPM_discrete_fit(FRET,P_table,F_grid,E_grid,flx_grid,opt,plt)
% [E_post,flx_post,w_post,E_spl,flx_spl,w_spl,plt] = ...
%    DPPM_discrete_fit(FRET,P_table,F_grid,E_grid,flx_grid,opt,plt)

% collect inference parameters and constants
prm = opt.prm;
const = opt.const;

% Get grid dimensions
M = length(E_grid);
P = length(flx_grid);

% Build prior distribution p(E)
switch prm.prior_E
    case 'uniform' % uniform prior
        E_prior = ones(M,1); 
    case 'beta' % beta prior
        [prm.hyp.alpha_beta,prm.hyp.beta_beta] = fit_beta_prior(FRET);
        E_prior = generate_beta_prior(E_grid, prm.hyp.alpha_beta, ...
            prm.hyp.beta_beta)'; 
    otherwise
        error('prior_E must be ''uniform'' or ''beta''');
end

% Build prior distribution p(flx)
switch prm.prior_flx
    case 'uniform'
        flx_prior = ones(1,P); % uniform prior 
    case 'normal'
        flx_prior = normpdf(flx_grid,0.5,0.5); % Normal prior centered on 0.5
    otherwise
        error('prior_flx must be ''uniform'' or ''normal''');
end

% Compose prior
G0_prior = repmat(flx_prior/sum(flx_prior),[M,1]).*...
    repmat(E_prior/sum(E_prior),[1,P]);

% Normalize prior
G0_prior = G0_prior/sum(G0_prior,[1,2]);

% Plot prior
if ~isempty(plt)
    plot_prior(plt,E_grid,flx_grid,G0_prior);
end

% Print dummy line
[nb_max,~] = get_nb_char_to_print(prm);
fprintf(['\n',repmat(' ',1,nb_max)]);

% Parallel Gibbs sampler
K_trace = cell(1,prm.nruns);
LP_trace = cell(1,prm.nruns);
E_post = cell(1,prm.nruns);
flx_post = cell(1,prm.nruns);
w_post = cell(1,prm.nruns);
E_spl = cell(1,prm.nruns);
flx_spl = cell(1,prm.nruns);
w_spl = cell(1,prm.nruns);
LP = zeros(1,prm.nruns);
Kopt = zeros(1,prm.nruns);
parfor i = 1:prm.nruns
    [E_samples,flx_samples,w_samples,K_trace{i},LP_trace{i},~] = ...
        gibbs_sampler(FRET,P_table,E_grid,flx_grid,prm,prm.hyp,const,...
        G0_prior,F_grid,[],false,i);

    % Determine mode of posterior nb. of cluster
    Kopt(i) = round(mode(K_trace{i}((prm.burnin+1):prm.niter)));

    % Determine median of posterior probability
    LP(i) = median(LP_trace{i}((prm.burnin+1):prm.niter));
    
    % Sort parameters according to mean FRET efficiency
    [E_spl{i},colid] = sort(E_samples{Kopt(i)},2);
    matid = (colid-1)*size(E_spl{i},1) + ...
        repmat((1:size(E_spl{i},1))',1,Kopt(i));
    w_spl{i} = w_samples{Kopt(i)}(matid);
    flx_spl{i} = flx_samples{Kopt(i)}(matid);
    E_post{i} = median(E_spl{i},1);
    flx_post{i} = median(flx_spl{i},1);
    w_post{i} = median(w_spl{i},1);
    w_post{i} = w_post{i}/sum(w_post{i});
end
    
% Plot best out and traces
if opt.plotfinal
    for i = 1:prm.nruns
        plt{i}.objktr = plot_K_trace(plt{i}.objktr,K_trace{i});
        plt{i}.objpost = plot_logpost_trace(plt{i}.objpost,LP_trace{i});
        plt{i}.objkhist = plot_K_hist(plt{i}.objkhist,...
            K_trace{i}((prm.burnin+1):prm.niter));
        plt{i}.objposthist = plot_logpost_hist(plt{i}.objposthist,...
            LP_trace{i}((prm.burnin+1):prm.niter));

        PMF = calc_mixture_pmf(w_post{i},...
            {E_post{i},flx_post{i}},{E_grid,flx_grid},P_table);
        plt{i}.objdpmm0 = plot_mixture(plt{i}.objdpmm0,F_grid,PMF);
        plt{i}.objdpmmk = plot_mixture_components(plt{i}.objdpmmk,F_grid,...
            E_post{i},PMF);

        plt{i}.objEhist = plot_hist(plt{i}.objEhist,6,2,8,E_spl{i},'E');
        plt{i}.objflxhist = plot_hist(plt{i}.objflxhist,6,2,10,...
            flx2sig(flx_spl{i},const.k,const.S),'\\sigma r');
        plt{i}.objwhist = plot_hist(plt{i}.objwhist,6,2,12,w_spl{i},'w');
    end
    drawnow;
end

% Determine best run
[~,bestrun] = max(LP);
fprintf('Best run: %i\n',bestrun);
E_post = E_post{bestrun};
flx_post = flx_post{bestrun};
w_post = w_post{bestrun};
E_spl = E_spl{bestrun};
flx_spl = flx_spl{bestrun};
w_spl = w_spl{bestrun};
Kopt = Kopt(bestrun);

% Calculate credible intervals of posterior
dE_post = zeros(2,Kopt);
dflx_post = zeros(2,Kopt);
dw_post = zeros(2,Kopt);
for k = 1:Kopt
    [dE_k,~] = HPD_intervalle(E_spl(:,k)',prm.clvl);
    dE_post(:,k) = dE_k';
    [dflx_k,~] = HPD_intervalle(flx_spl(:,k)',prm.clvl);
    dflx_post(:,k) = dflx_k';
    [dw_k,~] = HPD_intervalle(w_spl(:,k)',prm.clvl);
    dw_post(:,k) = dw_k';
end

% Concatenate results
E_post = [E_post;dE_post];
flx_post = [flx_post;dflx_post];
w_post = [w_post;dw_post];


function [E_spl,flx_spl,w_spl,K_trace,LP_trace,plt] = gibbs_sampler(FRET,...
    P_table,E_grid,flx_grid,prm,hyp,const,G0_prior,F_grid,plt,plotiter,...
    runid)

% default
iv_iter = 50; % nb. of iterations to skip between each display

% get dimensions
N = length(FRET);
[L,M,P] = size(P_table);
MP = M*P;

% pre-calculations and vectorisation (for speed)
log_P_table = log(P_table + eps); 
log_P_flat = reshape(log_P_table, L, MP);
G0_flat = G0_prior(:);
log_G0_flat = log(G0_flat);
marginal_per_bin = reshape(P_table, L, MP) * G0_flat; 
log_marginal_per_bin = log(marginal_per_bin + eps);

% map FRET to grid indexes
[~,x_idx_all] = ismember(FRET,F_grid);
if any(x_idx_all==0)
    error('some FRET values lay outside of F_grid.');
end

% initialize clusters
K = prm.K_init; 
try
    % use k-means clustering to initialize cluster assignment (for 
    % reproducibility)
    [Z,E_kmeans] = kmeans(FRET(:),K,'Replicates',3,'MaxIter',100);
    
    param_indices = zeros(K,1);
    for k = 1:K
        % initialize E means to k-means centers
        [~,E_id] = min(abs(E_grid-E_kmeans(k)));
        
        % initialize flexibility to random value
        flx_id = randi(P);
        
        % vectorize parameter indexes
        param_indices(k) = sub2ind([M,P],E_id,flx_id);
    end

catch
    % initialize cluster assignment and parameters to random
    warning('K-means initialization failed: fallback to random.');
    Z = randi(K,1,N); 
    param_indices = randi(MP,K,1);
end

% compute sizes of clusters
num_samples = zeros(1,K);
for k = 1:K
    num_samples(k) = sum(Z==k);
end

% initialize results
K_trace = zeros(1,prm.niter);
LP_trace = zeros(1,prm.niter);
E_spl_iter = cell(prm.niter,1);
flx_spl_iter = cell(prm.niter,1);
w_spl_iter = cell(prm.niter,1);

% calculate nb. of printed characters to erase from display
[nb,str1,str2] = get_nb_char_to_print(prm);
for iter = 1:prm.niter
    
    % STEP 1: update cluster assignments (CRP)
    for i = 1:N
        current_k = Z(i);
        x_id = x_idx_all(i);
        
        % remove current observation
        num_samples(current_k) = num_samples(current_k) - 1;
        
        % delete cluster if empty (swap with last & pop)
        if num_samples(current_k)==0
            if current_k<K
                % swap cluster with greatest index with current
                Z(Z==K) = current_k;
                num_samples(current_k) = num_samples(K);
                param_indices(current_k) = param_indices(K);
            end
            % cluster with greatest index is deleted
            num_samples(K) = [];
            param_indices(K) = [];
            K = K - 1;
        end
        
        % calculate clusters log-probabilities
        log_prob_exist = log(num_samples(:)) + ...
            log_P_flat(x_id, param_indices)';
        
        % calculate log-probability for a new cluster
        log_prob_new = log(prm.alpha0) + log_marginal_per_bin(x_id);
        
        % normalize probabilities using stabilization with log-sum-exp
        log_probs = [log_prob_exist; log_prob_new];
        max_log = max(log_probs);
        probs = exp(log_probs - max_log);
        probs = probs / sum(probs);
        
        % sample cluster assignment
        new_k = randsample(length(probs), 1, true, probs);
        if new_k <= K
            % assign to existing cluster
            Z(i) = new_k;
            num_samples(new_k) = num_samples(new_k) + 1;
        else
            % create a new cluster
            K = K + 1;
            Z(i) = K;
            num_samples(K) = 1;
            
            % calculate posterior probabilities
            log_post_new = log_P_flat(x_id, :)' + log_G0_flat;
            p_new = exp(log_post_new - max(log_post_new));

            % sample new cluster's parameters from posterior probabilities
            param_indices(K) = randsample(MP, 1, true, p_new);
        end
    end
    
    % STEP 2: update cluster parameters
    for k = 1:K
        % Indices des données dans ce cluster
        idx_in_cluster = x_idx_all(Z == k); 
        
        if isempty(idx_in_cluster)
            continue; 
        end
        
        % Somme des log-vraisemblances pour toutes les données du cluster k
        % (1 x MP) = sum( (N_k x MP), 1 )
        sum_log_lik = sum(log_P_flat(idx_in_cluster, :), 1);
        
        % Ajout du prior
        log_posterior = sum_log_lik(:) + log(G0_flat);
        
        % Sampling
        post_prob = exp(log_posterior - max(log_posterior));
        param_indices(k) = randsample(MP, 1, true, post_prob);
    end
    
    % STEP 3: update alpha (optional)
    if isfield(hyp,'alpha_a')
         prm.alpha0 = update_alpha(prm.alpha0,N,K,hyp.alpha_a,hyp.alpha_b);
    end

    % calculate and store log-posterior
    LP_trace(iter) = calc_log_posterior(prm,log_P_flat,log_G0_flat,N,K,Z,...
        x_idx_all,num_samples,param_indices);
    
    % store nb. of clusters
    K_trace(iter) = K;
    
    % get cluster parameters for results storage and/or plot
    [E_id,flx_id] = ind2sub([M,P],param_indices);
    curr_E = E_grid(E_id);
    curr_flx = flx_grid(flx_id);
    curr_w = num_samples / N; 

    % sort parameters according to E
    [E_sort,ord] = sort(curr_E(:)');
    flx_sort = curr_flx(ord);
    w_sort = curr_w(ord);

    if iter>prm.burnin 
        % store parameters
        E_spl_iter{iter} = E_sort;
        flx_spl_iter{iter} = flx_sort;
        w_spl_iter{iter} = w_sort;
        
        % print progress
        if mod(iter,iv_iter)==0
            Kopt = round(mode(K_trace((prm.burnin+1):iter)));
            % fmt = repmat('%.2f ',1,K);
            % nb = eraseandfprintf(sprintf(['iter %i:, K=%i, Kopt=%i, ',...
            %    'alpha=%d\n\t E=',fmt,'\n\tflx=',fmt,'\n\t w=',fmt,'\n'],...
            %    iter,K,Kopt,prm.alpha0,E_sort,flx_sort,w_sort),nb);
            eraseandfprintf(sprintf(str1,runid,iter,prm.niter,K,Kopt,...
                prm.alpha0),nb);
        end

   elseif mod(iter,iv_iter)==0
       % print remaining iterations to burn
       % nb = eraseandfprintf(sprintf('iter %i: burnin phase %i/100\n',...
       %     iter,round(100*iter/prm.burnin)),nb);
       eraseandfprintf(sprintf(str2,runid,iter,round(100*iter/prm.burnin)),...
           nb);
   end

    % plot iteration
    if plotiter && mod(iter,iv_iter)==0
        PMF = calc_mixture_pmf(w_sort,{E_sort,flx_sort}, ...
            {E_grid,flx_grid},P_table);
        plt.objdpmm0 = plot_mixture(plt.objdpmm0,F_grid,PMF);
        plt.objdpmmk = ...
            plot_mixture_components(plt.objdpmmk,F_grid,E_sort,PMF);
        plt.objktr = plot_K_trace(plt.objktr,K_trace);
        plt.objpost = plot_logpost_trace(plt.objpost,LP_trace);
        if iter>prm.burnin
            plt.objkhist = ...
                plot_K_hist(plt.objkhist,K_trace((prm.burnin+1):iter));
            plt.objposthist = plot_logpost_hist(plt.objposthist,...
                LP_trace((prm.burnin+1):iter));
            % valid_idx = (1:prm.niter)>prm.burnin & (1:prm.niter)<=iter & ...
            %     K_trace==K;
            % E_spl_K = cell2mat(E_spl_iter(valid_idx));
            % flx_spl_K = cell2mat(flx_spl_iter(valid_idx));
            % w_spl_K = cell2mat(w_spl_iter(valid_idx));
            % 
            % plt.objEhist = plot_hist(plt.objEhist,6,2,8,E_spl_K,'E');
            % plt.objflxhist = plot_hist(plt.objflxhist,6,2,10,...
            %     flx2sig(flx_spl_K,const.k,const.S),'\\sigma r');
            %     plt.objwhist = plot_hist(plt.objwhist,6,2,12,w_spl_K,'w');
        end
        drawnow;
    end
end

% reshape output
K_values = unique(K_trace);
Kmax = max(K_values);
E_spl = cell(1,Kmax);
flx_spl = cell(1,Kmax);
w_spl = cell(1,Kmax);
for K = K_values
    valid_idx = (1:prm.niter)>prm.burnin & K_trace==Kopt;
    E_spl{K} = cell2mat(E_spl_iter(valid_idx));
    flx_spl{K} = cell2mat(flx_spl_iter(valid_idx));
    w_spl{K} = cell2mat(w_spl_iter(valid_idx));
end


function [alpha_E, beta_E] = fit_beta_prior(FRET)
% [alpha_E, beta_E] = fit_beta_prior(E_fret_data)
% Determine Beta hyperparameters (alpha, beta) using moments method from
% FRET data.
%
% FRET: observed FRET data

FRET = FRET(:);

% replace exterme values to prevent Inf/NaN.
FRET(FRET<=0) = eps;
FRET(FRET>=1) = 1 - eps;

% calculate moments
mu = mean(FRET);
sigma2 = var(FRET);

% calulate concentration parameter
nu = mu*(1-mu)/sigma2 - 1;

% check pathlogical case (too small variance)
if nu<=0
    warning(['Empirical variance is too small: weakly informative prior ',...
        'is used instead.']);
    nu = 1; 
end

alpha_E = mu*nu;
beta_E = (1-mu)*nu;


function log_post = calc_log_posterior(prm,log_P_flat,log_G0_flat,N,K,Z,...
    x_idx_all,num_samples,param_indices)
% log_post = calc_log_posterior(prm,log_P_flat,log_G0_flat,N,K,Z,...
%   x_idx_all,num_samples,param_indices)
%
% Calculate posterior log-probabilities of current DPMM.

% calculate log-likelihood
log_lik = 0;
for i = 1:N
    current_k = Z(i);
    x_id = x_idx_all(i);
    param_idx = param_indices(current_k);
    
    log_lik = log_lik + log_P_flat(x_id, param_idx);
end

% calculate log-prior for parameters
log_prior_param = 0;
for k = 1:K
    param_idx = param_indices(k);

    log_prior_param = log_prior_param + log_G0_flat(param_idx);
end

% calculate log-prior for cluster partition (given by CRP for DP)
% K*log(alpha) + sum_{k=1}^K log(Gamma(n_k)) - log(Gamma(N + alpha)) + 
%     log(Gamma(alpha))

alpha0 = prm.alpha0;
term_k_alpha = K * log(alpha0);
term_gamma_nk = sum(gammaln(num_samples)); % gammaln(n) = log(Gamma(n))
term_gamma_total = gammaln(alpha0) - gammaln(N + alpha0);

log_prior_partition = term_k_alpha + term_gamma_nk + term_gamma_total;

% calculate total log-posterior
log_post = log_lik + log_prior_param + log_prior_partition;


function plot_prior(plt,E_grid,flx_grid,G0_prior)

if ~iscell(plt)
    plt = {plt};
end
for i = 1:numel(plt)
    histogram2(plt{i}.axprior,'XBinEdges',center2edg(E_grid),'YBinEdges',...
        center2edg(flx_grid),'BinCounts',G0_prior,'DisplayStyle','tile');
    xlabel(plt{i}.axprior,'E mean');
    ylabel(plt{i}.axprior,'flexibility');
    title(plt{i}.axprior,'Prior distribution');
    
    % add color bar
    colormap(plt{i}.axprior,'turbo');
    colorbar(plt{i}.axprior);
    
    % adjust axes position to prevent axis label overlay
    plt{i}.axprior.OuterPosition = plt{i}.axprior.Position;
end

% refresh display
drawnow;



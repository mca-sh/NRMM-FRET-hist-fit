function [E_post,flx_post,w_post,E_spl,flx_spl,w_spl,plt] = ...
    EM_discrete_fit(FRET,F_counts,P_table,F_grid,E_grid,flx_grid,opt,plt)
% [E_post,flx_post,w_post,E_spl,flx_spl,w_spl,plt] = ...
%     EM_discrete_fit(FRET,F_counts,P_table,F_grid,E_grid,flx_grid,opt,plt)

% Get nb. of observations and grid dimensions
N = sum(F_counts);

% collect inference parameters and constants
prm = opt.prm;
const = opt.const;

% Initialize results
cell_Kmax = cell(1,prm.Kmax);
E = cell_Kmax;
flx = cell_Kmax;
w = cell_Kmax;
LL = -Inf(1,prm.Kmax);
BIC = Inf(1,prm.Kmax);
nb = 0;
for K = 1:prm.Kmax
    for i = 1:prm.ninit
        % parameter initialization
        [E_0,flx_0,w_0] = init_prm(K, F_grid, F_counts, E_grid, flx_grid);
        if any(w_0==0 )|| size(unique([E_0;flx_0]','rows'),1)<K
            continue
        end
        
        % EM inference of PMF mixture
        [E_i,flx_i,w_i,LL_i,plt] = em_discrete_mixture(F_counts,P_table,...
            E_grid,flx_grid,F_grid,K,E_0,flx_0,w_0,plt,prm,opt.plotiter);
        if any(w_i==0 )|| size(unique([E_i;flx_i]','rows'),1)<K
            continue
        end
        
        % calculate BIC for inferred parameters
        num_parameters = 3*K-1;
        BIC_i = calculate_bic(LL_i, num_parameters, N);
    
        % store ML estimate
        if LL_i>LL(K)
            [~,ord] = sort(E_i); % sort parameters according to E values
            E{K} = E_i(ord);
            flx{K} = flx_i(ord);
            w{K} = w_i(ord);
            LL(K) = LL_i;
            BIC(K) = BIC_i;
        end
        
        % print ML estimate
        if mod(i,10)==0
            nb0 = eraseandfprintf(sprintf(['Optimum for %i components ',...
                '(random init %i/%i):\n'],K,i,prm.ninit),nb);
            fmt = repmat('%.3f ',1,K);
            nb1 = fprintf(['\t  E=',fmt,'\n\tflx=',fmt,'\n\t  w=',fmt,'\n'],...
                E{K},flx{K},w{K});
            nb2 = fprintf('\tLog-likelihood: %.1E\n', LL(K));
            nb3 = fprintf('\tBIC: %.1E\n', BIC(K));
            nb = nb0+nb1+nb2+nb3;
        end
            
        % plot ML estimate
        if opt.plotfinal
            PMF = calc_mixture_pmf(w{K},{E{K},flx{K}},{E_grid,flx_grid},...
                P_table);
            plt.objdpmm0 = plot_mixture(plt.objdpmm0,F_grid,PMF);
            plt.objdpmmk = ...
                plot_mixture_components(plt.objdpmmk,F_grid,E{K},PMF);
            drawnow;
        end
    end
end

% Model selection (min. BIC)
[~,Kopt] = min(BIC);
eraseandfprintf(sprintf('Best out (minimum BIC): %d components\n',Kopt),...
    nb);

% Plot model selection
if ~isempty(plt)
    plot_BIC(plt.axprior,1:prm.Kmax,BIC,LL,Kopt);
end

% Determine uncertainty on parameters with bootstrapping
fprintf('Determine uncertainty with bootstrapping ');
E_spl = [];
flx_spl = [];
w_spl = [];
nb = 0;
F_edg = center2edg(F_grid);
num_parameters_opt = 3*Kopt-1;
for s = 1:prm.nSpl

    % resample observations
    FRET_s = randsample(FRET,N,true);
    F_counts_s = histcounts(FRET_s,F_edg);
    N_s = sum(F_counts_s);
   
    % parameter initialization (only once)
    if s==1
        E_0 = E{Kopt};
        flx_0 = flx{Kopt};
        w_0 = w{Kopt};
    end
    
    % EM inference
    [E_s,flx_s,w_s,LL_s,plt] = em_discrete_mixture(F_counts_s,P_table,...
        E_grid,flx_grid,F_grid,Kopt,E_0,flx_0,w_0,plt,prm,opt.plotiter);
    
    % calculate BIC for inferred model
    BIC_s = calculate_bic(LL_s, num_parameters_opt, N_s);
    
    % sort parameters according to E values
    [E_sort,ord] = sort(E_s);
    flx_sort = flx_s(ord);
    w_sort = w_s(ord);

    % store ML estimate for current sample
    E_spl = cat(1,E_spl,E_sort);
    flx_spl = cat(1,flx_spl,flx_sort);
    w_spl = cat(1,w_spl,w_sort);

    % print sample's ML estimate
    nb0 = eraseandfprintf(sprintf('(sample %i/%i):\n',s,prm.nSpl),nb);
    fmt = repmat('%.3f ',1,Kopt);
    nb1 = fprintf(['\t  E=',fmt,'\n\tflx=',fmt,'\n\t  w=',fmt,'\n'],E_sort,...
        flx_sort,w_sort);
    nb2 = fprintf('\tLog-likelihood: %.1E\n', LL_s);
    nb3 = fprintf('\tBIC: %.1E\n', BIC_s);
    nb = nb0+nb1+nb2+nb3;

    % plot sample's ML estimate
    if opt.plotiter
        PMF = calc_mixture_pmf(w_sort,{E_sort,flx_sort},{E_grid,flx_grid},...
            P_table);
        plt.objdpmm0 = plot_mixture(plt.objdpmm0,F_grid,PMF);
        plt.objdpmmk = ...
            plot_mixture_components(plt.objdpmmk,F_grid,E_sort,PMF);
        plt.objEhist = plot_hist(plt.objEhist,6,2,8,E_spl,'E');
        plt.objflxhist = plot_hist(plt.objflxhist,6,2,10,...
            flx2sig(flx_spl,const.k,const.S),'\\sigma r');
        plt.objwhist = plot_hist(plt.objwhist,6,2,12,w_spl,'w');

        drawnow;
    end
end

% Calculate bootstrap means and standard deviations
E_post = mean(E_spl,1);
flx_post = mean(flx_spl,1);
w_post = mean(w_spl,1);
dE_post = repmat(E_post,2,1) + ...
    repmat([-prm.nsig;prm.nsig],1,Kopt).*repmat(std(E_spl,[],1),[2,1]);
dflx_post = repmat(flx_post,2,1) + repmat([-prm.nsig;prm.nsig],1,Kopt)...
    .* repmat(std(flx_spl,[],1),[2,1]);
dw_post = repmat(w_post,2,1) + ...
    repmat([-prm.nsig;prm.nsig],1,Kopt).*repmat(std(w_spl,[],1),[2,1]);

% Plot model using bootstrap model parameters
if opt.plotfinal
    PMF = calc_mixture_pmf(w_post,{E_post,flx_post},{E_grid,flx_grid},...
        P_table);
    plt.objdpmm0 = plot_mixture(plt.objdpmm0,F_grid,PMF);
    plt.objdpmmk = plot_mixture_components(plt.objdpmmk,F_grid,E_post,PMF);
    plt.objEhist = plot_hist(plt.objEhist,6,2,8,E_spl,'E');
    plt.objflxhist = plot_hist(plt.objflxhist,6,2,10,...
        flx2sig(flx_spl,const.k,const.S),'\\sigma r');
    plt.objwhist = plot_hist(plt.objwhist,6,2,12,w_spl,'w');
end

% Concatenate output
E_post = [E_post;dE_post];
flx_post = [flx_post;dflx_post];
w_post = [w_post;dw_post];


function [E,flx,w,LL,plt] = em_discrete_mixture(F_counts,P_table,E_grid,flx_grid,...
    F_grid,K,E_0,flx_0,w_0,plt,prm,plotiter)
% EM_DISCRETE_MIXTURE EM algorithm optimized for a mixture of FRET PMFs

% default
iv_iter = 10; % nb. of iterations to ski^between each display

F_counts = F_counts(:); % Force column vector (L x 1)
[L, M, P] = size(P_table);

N = sum(F_counts);

% Add eps to prevent log(0) = -Inf
log_P_table = log(P_table + eps); 

% Flatten table to vectorize M-step (for speed)
log_P_table_flat = reshape(log_P_table, L, M*P);

% Initialization
E = E_0;
flx = flx_0;
w = w_0(:)';
LL = -Inf;

for iter = 1:prm.niter
    LL_prev = LL;

    % E-STEP : Calculate responsabilities
    PMF_k = zeros(L, K);
    
    for k = 1:K
        % Find grid indexes for current parameters
        E_id = get_nearest(E(k), E_grid);
        flx_id = get_nearest(flx(k), flx_grid);
        
        % Extract PMF for current parameters
        PMF_k(:,k) = w(k) * P_table(:, E_id, flx_id);
    end
    
    % add epsilon to prevent artificial explosion of LL for one component
    PMF_sum = eps + sum(PMF_k, 2); 
    
    % prevents division by 0
    valid_idx = PMF_sum > 0;
    gamma_k = zeros(L, K);
    gamma_k(valid_idx, :) = PMF_k(valid_idx, :) ./ PMF_sum(valid_idx);
    
    % calculate log-likelihood
    LL = sum(F_counts(valid_idx) .* log(PMF_sum(valid_idx)));
    
    % check for convergence
    if iter > 1
        delta_LL = (LL - LL_prev) / abs(LL_prev); 
        if abs(delta_LL) < prm.tol
            return
        end
    end

    % M-STEP: recalculate model parameters
    % update cluster weights
    N_k = sum(gamma_k .* F_counts, 1); 
    w = N_k / N; 
    
    % update parameters (E, flx)
    for k = 1:K
        if w(k) < 1e-6 
            % skip dead cluster
            continue; 
        end

        eff_weights = F_counts .* gamma_k(:, k);
       
        % Find (m, p) index that maxmizes Q
        log_Q = eff_weights' * log_P_table_flat;
        [max_log_Q, Q_id] = max(log_Q);
        
        if isfinite(max_log_Q)
            [E_id, flx_id] = ind2sub([M, P], Q_id);
            E(k) = E_grid(E_id);
            flx(k) = flx_grid(flx_id);
        end
    end
    
    % plot mixture
    if plotiter && mod(iter,iv_iter)==0
        PMF_plot = calc_mixture_pmf(w, {E,flx}, {E_grid,flx_grid}, P_table);
        plt.objdpmm0 = plot_mixture(plt.objdpmm0, F_grid, PMF_plot);
        plt.objdpmmk = plot_mixture_components(plt.objdpmmk, F_grid, E, PMF_plot);
        drawnow;
    end
end


function [E_0,flx_0,w_0] = init_prm(K,F_grid,F_counts,E_grid,flx_grid)

% draw K random indexes in E grid
F_0 = randsample(F_grid(F_counts>0),K,true);
E_0 = zeros(1,K);
for k = 1:K
    E_id_k = get_nearest(F_0(k),E_grid);
    E_0(k) = E_grid(E_id_k);
end

% 0.5 flexibility
flx_id = get_nearest(0.5,flx_grid);
flx_0 = ones(1,K)*flx_grid(flx_id);

w_0 = ones(1,K) / K; % uniform weights




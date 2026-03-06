function [E_post,sig_post,w_post,E_spl,sig_spl,w_spl,plt] = EM_GMM_fit(...
    FRET,F_counts,F_grid,opt,plt)

N = sum(F_counts);

% collect inference parameters and constants
prm = opt.prm;

cell_Kmax = cell(1,prm.Kmax);
E = cell_Kmax; 
sig = cell_Kmax;
w = cell_Kmax;
LL = -Inf(1,prm.Kmax);
BIC = Inf(1,prm.Kmax);
nb = 0;
for K = 1:prm.Kmax
    for i = 1:prm.ninit

        % parameter initialization
        [E_0, sig_0, w_0] = init_gmm_prm(K,F_counts,F_grid); 
        if any(w_0==0 )|| size(unique([E_0;sig_0]','rows'),1)<K
            continue
        end
        
        % EM inference of GMM
        [E_i,sig_i,w_i,LL_i,plt] = em_gaussian_mixture(F_grid,F_counts,K,...
            E_0,sig_0,w_0,plt,prm,opt.plotiter);
        if any(w_i==0 )|| size(unique([E_i;sig_i]','rows'),1)<K
            continue
        end
        
        % calculate BIC for ML estimate
        num_parameters = 3*K-1;
        BIC_i = calculate_bic(LL_i, num_parameters, N);
        
       % store ML estimate
        if LL_i>LL(K)
            [~,ord] = sort(E_i); % sort parameters according to E values
            E{K} = E_i(ord);
            sig{K} = sig_i(ord);
            w{K} = w_i(ord);
            LL(K) = LL_i;
            BIC(K) = BIC_i;
        end
        
        % print ML estimate
        if mod(i,10)==0
            nb0 = eraseandfprintf(sprintf(['Optimum for %i components ',...
                '(init %i/%i):\n'],K,i,prm.ninit),nb);
            fmt = repmat('%.3f ',1,K);
            nb1 = fprintf(['\t  E=',fmt,'\n\tsig=',fmt,'\n\t  w=',fmt,'\n'],...
                E{K},sig{K},w{K});
            nb2 = fprintf('\tLog-likelihood: %.1E\n', LL(K));
            nb3 = fprintf('\tBIC: %.1E\n', BIC(K));
            nb = nb0+nb1+nb2+nb3;
        end
            
        % plot ML estimate
        if opt.plotfinal
            PMF = calc_gmm_pmf(w{K},E{K},sig{K},F_grid);
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
sig_spl = [];
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
        sig_0 = sig{Kopt};
        w_0 = w{Kopt};
    end
    
    % EM inference of GMM
    [E_s,sig_s,w_s,LL_s,plt] = em_gaussian_mixture(F_grid,F_counts_s,Kopt,...
        E_0,sig_0,w_0,plt,prm,opt.plotiter);
    
    % calculate BIC for inferred model
    BIC_s = calculate_bic(LL_s, num_parameters_opt, N_s);
    
    % sort parameters according to E values
    [E_sort,ord] = sort(E_s);
    sig_sort = sig_s(ord);
    w_sort = w_s(ord);

    % store ML estimate for current sample
    E_spl = cat(1,E_spl,E_sort);
    sig_spl = cat(1,sig_spl,sig_sort);
    w_spl = cat(1,w_spl,w_sort);

    % print sample's ML estimate
    nb0 = eraseandfprintf(sprintf('(sample %i/%i):\n',s,prm.nSpl),nb);
    fmt = repmat('%.3f ',1,Kopt);
    nb1 = fprintf(['\t  E=',fmt,'\n\tsig=',fmt,'\n\t  w=',fmt,'\n'],E_sort,...
        sig_sort,w_sort);
    nb2 = fprintf('\tLog-likelihood: %.1E\n', LL_s);
    nb3 = fprintf('\tBIC: %.1E\n', BIC_s);
    nb = nb0+nb1+nb2+nb3;

    % plot sample's ML estimate
    if opt.plotiter
        PMF = calc_gmm_pmf(w_sort,E_sort,sig_sort,F_grid);
        plt.objdpmm0 = plot_mixture(plt.objdpmm0,F_grid,PMF);
        plt.objdpmmk = ...
            plot_mixture_components(plt.objdpmmk,F_grid,E_sort,PMF);
        plt.objEhist = plot_hist(plt.objEhist,6,2,8,E_spl,'E');
        plt.objflxhist = plot_hist(plt.objflxhist,6,2,10,sig_spl,...
            '\\sigma');
        plt.objwhist = plot_hist(plt.objwhist,6,2,12,w_spl,'w');

        drawnow;
    end
end

% Calculate bootstrap means and standard deviations
E_post = mean(E_spl,1);
sig_post = mean(sig_spl,1);
w_post = mean(w_spl,1);
dE_post = repmat(E_post,2,1) + ...
    repmat([-prm.nsig;prm.nsig],1,Kopt).*repmat(std(E_spl,[],1),[2,1]);
dsig_post = repmat(sig_post,2,1) + repmat([-prm.nsig;prm.nsig],1,Kopt)...
    .* repmat(std(sig_spl,[],1),[2,1]);
dw_post = repmat(w_post,2,1) + ...
    repmat([-prm.nsig;prm.nsig],1,Kopt).*repmat(std(w_spl,[],1),[2,1]);

% Plot model using bootstrap model parameters
if opt.plotfinal
    PMF = calc_gmm_pmf(w_post,E_post,sig_post,F_grid);
    plt.objdpmm0 = plot_mixture(plt.objdpmm0,F_grid,PMF);
    plt.objdpmmk = plot_mixture_components(plt.objdpmmk,F_grid,E_post,PMF);
    plt.objEhist = plot_hist(plt.objEhist,6,2,8,E_spl,'E');
    plt.objflxhist = plot_hist(plt.objflxhist,6,2,10,sig_spl,'\\sigma');
    plt.objwhist = plot_hist(plt.objwhist,6,2,12,w_spl,'w');
end

% Concatenate output
E_post = [E_post;dE_post];
sig_post = [sig_post;dsig_post];
w_post = [w_post;dw_post];


function [mu, sig, w, LL, plt] = em_gaussian_mixture(x_grid,x_counts,K,...
    mu_0,sig_0,w_0,plt,prm,plotiter)

% default
iv_iter = 10;     % nb. of iterations to skip between each display
min_sigma = 1E-5; % minimum peak stndard deviation

x_counts = x_counts(:); % Force column vector (L x 1)
N = sum(x_counts);
L = length(x_grid);

% Initialization
mu = mu_0;
sig = sig_0;
w = w_0(:)'; % Force row vector (1 x K)
LL = -Inf;

for iter = 1:prm.niter
    LL_prev = LL;
    
    % E-STEP: calculate responsibilities (gamma_k)
    % calculate PDFs of Gaussian components
    pmf_k = calc_gmm_pmf(w,mu,sig,x_grid)';
    pmf_sum = sum(pmf_k,2);
    
    % calculate responsibilities
    valid_idx = pmf_sum > 0;
    gamma_k = zeros(L, K); 
    gamma_k(valid_idx, :) = pmf_k(valid_idx, :) ./ pmf_sum(valid_idx);
    
    % calculate log-likelihood
    LL = sum(x_counts(valid_idx) .* log(pmf_sum(valid_idx)));
    
    % check for convergence
    if iter > 1
        delta_LL = (LL - LL_prev) / abs(LL_prev); 
        if abs(delta_LL) < prm.tol
            return
        end
    end
    
    % M-STEP : update Gaussian parameters
    % update weights of components
    N_k = sum(gamma_k.*x_counts, 1);
    w = N_k/N; 
    
    % update Gaussian means and standard deviations
    for k = 1:K
        if w(k) < 1e-6 
            % skip dead cluster
            continue
        end
        
        % update means
        mu(k) = sum(x_counts.*gamma_k(:,k).*x_grid')/N_k(k);
        
        % update standard deviations
        sig(k) = ...
            sqrt(sum(x_counts.*gamma_k(:, k).*(x_grid'-mu(k)).^2)/N_k(k));
        
        % prevent convergence to singularity (null standard deviation)
        sig(k) = max(sig(k),min_sigma);
    end
    
    % plot mixture
    if plotiter && mod(iter,iv_iter)==0
        PMF_plot = calc_gmm_pmf(w,mu,sig,x_grid);
        plt.objdpmm0 = plot_mixture(plt.objdpmm0,x_grid,PMF_plot);
        plt.objdpmmk = plot_mixture_components(plt.objdpmmk,x_grid,mu,...
            PMF_plot);
        drawnow;
    end
end


function [mu_0, sig_0, w_0] = init_gmm_prm(K,F_counts,F_grid)
FRET = F_grid(F_counts>0);
mu_0 = FRET(randi(length(FRET),1,K));

% standard deviation of about 10% of data range
sig_0 = ones(1,K)*(max(FRET)-min(FRET))/(4*K);

w_0 = ones(1,K)/K; % poids uniformes
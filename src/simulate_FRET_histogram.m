function [F_hist,E_k,flx_k,w_k] = simulate_FRET_histogram(N,K,const)

% CONSTANTS
BIN_EDGES = 0:0.02:1; % edges of FRET histogram bins
W_MIN = 0.05;         % minimum weigh of component in mixture model

% initialization
E_k = randsample(linspace(0,1,max([2*K,10])),K,'false');
flx_k = rand(1,K);
w_k = [];
while isempty(w_k) || any(w_k<W_MIN)
    w_k = rand(1,K);
    w_k = w_k/sum(w_k);
end
N_k = round(N*w_k);
FRET = [];

% compute bg means (use same N_sim to reduce variance)
bgA_tr = sum(noisedistrib(ones(const.npix,N)*const.bgA, const),1);
bgD_tr = sum(noisedistrib(ones(const.npix,N)*const.bgD, const),1);

% simulation
for k = 1:K
    % distribute E_m
    r_k = const.R0*((1/E_k(k))-1)^(1/6);
    sigr_k = flx2sig(flx_k(k),const.k,const.S);
    r = normrnd(r_k,sigr_k,[1,N_k(k)]);
    E = 1./(1+(r/const.R0).^6);

    % simulate noisy I_A and I_D with background and correction
    I_A = repmat(E*const.I0,const.npix,1);
    I_D = repmat((1-E)*const.I0/const.gamma,const.npix,1);
    
    n_ic_A = sum(noisedistrib((I_A+const.bgA),const),1)-bgA_tr(1:N_k(k));
    n_ic_D = sum(noisedistrib((I_D+const.bgD),const),1)-bgD_tr(1:N_k(k));

    FRET = cat(2,FRET,n_ic_A./(n_ic_A + const.gamma*n_ic_D));
end

% clean FRET data
if strcmp(const.noisetype,'P')
    FRET(FRET<0 | FRET>1) = [];
end

% build histogram
F_counts = histcounts(FRET,BIN_EDGES);
F_bins = mean([BIN_EDGES(1:end-1);BIN_EDGES(2:end)],1);
F_hist = [F_bins;F_counts];

[E_k,ord] = sort(E_k);
flx_k = flx_k(ord);
w_k = w_k(ord);

fprintf('Simulation done:\n');
fmt = repmat('%.2f ',1,K);
fprintf(['\t  E=',fmt,'\n\tflx=',fmt,'\n\t  w=',fmt,'\n'],E_k,flx_k,w_k);

function PMF = calc_gmm_pmf(w,mu,sig,bincenters)

% get bin edges
binedges = center2edg(bincenters);

% initialize PMF
K = length(mu);
PMF = zeros(K,length(bincenters));

for k = 1:K
    % calculate cum. PDF at each edge
    CDF_k = normcdf(binedges, mu(k), sig(k));
    
    % calculae PMF as difference between edges' CDFs
    PMF_k = diff(CDF_k);

    % normalize and weight PMF of component k
    PMF(k,:) = w(k)*PMF_k/sum(PMF_k);
end

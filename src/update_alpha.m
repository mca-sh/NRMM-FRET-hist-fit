function alpha_new = update_alpha(alpha, N, K, a, b)
    % Auxiliaire standard pour Gibbs DPMM (Escobar & West 1995)
    eta = betarnd(alpha + 1, N);
    pi_eta = (a + K - 1) / (a + K - 1 + N * (b - log(eta)));
    if rand < pi_eta
        alpha_new = gamrnd(a + K, 1 / (b - log(eta)));
    else
        alpha_new = gamrnd(a + K - 1, 1 / (b - log(eta)));
    end
end
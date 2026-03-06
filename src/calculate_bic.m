function BIC = calculate_bic(LL_max, nprm, nobs)

if nargin < 3
    error('La fonction calculate_bic nécessite la log-vraisemblance, le nombre de paramètres et le nombre d''observations.');
end

N = nobs;
p = nprm;

% Formule du BIC : BIC = -2 * LL_max + p * ln(N)
BIC = -2 * LL_max + p * log(N);

function x = dirichlet_sample(alpha)
    % sample from Dirichlet(alpha) where alpha is a 1 x K vector
    y = gamrnd(alpha, 1);
    x = y / sum(y);
end
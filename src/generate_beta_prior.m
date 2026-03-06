function G0 = generate_beta_prior(E_grid, alpha_beta, beta_beta)    
    % Convertir E_grid en colonne
    E = E_grid(:); 
    
    beta_density = betapdf(E, alpha_beta, beta_beta);
    
    % idx_zero = find(E == 0, 1);
    % if ~isempty(idx_zero)
    %     if alpha_beta < 1 && isinf(beta_density(idx_zero))
    %         % Utiliser la valeur du point voisin (approximation de la limite)
    %         beta_density(idx_zero) = betapdf(E(idx_zero + 1) / 2, alpha_beta, beta_beta) * 1e3; 
    %     elseif alpha_beta < 1 && beta_density(idx_zero) == 0 % S'assurer que ce n'est pas zéro
    %          beta_density(idx_zero) = max(beta_density) * 1e-6;
    %     end
    % end
    % 
    % % Trouver l'indice de E=1 (normalement le dernier)
    % idx_one = find(E == 1, 1);
    % if ~isempty(idx_one)
    %      % Si beta_beta < 1, la densité théorique tend vers l'infini.
    %     if beta_beta < 1 && isinf(beta_density(idx_one))
    %         % Utiliser la valeur du point voisin (approximation de la limite)
    %         beta_density(idx_one) = betapdf(1 - (1 - E(idx_one - 1)) / 2, alpha_beta, beta_beta) * 1e3;
    %     elseif beta_beta < 1 && beta_density(idx_one) == 0
    %         beta_density(idx_one) = max(beta_density) * 1e-6;
    %     end
    % end
    % 
    % 
    % % --- Discrétisation et Normalisation ---
    % % Pour obtenir une PMF discrète, il suffit de normaliser la densité (PDF) 
    % % pour que la somme totale soit 1.
    % 
    % % 1. S'assurer qu'aucune valeur n'est nulle (lissage implicite pour la MCMC)
    % beta_density = max(beta_density, 1e-10);
    % 
    % 2. Normaliser la somme pour obtenir la PMF (G0_pmf)
    G0 = beta_density' / sum(beta_density);

end



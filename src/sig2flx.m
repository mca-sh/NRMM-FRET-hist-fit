% function flx = sig2flx(sig,k,S)
% % flx = sig2flx(sig,k,S)
% %
% % Calculate flexibility (bounds [0;1]) of an RNA structure based on the 
% % standard deviation of the distribution of FRET distances. Flexibility is
% % modeled by sigmoïde: Flexibility = 1 / (1 + exp(-k * (sig - S)))
% %
% % sig: standard deviaition including RNA and fluorhpore linkers 
% %      structural flexibility.
% % k: slope factor (0.24 per Angstrom for a sharp slope)
% % S: tipping point in Angstrom, i.e., standard deviaition for 0.5 
% %    flexibility (6 Angstrom for typical linker flexibility)
% 
%     % Calcul du score de flexibilité (fonction sigmoïde)
% 
%     flx = 1 ./ (1 + exp(-k * (sig - S)));
% end

function flx = sig2flx(sig_r, k, S)
    % On évite log(0) en s'assurant que sig_r est strictement positif si nécessaire
    % Le modèle Log-Logistique force flx -> 0 quand sig_r -> 0
    
    % Formule modifiée : on travaille sur le log des distances
    flx = 1 ./ (1 + exp(-k * (log(sig_r) - log(S))));
    
    % Note : Cela équivaut mathématiquement à l'équation de Hill :
    % flx = (sig_r.^k) ./ (S^k + sig_r.^k);
end
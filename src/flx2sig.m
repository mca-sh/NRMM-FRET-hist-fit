% function sig = flx2sig(flx,k,S)
% % sig = flx2sig(flx,k,S)
% %
% % Calculate standard deviation of the distribution of FRET distances based 
% % on RNA structure flexibility (bounds [0;1]). Flexibility is
% % modeled by sigmoïde: Flexibility = 1 / (1 + exp(-k * (sig - S)))
% %
% % flx: RNA structure flexibility
% % k: slope factor (0.24 per Angstrom for a sharp slope)
% % S: tipping point in Angstrom, i.e., standard deviaition for 0.5 
% %    flexibility (6 Angstrom for typical linker flexibility)
% 
%     sig = S - log(1./flx - 1)/k;
% 
% end

function sig_r = flx2sig(flx, k, S)
    % Formule inverse adaptée au modèle Log-Logistique
    % Quand flx -> 0, le terme log(...) -> inf, donc l'exponentielle -> 0
    
    term = log(1./flx - 1);
    sig_r = exp(log(S) - term/k);
    
    % Version simplifiée équivalente :
    % sig_r = S * (flx ./ (1 - flx)).^(1/k);
end
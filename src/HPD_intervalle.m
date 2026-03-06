function [hpd_inter, hpd_longueur] = HPD_intervalle(echantillons, ...
    niveau_credibilite)
% HPD_INTERVALLE Calcule l'Intervalle de Plus Haute Densité Postérieure (HPD)
% à partir d'échantillons MCMC.
%
% ENTRÉES:
%   echantillons: Vecteur contenant les échantillons MCMC du posterior.
%   niveau_credibilite: Le niveau de crédibilité souhaité (ex: 0.95 pour 95%).
%
% SORTIES:
%   hpd_inter: Un vecteur [borne_inf, borne_sup] de l'intervalle HPD.
%   hpd_longueur: La longueur de l'intervalle HPD.

    % Vérification des entrées
    if niveau_credibilite <= 0 || niveau_credibilite >= 1
        error('Le niveau de crédibilité doit être entre 0 et 1.');
    end

    N = length(echantillons);
    
    % Calcul du nombre d'échantillons à inclure dans l'intervalle
    M = floor(N * niveau_credibilite); 
    
    if M < 2
        error('Pas assez d''échantillons pour calculer l''intervalle HPD.');
    end

    % 1. Trier les échantillons
    echantillons_tries = sort(echantillons);
    
    % Initialisation de la longueur minimale et de l'indice de départ
    longueur_min = inf;
    indice_debut_hpd = 1;
    
    % 2. Parcourir tous les sous-ensembles de taille M
    % Les sous-ensembles vont de l'indice i jusqu'à l'indice i + M - 1
    for i = 1 : (N - M + 1)
        % Le sous-ensemble est echantillons_tries(i) à echantillons_tries(i + M - 1)
        
        % Borne inférieure de l'intervalle actuel
        borne_inf = echantillons_tries(i);
        
        % Borne supérieure de l'intervalle actuel
        borne_sup = echantillons_tries(i + M - 1);
        
        % Calcul de la longueur de cet intervalle
        longueur_actuelle = borne_sup - borne_inf;
        
        % 3. Mettre à jour si cette longueur est la plus petite trouvée
        if longueur_actuelle < longueur_min
            longueur_min = longueur_actuelle;
            indice_debut_hpd = i;
        end
    end

    % 4. Construction de l'intervalle HPD final
    hpd_inter = [echantillons_tries(indice_debut_hpd), echantillons_tries(indice_debut_hpd + M - 1)];
    hpd_longueur = longueur_min;

end


%% Script MATLAB pour K=10 Composantes SANS Optimization Toolbox

% Nécessite : Statistics and Machine Learning Toolbox

% --- 0. Nettoyage et Initialisation ---
clear; close all;
rng(42); % Fixer le générateur de nombres aléatoires pour la reproductibilité

K = 2;

load(['/home/mimi/Documents/MATLAB/FRET-hist-fit-with-chatGPT/',...
    'MH45-FRET-data.mat'],'F_obs');

% Tronquer légèrement pour éviter 0 et 1 exacts pour le PDF
F_obs = max(0.001, min(0.999, F_obs')); 

% --- 2. Ajustement par Mélange Gaussien (GMM) ---
% Le GMM est utilisé pour approximer les moyennes et isoler les données.
disp('Ajustement du Mélange Gaussien (GMM)...');
GMModel = fitgmdist(F_obs, K, 'RegularizationValue', 0.001); % Régularisation pour la stabilité
P = posterior(GMModel, F_obs); % Probabilités d'appartenance

% Identifier les groupes par "hard assignment"
[~, cluster_indices] = max(P, [], 2);

% --- 3. Ajustement de la Distribution Beta sur les Données Isolées ---
fitted_params = struct('alpha', cell(1, K), 'beta', cell(1, K), 'p', cell(1, K));

disp('Ajustement des distributions Beta...');
for k = 1:K
    Data_GMM_k = F_obs(cluster_indices == k);
    
    if ~isempty(Data_GMM_k)
        try
            pd_beta_k = fitdist(Data_GMM_k, 'Beta');
            fitted_params(k).alpha = pd_beta_k.a;
            fitted_params(k).beta = pd_beta_k.b;
        catch
            warning(['Erreur Fit Beta pour l''état ', num2str(k)]);
            fitted_params(k).alpha = Alphas(k); % Utilisation de la valeur réelle en cas d'échec
            fitted_params(k).beta = Betas(k);
        end
        fitted_params(k).p = GMModel.ComponentProportion(k);
    else
        % Si un cluster GMM est vide (peu probable ici)
        fitted_params(k).alpha = 0; 
        fitted_params(k).beta = 0;
        fitted_params(k).p = 0; 
    end
end

% Tri des composantes ajustées par Efficacité moyenne (E = alpha / (alpha+beta))
E_fits = arrayfun(@(s) s.alpha / (s.alpha + s.beta), fitted_params);
[~, sorted_indices] = sort(E_fits);
fitted_params = fitted_params(sorted_indices);

% --- 4. Plot des Résultats ---
figure('Position', [100, 100, 800, 600]);
nbins = 80; 
histogram(F_obs, nbins, 'Normalization', 'pdf', 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
hold on;

E_range = linspace(0.001, 0.999, 500);
y_fit_total = zeros(size(E_range));
colors = jet(K); % Utiliser une palette de couleurs

disp('Paramètres ajustés (triés par Efficacité moyenne) :');
for k = 1:K
    alpha = fitted_params(k).alpha;
    beta = fitted_params(k).beta;
    p = fitted_params(k).p;
    
    % Calcul du mode Beta (E le plus probable)
    mode_E = (alpha - 1) / (alpha + beta - 2); 
    
    % Affichage des résultats
    fprintf('  État %2d: E mode = %.3f, (a=%.1f, b=%.1f), P = %.2f\n', k, mode_E, alpha, beta, p);
    
    % Plot de la composante Beta
    y_fit_k = p * betapdf(E_range, alpha, beta);
    plot(E_range, y_fit_k, 'LineWidth', 1.5, 'Color', colors(k, :));
    y_fit_total = y_fit_total + y_fit_k;
end

% Plot du Fit Total
plot(E_range, y_fit_total, 'LineWidth', 3, 'Color', 'k', 'LineStyle', '--');

title(['Ajustement Hybride de K=', num2str(K), ' Mélanges FRET (10 États Équidistants)']);
xlabel('Efficacité de FRET (E)');
ylabel('Densité de Probabilité');
legend_entries = [arrayfun(@(k) ['Comp. ', num2str(k)], 1:K, 'UniformOutput', false), {'Fit Total'}];
% Affichage réduit des légendes pour clarté
legend(legend_entries{[1, round(K/2), K, K+1]}, 'Location', 'NorthWest'); 
xlim([0 1]);
grid on;
hold off;
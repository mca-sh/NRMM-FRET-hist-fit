function mixture_pmf = calc_mixture_pmf(w, param, grids, P_table)
% mixture_pmf = calc_mixture_pmf(w, param, grids, P_table)
%
% Calculate PMF of individual K components in a mixture model. The model of 
% one component is characterized by D parameters. PMF values are read from 
% a table pre-calculated for series of values of each parameter.
%
% w: [1-by-K] relative weights of components in the mixture.
% param: {1-by-D}[1-by-K] parameters for each components.
% grids: {1-by-D+1} grids of x and parameter values for which PMFs were
%        calculated in the P_table.
% P_table: [N_x,N1,N2,...,ND] pre-calculated PMF(x|d1,d2,...,dD) for
%          ranges of x and parameter values.

    % get data dimensions
    K = length(param{1});
    D = numel(param);
    sz = size(P_table);

    % ensure normalization of weights
    w = w/sum(w);

    % initialize components' PMFs
    mixture_pmf = zeros(K,sz(1));

    % calculate the PMF for each component based on the weights and 
    % pre-calculated table
    for k = 1:K

        % build expression of indexes in table for a number D of dimensions
        str_var = '1:sz(1)';
        for d = 1:D
            str_var = cat(2,str_var,...
                [',',num2str(get_nearest(param{d}(k), grids{d}))]);
        end

        % evaluate expression
        id_k = eval(['sub2ind(sz,',str_var,')']);

        % collect PMF values for specified indexes in the table
        mixture_pmf(k,:) = w(k) * P_table(id_k);
    end

    % normalize components' PMFs by the total sum
    mixture_pmf = mixture_pmf ./ sum(mixture_pmf(:));
end



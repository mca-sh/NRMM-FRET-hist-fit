function plot_inference_success(ax,res,axis_grid)

if ~any(ndims(res)==[2,3])
        error('Dimensions of input argument 2 must be 2 or 3.');
end

% collect inferred nb. of components
sz = size(res,[1,2,3]);
K_inferred = zeros(sz);
for i_K = 1:sz(1)
    for i_d2 = 1:sz(2)
        for i_d3 = 1:sz(3)
            if ~isfield(res{i_K,i_d2,i_d3},'E')
                continue
            end
            K_inferred(i_K,i_d2,i_d3) = size(res{i_K,i_d2,i_d3}.E,2);
        end
    end
end

% plot inference success as a histogram
K_grid = axis_grid{1};
if sz(3)==1
    success = sum(K_inferred==repmat(K_grid',1,sz(2)),2);
    histogram(ax,'binedges',center2edg(K_grid),'bincounts',success/sz(2));
    xlabel(ax,'GT nb. of components');
    xticks(ax,K_grid);
    ylabel(ax,'success');
else
    alpha_grid = axis_grid{2};
    success = sum(K_inferred==repmat(K_grid',1,sz(2),sz(3)),[1,3]);
    histogram(ax,'binedges',center2edg(alpha_grid),'bincounts',...
        success/(sz(1)*sz(3)));
    xlabel(ax,'alpha');
    xscale(ax,'log');
    xticks(ax,alpha_grid);
    ylabel(ax,'success');
end

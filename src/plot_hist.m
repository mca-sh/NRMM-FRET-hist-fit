function ax = plot_hist(ax_prev,R,C,n,X,xlbl)
    nbins = 10;
    cnt = [];
    edg = [];
    K = size(X,2);
    clr = get_component_color(K);
    for k = 1:K
        [cnt_k,edg_k] = histcounts(X(:,k)',nbins);
        cnt = cat(1,cnt,cnt_k);
        edg = cat(1,edg,edg_k);
    end

    for k_prev = 1:numel(ax_prev)
        if k_prev==1
            fig = ax_prev(k_prev).Parent;
        end
        delete(ax_prev(k_prev));
    end
    
    ax = [];
    r = ceil(n/C);
    c = n-(r-1)*C;
    for k = 1:K
        ax_k = create_subaxes(fig,[R,C*K,(r-1)*C*K+(c-1)*K+k]);
        ax = cat(2,ax,ax_k);

        histogram(ax_k,'bincounts',cnt(k,:),'binedges',edg(k,:),...
            'facecolor',clr(k,:));
        X_post_k = median(X(:,k));
        plot(ax_k,[X_post_k,X_post_k],ax_k.YLim,'-k','linewidth',2);
        
        xlbl_k = sprintf([xlbl,'_%i'],k);
        xlabel(ax_k,xlbl_k);
        if k==1
            ylabel(ax_k,'count');
        end
        nonemptybins = find(cnt(k,:)>0);
        n_nonempty = numel(nonemptybins);
        nonemptyedg = unique(reshape(...
            [edg(k,nonemptybins);edg(k,nonemptybins+1)],1,2*n_nonempty));
        x1 = nonemptyedg(1);
        x2 = nonemptyedg(end);
        
        binsz = edg(k,2)-edg(k,1);
        xlim(ax_k,[x1-binsz,x2+binsz]);

        % adjust axes position to prevent overlay of axis labels
        ax_k.OuterPosition = get_subplot_position(ax_k.UserData);
    end
end

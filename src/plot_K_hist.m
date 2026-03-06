function obj = plot_K_hist(obj0,K_trace)
    edg = 0.5:max(K_trace(K_trace>0))+0.5;
    cnt = histcounts(K_trace(K_trace>0),edg);
    if isa(obj0,'matlab.graphics.axis.Axes')
        ax = obj0;
        obj = histogram(ax,'bincounts',cnt,'binedges',edg);
        xlabel(ax,'nb. of components');
        ylabel(ax,'count');
        xlim(ax,[0,max(K_trace)+1]);

        % adjust axes position to prevent overlay of axis labels
        ax.OuterPosition = get_subplot_position(ax.UserData);
    else
        obj = obj0;
        ax = obj.Parent;
        set(obj,'bincounts',cnt,'binedges',edg);
        xlim(ax,[0,max(K_trace)+1]);
    end
end
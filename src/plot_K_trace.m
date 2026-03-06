function obj = plot_K_trace(obj0,K_trace)
    if isa(obj0,'matlab.graphics.axis.Axes')
        ax = obj0;
        yyaxis(ax,'left');
        obj = plot(ax,1:find(K_trace>0,1,'last'),K_trace(K_trace>0));
        xlabel(ax,'iteration');
        ylabel(ax,'nb. of components');
        xlim(ax,[0,numel(K_trace)+1]);
        ylim(ax,[0,max(K_trace)+1]);

        % adjust axes position to prevent overlay of axis labels
        ax.OuterPosition = get_subplot_position(ax.UserData);
    else
        obj = obj0;
        ax = obj.Parent;
        yyaxis(ax,'left');
        set(obj,'XData',1:find(K_trace>0,1,'last'),'YData',...
            K_trace(K_trace>0))
        xlim(ax,[0,numel(K_trace)+1]);
        ylim(ax,[0,max(K_trace)+1]);
    end
end
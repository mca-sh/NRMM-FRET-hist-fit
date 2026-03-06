function obj = plot_logpost_trace(obj0,logpost_trace)
    posprob = logpost_trace~=0;
    if isa(obj0,'matlab.graphics.axis.Axes')
        ax = obj0;
        yyaxis(ax,'right');

        obj = plot(ax,1:find(posprob,1,'last'),logpost_trace(posprob));
        ylabel(ax,'log-posterior');
        xlim(ax,[0,numel(logpost_trace)+1]);

        % adjust axes position to prevent overlay of axis labels
        ax.OuterPosition = get_subplot_position(ax.UserData);
    else
        obj = obj0;
        ax = obj.Parent;
        yyaxis(ax,'right');
        set(obj,'XData',1:find(posprob,1,'last'),'YData',...
            logpost_trace(posprob));
    end
end
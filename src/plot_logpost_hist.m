function obj = plot_logpost_hist(obj0,LP_trace)
    posprob = LP_trace~=0;
    [cnt,edg] = histcounts(LP_trace(posprob));
    LP_median = median(LP_trace(posprob));
    if all(isa(obj0,'matlab.graphics.axis.Axes'))
        ax = obj0;
        obj1 = histogram(ax,'bincounts',cnt,'binedges',edg);
        obj2 = line(ax,[LP_median,LP_median],ax.YLim,'linestyle','--',....
            'linewidth',1.5);
        xlabel(ax,'log-posterior');
        ylabel(ax,'count');
        
        % print mode of log-posterior
        obj3 = text(ax,mean(ax.XLim),mean(ax.YLim),0,['LP_{mode}=',...
            num2str(LP_median)],'horizontalalignment','center');
        obj = [obj1,obj2,obj3];

        % adjust axes position to prevent overlay of axis labels
        ax.OuterPosition = get_subplot_position(ax.UserData);
    else
        obj = obj0;
        ax = obj(1).Parent;
        set(obj(1),'bincounts',cnt,'binedges',edg);
        set(obj(2),'xdata',[LP_median,LP_median],'ydata',ax.YLim);
        obj(3).String = num2str(LP_median);
    end
end
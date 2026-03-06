function obj = plot_mixture(obj0, x_bins, mixture_pmf)
    if isa(obj0,'matlab.graphics.axis.Axes')
        ax = obj0;
        obj = plot(ax,x_bins, sum(mixture_pmf,1),'-r','linewidth',1.5);

        xlim(ax,[0,1]);
        xlabel(ax,'apparent FRET');
        ylabel(ax,'Normalized counts');
        title(ax,'Mixture');

        % adjust axes position to prevent overlay of axis labels
        ax.OuterPosition = get_subplot_position(ax.UserData);
    else
        obj = obj0;

        % update mixture plot
        obj.YData = sum(mixture_pmf,1);
    end
end



function obj = plot_mixture_components(obj0, x_bins, E, mixture_pmf)
    if isa(obj0,'matlab.graphics.axis.Axes')
        ax = obj0;
        obj = plot_components(ax,E,x_bins,mixture_pmf);

        xlim(ax,[0,1]);
        xlabel(ax,'apparent FRET');
        ylabel(ax,'Normalized counts');
        title(ax,'Components');

        % adjust axes position to prevent overlay of axis labels
        ax.OuterPosition = get_subplot_position(ax.UserData);

    else
        obj = obj0;
        ax = obj.Parent;

        % delete previous mixture's components
        delete(ax.Children(1:end-1)); % skip experimental data (last)

        % plot new mixture's components
        hline = plot_components(ax,E,x_bins,mixture_pmf);
        obj = hline(1);
    end

end


function h = plot_components(ax,mus,x_bins,mixture_pmf)
    K = length(mus);
    clr = get_component_color(K);
    h = [];
    for k = 1:length(mus)
        h = cat(2,h,plot(ax,x_bins, mixture_pmf(k,:),'color',clr(k,:),...
            'LineWidth',1.5));
    end
    for k = 1:length(mus)
        h = cat(2,h,plot(ax,[mus(k),mus(k)], ax.YLim,'color',clr(k,:),...
            'linestyle','--','LineWidth',1.5));
    end
end

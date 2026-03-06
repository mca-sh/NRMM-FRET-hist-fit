function plotPMFlookup(ax,edgx,edgy,P_table)
    % heatmap
    histogram2(ax(1),'xbinedges',edgx,'ybinedges',edgy,'bincounts',...
        P_table','displaystyle','tile');
    xlim(ax(1),edgx([1,end]));
    ylim(ax(1),edgy([1,end]));
    xlabel(ax(1),'apparent FRET')
    ylabel(ax(1),'E mean')
    title(ax(1),'Pre-calculated PMFs (0.5 flexibility)');

    % adjust axes position to prevent overlay of axis labels
    ax(1).OuterPosition = get_subplot_position(ax(1).UserData);

    % add colorbar
    cb = colorbar(ax(1));
    colormap(ax(1),'turbo');

    % ax(1).OuterPosition = tightPosition(ax(1),'includelabels',true,'units',...
    %     ax(1).Units);

    % 2D plot
    szY = size(P_table,2);
    xmeans = mean([edgx(1:end-1);edgx(2:end)],1);
    id = round(linspace(1,szY,9));
    n = 0;
    for y = id
        n = n+1;
        clr = 0.4+zeros(1,3)+0.3*mod(n,2);
        plot(ax(2),xmeans',P_table(:,y),'color',clr,'linewidth',1.5);
    end
    xlabel(ax(2),'apparent FRET')
    ylabel(ax(2),'Probability')
    ylim(ax(2),[0,max(P_table(:))/2]);

    % adjust axes position to prevent overlay of axis labels
    ax(2).OuterPosition = get_subplot_position(ax(2).UserData);
end
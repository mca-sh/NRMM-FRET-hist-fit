function plot_BIC(ax,Ks,BIC,LL,Kopt)

% plot BIC against nb. of components
yyaxis(ax,'left');
scatter(ax,Ks,BIC,'marker','+','sizedata',50,'linewidth',1);
scatter(ax,Kopt,BIC(Kopt),'marker','o','markeredgecolor','red','sizedata',...
    100,'linewidth',1.5);
ylabel(ax,'BIC');

% plot log-likelihood against nb. of components
yyaxis(ax,'right');
scatter(ax,Ks,LL,'marker','+','sizedata',50,'linewidth',1,...
    'markeredgecolor',[1,0.5,0]);
ylabel(ax,'log-likelihood');

xlim(ax,[0.5,max(Ks)+0.5]);
xticks(ax,Ks);
xlabel(ax,'nb. of components');

% print minimum BIC value
text(ax,mean(ax.XLim),mean(ax.YLim),0,...
    ['BIC_{opt}=',num2str(BIC(Kopt))],'horizontalalignment','center');

% adjust axes position to prevent label overlay
ax.OuterPosition = ax.Position;
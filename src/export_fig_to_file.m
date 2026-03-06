function export_fig_to_file(fig,fle,meth)

[dest,fnm,~] = fileparts(fle);
prevtheme = fig.Theme;
fig.Theme = 'light';
savefig(fig,addnbtofilename([dest,filesep,fnm,'_',meth,'.fig']));
print(fig,addnbtofilename([dest,filesep,fnm,'_',meth,'.png']),'-dpng');
fig.Theme = prevtheme;
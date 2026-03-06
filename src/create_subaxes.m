function ax = create_subaxes(fig,pos_on_grid)

ax = axes(fig,'position',get_subplot_position(pos_on_grid),'nextplot',...
    'add','userdata',pos_on_grid);
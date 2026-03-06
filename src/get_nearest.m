function E_id = get_nearest(E,E_grid)
    [~,E_id] = min(abs(E_grid-E));
end
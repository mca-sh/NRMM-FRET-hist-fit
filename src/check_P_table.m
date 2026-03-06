function [P_table,E_grid,F_hist,fprm] = check_P_table(P_table,F_hist,Nsim,...
    fprm,flx_grid,const,fle)

% clean FRET histogram counts
F_hist(:,F_hist(1,:,1)<0 | F_hist(1,:,1)>1) = [];
F_bins = F_hist(1,:,1);

% check validity of pre-calculated PMF grid
E_grid = F_bins;
if ~isequal({fprm.E_grid0,fprm.flx_grid0,fprm.xbins0,fprm.const0},...
        {E_grid,flx_grid,F_bins,const})
    % pre-calculate PMF grid
    P_table = calc_E_PMF(E_grid,flx_grid,F_bins,Nsim,const);
    
    % save table to file
    flx_grid0 = flx_grid;
    E_grid0 = E_grid;
    xbins0 = F_bins;
    const0 = const;
    save(fle,'P_table','flx_grid0','E_grid0','xbins0','const0','-append');

    fprm.E_grid0 = E_grid0;
    fprm.flx_grid0 = flx_grid0;
    fprm.xbins0 = xbins0;
    fprm.const0 = const0;
end
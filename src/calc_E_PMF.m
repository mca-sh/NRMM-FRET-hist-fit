function P_table = calc_E_PMF(E_grid, flx_grid, x_bins, Nsim, const)
% Simulate PMF for each FRET means and flexibility values in grids

% default
plotit = false; % plot simulation (if true, replace "parfor" by "for" loop)

% turn off warning about temporary variable
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');

% get grid dimensions
M = length(E_grid);
N = length(flx_grid);
L = length(x_bins);
P_table = zeros(L, M, N);

% compute FRET histogram bin edges
edg = center2edg(x_bins);

% initialize plot (only without parfor)
if plotit
    fig = figure('WindowStyle','docked');
    ax_r = create_subaxes(fig,[3,3,1]);
    ax_E = create_subaxes(fig,[3,3,4]);
    ax_fret = create_subaxes(fig,[3,3,7]);
    ax_pcA = create_subaxes(fig,[3,3,2]);
    ax_icA = create_subaxes(fig,[3,3,5]);
    ax_pcD = create_subaxes(fig,[3,3,3]);
    ax_icD = create_subaxes(fig,[3,3,6]);

    n_ic_max = 1.5*const.npix*(pc2ic(const.I0,const)-const.b);
    n_pc_max = const.npix*const.I0;

    clr = orderedcolors("gem");
    n_ref = get_nearest(0.5,flx_grid);
    
else % dummy plot variables need to be defined when using parfor, even if 
     % not used
    ax_r = [];
    ax_E = [];
    ax_fret = [];
    ax_pcA = [];
    ax_icA = [];
    ax_pcD = [];
    ax_icD = [];

    n_ic_max = [];
    n_pc_max = [];
    
    clr = [];
    n_ref = [];
end

% compute bg means (use same N_sim to reduce variance)
bgA_tr = sum(noisedistrib(ones(const.npix,Nsim)*const.bgA, const),1);
bgD_tr = sum(noisedistrib(ones(const.npix,Nsim)*const.bgD, const),1);

if plotit
    c = 0; % initialize plot color running index
end
parfor m = 1:M % to be replaced by for m = 1:M if plotit==true
    E_m = E_grid(m);

    for n = 1:N
        flx_m = flx_grid(n);

        % distribute E_m
        r_m = const.R0*(((1/E_m)-1)^(1/6));
        sigr_m = flx2sig(flx_m,const.k,const.S);
        r = normrnd(r_m,sigr_m,[1,Nsim]);
        E = 1./(1+(r/const.R0).^6);

        % simulate noisy I_A and I_D with background and correction
        n_pc_A = repmat(E.*const.I0,const.npix,1);
        n_pc_D = repmat((1-E).*const.I0/const.gamma,const.npix,1);
        
        n_ic_A = sum(noisedistrib((n_pc_A+const.bgA),const),1)-bgA_tr;
        n_ic_D = sum(noisedistrib((n_pc_D+const.bgD),const),1)-bgD_tr;
        
        % calculate FRET ratio
        FRET_sim = n_ic_A./(n_ic_A + const.gamma*n_ic_D);
        if strcmp(const.noisetype,'P')
            FRET_sim(FRET_sim<0 | FRET_sim>1) = [];
        end
        
        % compute normalized FRET histogram (PMF)
        cnts = histcounts(FRET_sim,edg);
        if sum(cnts)>0
            P_table(:,m,n) = cnts/sum(cnts);
        end

        % plot simulation
        if plotit && n==n_ref
            c = c+1;
            if c>size(clr,1)
                c = 1;
            end

            % histogram FRET distances
            histogram(ax_r,r,linspace(0,2*const.R0,100),'facecolor',...
                clr(c,:));
            plot(ax_r,[r_m,r_m],ax_r.YLim,'linestyle','--','color',...
                clr(c,:));
            xlabel(ax_r, 'r (A)');

            % histogram mean FRET efficiencies
            histogram(ax_E,E,linspace(0,1,100),'facecolor',clr(c,:));
            plot(ax_E,[E_m,E_m],ax_E.YLim,'linestyle','--','color',...
                clr(c,:));
            xlabel(ax_E, 'E');

            % histogram apparent FRET
            histogram(ax_fret,'binedges',edg,'bincounts',cnts,'facecolor',...
                clr(c,:));
            plot(ax_fret,[E_m,E_m],ax_fret.YLim,'linestyle','--','color',...
                clr(c,:));
            xlabel(ax_fret, 'app. FRET');

            % histogram photon counts
            histogram(ax_pcA,sum(n_pc_A,1),linspace(0,n_pc_max,100),...
                'facecolor',clr(c,:));
            histogram(ax_pcD,sum(n_pc_D,1),linspace(0,n_pc_max,100),...
                'facecolor',clr(c,:));
            xlabel(ax_pcA, 'acc PC');
            xlabel(ax_pcD, 'don PC');

            % histogram image counts
            histogram(ax_icA,n_ic_A,linspace(0,n_ic_max,100),'facecolor',...
                clr(c,:));
            histogram(ax_icD,n_ic_D,linspace(0,n_ic_max,100),'facecolor',...
                clr(c,:));
            xlabel(ax_icA, 'acc IC');
            xlabel(ax_icD, 'don IC');

            ylabel([ax_r,ax_E,ax_fret,ax_pcA,ax_pcD,ax_icA,ax_icD],...
                'counts');
        end
    end
    fprintf('  grid index %d / %d done\n',m*N,N*M);
end

% adjust axes position to rpevent label overlay
if plotit
    ax_r.OuterPosition = ax_r.Position;
    ax_E.OuterPosition = ax_E.Position;
    ax_fret.OuterPosition = ax_fret.Position;
    ax_pcA.OuterPosition = ax_pcA.Position;
    ax_icA.OuterPosition = ax_icA.Position;
    ax_pcD.OuterPosition = ax_pcD.Position;
    ax_icD.OuterPosition = ax_icD.Position;
end

fprintf('Probability grid was successfully computed!\n');
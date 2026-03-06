function [res,plt] = FRET_hist_discrete_fit(F_hist,P_table,E_grid,flx_grid,...
    opt,plt)
% [res,plt] = FRET_hist_discrete_fit(F_hist,P_table,E_grid,flx_grid,opt,plt)
%
% Cluster histogrammed data with an optimum nb of clusters. Data is modeled
% with a mixture of pre-calculated simulated discrete distributions.
% Simulation include structure flexibility and realistic EMCCD or 
% single-photon APD 
% noise.

% initialize output
res = [];

% get histogram data
F_grid = F_hist(1,:);
F_counts = rescale_hist_counts(F_hist(2,:),opt.Nmax);

% initialize plot
if opt.plotiter || opt.plotfinal
    if strcmp(opt.inference_method,'DPMM')
        plt = cell(1,opt.prm.nruns);
        for i = 1:opt.prm.nruns
            plt{i} = init_plot(opt,[F_grid;F_counts],flx_grid,E_grid,...
                P_table);
            % plt{i}.fig.WindowStyle = 'normal'; % parallel computing ignore 
            %                                    % docked windows and create 
            %                                    % new ones
        end
    else
        plt = init_plot(opt,[F_grid;F_counts],flx_grid,E_grid,P_table);
    end
end

% unfold histogram for inference
FRET = [];
for i = 1:length(F_grid)
    FRET = cat(1,FRET,repmat(F_grid(i), F_counts(i), 1)); 
end

% perform inference
switch opt.inference_method

    case 'EM'
        [E,flx,w,E_spl,flx_spl,w_spl,plt] = EM_discrete_fit(FRET,F_counts,...
            P_table,F_grid,E_grid,flx_grid,opt,plt);

    case 'DPMM'
        [E,flx,w,E_spl,flx_spl,w_spl,plt,bestrun] = DPPM_discrete_fit(FRET,...
            P_table,F_grid,E_grid,flx_grid,opt,plt);
        res.bestrun = bestrun;

    case 'EM-GMM'
        [E,sig,w,E_spl,sig_spl,w_spl,plt] = EM_GMM_fit(FRET,F_counts,...
            F_grid,opt,plt);

    otherwise
        error('Inference method must be ''EM-GMM'', ''EM'', or ''DPMM''.');
end

% collect results to return
res.E = E;
res.w = w;
res.E_spl = E_spl;
res.w_spl = w_spl;
switch opt.inference_method
    case {'EM','DPMM'}
        prm2_str = 'flx';
        prm2 = flx;
        res.flx = flx;
        res.flx_spl = flx_spl;
    case 'EM-GMM'
        prm2_str = 'sig';
        prm2 = sig;
        res.sig = sig;
        res.sig_spl = sig_spl;
end

% print best out
Kopt = size(E,2);
fprintf('Optimum nb. of components: %d\n', Kopt);
fmt = repmat('%.3f ',1,Kopt);
fprintf(['\t  E=',fmt,'\n\t',prm2_str,'=',fmt,'\n\t  w=',fmt,'\n'],E(1,:),...
    prm2(1,:),w(1,:));
fmt = repmat('(%.3f;%.3f) ',1,Kopt);
switch opt.inference_method
    case {'EM','EM-GMM'}
        fprintf('Confidence intervals (1 bootstrap deviaition):\n');
    case 'DPMM'
        fprintf('Credibility intervals:\n');
end
fprintf(['\t  E=',fmt,'\n\t',prm2_str,'=',fmt,'\n\t  w=',fmt,'\n'],E([2,3],:),...
    prm2([2,3],:),w([2,3],:));

fprintf('Process successfully terminated!\n');


function F_counts = rescale_hist_counts(F_counts,Nmax)
    N = sum(F_counts);
    if N>Nmax
        N = Nmax;
    end
    F_counts = round(N*F_counts/sum(F_counts));


function plt = init_plot(opt,F_hist,flx_grid,E_grid,P_table)

% build figure and axes
plt.fig = figure('windowstyle','docked','theme', 'light');
plt.axpmfmap = create_subaxes(plt.fig,[4,3,1]);
plt.axpmfhist = create_subaxes(plt.fig,[4,3,4]);
plt.axflxmdl = create_subaxes(plt.fig,[4,3,2]);
plt.axprior = create_subaxes(plt.fig,[4,3,5]);
plt.axdpmm0 = create_subaxes(plt.fig,[4,3,3]);
plt.axdpmmk = create_subaxes(plt.fig,[4,3,6]);
plt.axtr = create_subaxes(plt.fig,[4,2,5]);
plt.axkhist = create_subaxes(plt.fig,[4,4,13]);
plt.axposthist = create_subaxes(plt.fig,[4,4,14]);
plt.axEhist = create_subaxes(plt.fig,[6,2,8]);
plt.axflxhist = create_subaxes(plt.fig,[6,2,10]);
plt.axwhist = create_subaxes(plt.fig,[6,2,12]);

% initialize plot objects structure
plt.objdpmm0 = plt.axdpmm0;
plt.objdpmmk = plt.axdpmmk;
plt.objktr = plt.axtr;
plt.objpost = plt.axtr;
plt.objkhist = plt.axkhist;
plt.objposthist = plt.axposthist;
plt.objEhist = plt.axEhist;
plt.objflxhist = plt.axflxhist;
plt.objwhist = plt.axwhist;

% hide unused axes
if contains(opt.inference_method,'EM')
    set([plt.axkhist,plt.axtr,plt.axposthist],'visible','off');
end

% plot FRET histogram
F_counts = F_hist(2,:);
F_grid = F_hist(1,:);
F_normcounts = F_counts/sum(F_counts);
hb = bar(plt.axdpmm0,F_grid,F_normcounts,'hist');
hb.FaceColor = [0.8,0.8,0.8];
hb = bar(plt.axdpmmk,F_grid,F_normcounts,'hist');
hb.FaceColor = [0.8,0.8,0.8];
xlim([plt.axdpmm0,plt.axdpmmk],[-0.2,1.2]);


if strcmp(opt.inference_method,'EM-GMM')
    set([plt.axflxmdl,plt.axpmfmap,plt.axpmfhist],'visible','off');
else
    % plot flexibilty model
    plot(plt.axflxmdl,flx2sig(flx_grid,opt.const.k,opt.const.S),flx_grid,...
        'marker','o');
    xlabel(plt.axflxmdl,'\sigma_r (A)');
    ylabel(plt.axflxmdl,'Flexibility');
    
    % adjust axes position to prevent overlay of axis labels
    plt.axflxmdl.OuterPosition = plt.axflxmdl.Position;
    
    % Plot PMF grid
    id_flx05 = get_nearest(0.5,flx_grid);
    plotPMFlookup([plt.axpmfmap,plt.axpmfhist],center2edg(E_grid),...
        center2edg(F_grid),P_table(:,:,id_flx05));
end

% refresh display
drawnow;

% --- README --------------------------------------------------------------
% SCRIPT_FIT_HIST is desinged to fit FRET histograms with realistic-noise
% mixture models.
%
% Histogram data must first be formatted to a suitable .mat file. Run
% script "import_from_txt.m" to convert a text file to a suitable .mat file
% prior runing this script. Column 1 of the text file must contains the bin 
% centers and the column 2 the bin counts.
%
% Adjust user parameters in the section "USER PARAMETERS" below prior
% running this script.
%
% Because the noise model possesses no analytical form, model's 
% probabilitities as pre-calculated using data simulation for fixed grids 
% of parameters (peak means and structure flexibilities). This lookup table 
% is calculated once for the input noise model parameters (progress shown 
% in prompt with "grid index xxxx / xxxx done") and saved to the analysis 
% mat file. If the noise parameters change, the table is recalculated.
%
% For speed, table and DPMM calculations uses the parallel computing tool.
% It is highly recommended to install and use it.
%
% When chosen, several runs of DPMM inference are executed in parallel and 
% the best outcome is selected. A figure is created for each run and
% exported to png files once the inference is completed.
%
% Inference results are saved to a *_results_[method].mat file and a 
% graphical summary is exported to .png and .fig image files, all in a 
% subfolder named after the source file.
%
% The graphical summary shows from left to right and top to bottom:
% (1) a slice (at flexibility=0.5) of the pre-calculated PMF of an 
%     histogram peak in function of its FRET mean, 
% (2) the flexibility in function of the standard FRET distance
%     deviation of the structure, 
% (3) the experimental FRET histogram with the fitted mixture, 
% (4) a superposition of 2D representations of the pre-calculated PMF of an 
%     histogram peak for some values of the FRET mean and at 0.5
%     flexibility,
% (5) the prior distribution of the model parameters (used in DPMM)
% (6) the experimental FRET histogram with the individual components of 
%     fitted mixture, 
% (7) iteration traces of the inferred number of peaks (K) and the 
%     posterior log-probability of the inferred model
% (8) distribution of the means of the inferred FRET peaks accross the 
%     iteration where K=Kopt
% (9) distribution of K accross the iterations
% (10) distribution of posterior log-prob accross the iteration
% (11) distribution of the standard distance deviation of the inferred
%      peaks
% (12) distribution of the relative weight of the peaks in the inferred
%      mixture.
% -------------------------------------------------------------------------

% --- USER PARAMETERS (to be set prior running the script) ----------------
inference_method = 'DPMM'; % 'DPMM' (infinite noise-realist mixture models)
                           % 'EM' (EM+BIC of noise-realist mxture models)
                           % 'EM-GMM' (EM+BIC of Gaussian mixture models)
recordingtype = 'SPAD'; % 'SPAD' (confocal recording)
                        % 'EMCCD' (camera recording)
                        % 'none' (effect-less ideal recording)
gamma = 1; % experimental gamma factor

% For more control on parameters, check:
%   src/set_FRET_hist_fit_param.m
%   src/get_inference_param.m
% -------------------------------------------------------------------------

% ask user for histogram text file
[fle,src] = uigetfile('*.mat','Select an histogram file.',pwd);
if ~sum(src)
    return
end

% set MATLAB search path
codePath = fileparts(mfilename('fullpath'));
addpath(genpath(codePath));

% collect inference parameters
[flx_grid,opt] = set_FRET_hist_fit_param(recordingtype, gamma,...
    inference_method);

% import data from file
[outfle0,F_hist,P_table,fprm] = collect_exp_FRET_hist(src,fle);
[dest,fnm,~] = fileparts(outfle0);

% create folder for DPMM junk runs
if strcmp(opt.inference_method,'DPMM')
    dest_junk = [dest,filesep,'DPMM-junk-runs'];
    if ~exist(dest_junk,'dir')
        mkdir(dest_junk)
    end
end

% pre_calculate PMF table if necessary
[P_table,E_grid,F_hist] = check_P_table(P_table,F_hist,opt.Nsim,fprm,...
    flx_grid,opt.const,[src,filesep,fle]);

% cluster FRET histogram with mixture of discrete distributions
[res,plt] = FRET_hist_discrete_fit(F_hist,P_table,E_grid,flx_grid,opt);

if opt.plotiter || opt.plotfinal
    % export figure to .fig and .png file in light theme
    if ~iscell(plt)
        plt = {plt};
    end
    for i = 1:numel(plt) % one figure exists for each DPMM run
        outfle = outfle0;
        if strcmp(opt.inference_method,'DPMM')
            if res.bestrun~=i
                outfle = [dest_junk,filesep,fnm,sprintf('_run%i.png',i)];
            end
        end
        export_fig_to_file(plt{i}.fig,outfle,opt.inference_method);
        close(plt{i}.fig);
    end
end

% save results to file
save(addnbtofilename(...
    [dest,filesep,fnm,'_results_',opt.inference_method,'.mat']),'res');

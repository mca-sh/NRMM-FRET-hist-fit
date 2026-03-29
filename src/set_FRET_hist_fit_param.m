function [flx_grid,opt] = set_FRET_hist_fit_param(recordingtype,gamma,...
    inference_method)
% [flx_grid,opt] = set_FRET_hist_fit_param(recording, gamma, method)
%
% Define parameters for histogram peak model and inference method.
%
% -- INPUT: ---------------------------------------------------------------
% recordingtype: Recording support 'SAPD' for single-photon avalanche
%                photo diod, 'EMCCD' for electron-multiplying camera, 
%                'none' for ideal recording.
% gamma: experimental gamma factor.
% method: inferrence method 'DPMM' for non-parametric inference of 
%         noise-realist mixture models, 'EM' for expectation-maximization 
%         of noise-realist mxture models, or 'EM-GMM' for expectation-
%         maximization of Gaussian mixture models.
%
% --- OUTPUT: -------------------------------------------------------------
% flx_grid: grid values of structure flexibility used in table of 
%           pre-calculated PMF.
% opt: structure containing model and inference parameters.
%
% -------------------------------------------------------------------------

% defaults
opt.Nsim = 1E5;                  % nb. of simulation for one cell of PMF grid
opt.plotiter = false;            % update plot after each iteration
opt.plotfinal = true;            % update plot after inference

% save information about user's data and method
opt.inference_method = inference_method; % "EM", "EM-GMM" or "DPMM"
const.gamma = gamma;                     % gamma factor

% build flexibility grid
flx_edg = linspace(0,1,50);
flx_grid = mean([flx_edg(1:end-1);flx_edg(2:end)],1);

% FRET distance distribution constants
const.R0 = 57; % (in A) Forster radius
const.k = 2;   % (in per A) Sigmoïde slope factor for flex. modelling
const.S = 2;   % (in A) Sigmoïde tipping point for flex. modelling

% noise constants
const.q = 0.95;  % Detection efficiency
const.npix = 1;  % nb. of pixels integrated
const.g = 1;     % EM gain
const.f = 1;     % analog-to-digital factor
const.b = 0;     % camera offset (IC/pixel)
const.bgA = 0;   % background light in acceptor channel (PC/pixel)
const.bgD = 0;   % background light in donnor channel (PC/pixel)
const.I0 = 100;  % total donnor emission (PC/pixel)
switch recordingtype
    case 'SPAD'
        % confocal + APD noise parameters (photon counts from file
        % 20240405_eD135e_L43E_100mMMg_pre_4_APBS_2CnoMFD.bur)
        const.noisetype = 'P';
        I0_mean = 71.74; % log(2)*median(GG+GR) using 50PC min threshold
        thresh_I0_min = 50;
        num_I0 = 0;
        const.I0 = [];
        while num_I0<opt.Nsim
            const.I0 = cat(2,const.I0,...
                exprnd(I0_mean,[1,opt.Nsim-num_I0]));
            const.I0(const.I0<thresh_I0_min) = [];
            num_I0 = length(const.I0);
        end
        const.bgD = 0.04*I0_mean;
        const.bgA = 0.04*I0_mean;

    case 'none'
        % for discrete FRET histograms
        const.noisetype = 'none';

    case 'EMCCD'
        % EMCCD noise parameters determined from MH45-FRET-data.mat
        const.noisetype = 'PGN';
        const.c = 0.02;
        const.g = 300;
        const.f = 5.199;
        const.b = 113;
        const.npix = 8;
        const.r = 0.067;
        const.I0 = ic2pc(10000/const.npix + const.b,const); % exp mean Itot=10000 ic
        const.bgD = ic2pc(2500/const.npix,const); % exp bg=2500 ic
        const.bgA = ic2pc(1340/const.npix,const);  % exp bg=1340 ic

    otherwise
        disp('recording type must be ''SPAD'', ''none'' or ''EMCCD''.');
        return
end

% build parameter structure of inference method
[prm,opt.Nmax] = get_inference_param(opt.inference_method);

opt.const = const;
opt.prm = prm;

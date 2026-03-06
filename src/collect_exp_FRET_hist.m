function [outfle,F_hist,P_table,fprm] = collect_exp_FRET_hist(src,fle)
% [outfle,F_hist,P_table,fprm] = collect_exp_FRET_hist(src,fle)
%
% Import FRET histogram and metadata from experimental data file

% create analysis forlder if not existing
[~,fnm,fext] = fileparts(fle);
dest = cat(2,src,filesep,fnm);
if ~exist(dest,'dir')
    mkdir(dest);
end
outfle = [dest,filesep,fnm,fext];

% turn off warning about a variable not found in .mat file
warning('off','MATLAB:load:variableNotFound');

% import FRET histogram and pre-calculated PMF grid
load([src,filesep,fle],'F_hist','P_table','flx_grid0','E_grid0',...
    'xbins0','const0');
if ~exist('P_table','var')
    P_table = [];
end

% build structure with file parameters
if exist('E_grid0','var')
    fprm.E_grid0 = E_grid0;
else
    fprm.E_grid0 = [];
end
if exist('flx_grid0','var')
    fprm.flx_grid0 = flx_grid0;
else
    fprm.flx_grid0 = [];
end
if exist('xbins0','var')
    fprm.xbins0 = xbins0;
else
    fprm.xbins0 = [];
end
if exist('const0','var')
    fprm.const0 = const0;
else
    fprm.const0 = [];
end

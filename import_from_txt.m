% --- README --------------------------------------------------------------
% IMPORT_FROM_TXT format histogram data contained in a text file to a 
% .mat file suitable for analysis with "script_fit_hist.m".
%
% The 1st column  of the text file must contains the centers of the 
% histogram bins and the second, the histogram counts.
%
% The .mat file is exported in the source directory.
% -------------------------------------------------------------------------

% ask user for histogram text file
[fle_txt,src_txt] = uigetfile('*','Select an histogram file.',pwd);
if ~sum(src_txt)
    return
end

% extract histogram data and save to .mat analysis file
txt_to_mat_hist_file([src_txt,filesep,fle_txt]);


function fle_mat = txt_to_mat_hist_file(fle_txt)
    % fle_mat = txt_to_mat_hist_file(fle_txt)
    %
    % Read histogram data from an ASCII file and export them to a .mat file 
    % suitable for analysis with FRET-hist-fit.
    % ASCII file must contains at least two columns whith the first column
    % listing the centers of the histogram bins and the second the 
    % histogram counts.
    %
    % --- INPUT: ----------------------------------------------------------
    % fle_txt: path to source ASCII file.
    %
    % --- OUTPUT: ---------------------------------------------------------
    % fle_mat: path to destination MATLAB file.
    %
    % ---------------------------------------------------------------------
    
    % collect histogram data from text file and format
    fdat = importdata(fle_txt);
    F_hist = fdat.data(:,[1,2])';
    
    % save to .mat file
    [src,fname,~] = fileparts(fle_txt);
    fle_mat = [src,filesep,fname,'.mat'];
    save(fle_mat,'F_hist');

    disp(['Histogram was successfully saved to file ',fle_mat, ' !']);
end
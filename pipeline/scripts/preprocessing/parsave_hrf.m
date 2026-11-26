function parsave_hrf(mat_path, para, beta_hrf, hrfa, PARA, event_bold)
%PARSAVE_HRF  Save HRF estimation variables inside a parfor loop.
%   Wrapping save() in a function avoids the workspace transparency
%   violation that occurs when save() with string variable names is
%   called directly inside parfor.
    save(mat_path, 'para', 'beta_hrf', 'hrfa', 'PARA', 'event_bold', '-v7.3');
end

function out_st = eeglab_ica(data_st,params)
% Ensemble wrapper script for using EEGLAB to perform an ICA analysis on the
% data contained in the input data structure.
%
% The input data structure should contain an EEG variable which either contains
% the filenames of .set files that contain EEG structures that serve as input,
% or contains the EEG structures themselves. 

% 12/03/07 Petr Janata

% Make sure that EEGLAB variables and paths have been initialized. This will
% pop-up an EEGLAB window if none currently exists
eeglab('nogui')

out_st = ensemble_init_data_struct;

% Perform any filtering that the calling function wants to have performed
if isfield(params, 'filt')
  data_st = ensemble_filter(data_st, params.filt);
end

% Get the column ids
cols = set_var_col_const(data_st.vars);

% Get a list of unique subjects
subids = unique(data_st.data{cols.subject_id});
nsub = length(subids);

fprintf('%s: Processing data for %d subjects\n', mfilename, nsub);

% Loop over subjects
for isub = 1:nsub
  subid = subids{isub};
  fprintf('%s: Subject (%d/%d): %s\n', mfilename,isub,nsub, subid);

  % Deal with finding the file to process
  subject_path = fullfile(params.path.project_root,subid);
  set_path = fullfile(subject_path,'set');
  try fstub = params.ica.fstub; catch fstub = ''; end
  
  fname = fullfile(set_path, [subid fstub '.set']);
  if ~exist(fname)
    fprintf('Could not find file: %s\n', fname);
    continue
  end
  
  try 
    EEG = pop_loadset('filename',[subid fstub '.set'],'filepath',set_path);
    
    try max_samps = params.ica.max_samps; catch max_samps = EEG.pnts; end
    
    start_samp = 1;
    stop_samp = min(EEG.pnts, max_samps);
    tmpdata = EEG.data;
    EEG.data = EEG.data(:,start_samp:stop_samp);
    EEG.pnts = stop_samp-start_samp+1;
    keyvals = {};
    check_ica_params = {'lrate'};
    for ip = 1:length(check_ica_params)
      cparam = check_ica_params{ip};
      if isfield(params.ica, cparam)
	keyvals{end+1} = cparam;
	keyvals{end+1} = params.ica.(cparam);
      end
    end
    EEG = pop_runica(EEG,'runica',keyvals{:});
    
    % Paste the original EEG data back in
    EEG.data = tmpdata;
    EEG.pnts = size(tmpdata,2);
    
    % Check and save the ICA information to the same file
    EEG = eeg_checkset(EEG);
    EEG = pop_saveset (EEG,'savemode','resave');
  catch
    fprintf('Failed to complete ICA analysis for subject: %s\n', subid);
  end
end % for isub=

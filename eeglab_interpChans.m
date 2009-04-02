function out_st = eeglab_interpChans(data_st,params)
% Ensemble wrapper script for interpolating data at missing channels 
%
% Interpolations are currently spherical splines interpolations using Tom
% Ferree's code

% 03/30/09 PJ - adapted from eeglab_rmbadchans

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

% Load channel information
locfile = params.eeglab.chanlocs.locfile;
fprintf('Attaching electrode locations from file: %s\n', locfile);
EEGinterp.chanlocs = readlocs(locfile);
nchanlocs = length(EEGinterp.chanlocs);

%
% Deal with removing (un)wanted channels.
% The use_chans field lists channels to keep.
try use_chans = params.bdf.use_chans; catch use_chans = []; end
if ~isempty(use_chans)
  EEGinterp = remove_chans(EEGinterp,use_chans,'keep');
end


% Get a list of unique subjects
subids = unique(data_st.data{cols.subject_id});
nsub = length(subids);

fprintf('%s: Processing data for %d subjects\n', mfilename, nsub);

% Loop over subjects
for isub = 1:nsub
  subid = subids{isub};
  fprintf('\n%s: Subject (%d/%d): %s\n', mfilename,isub,nsub, subid);

  % Deal with finding the file to process
  subject_path = fullfile(params.path.project_root,subid);
  set_path = fullfile(subject_path,'set');
  try srcfstub = params.srcfstub; catch srcfstub = ''; end
  try destfstub = params.destfstub; catch destfstub = ''; end
  
  fname = fullfile(set_path, [subid srcfstub '.set']);
  if ~exist(fname)
    fprintf('Could not find file: %s\n', fname);
    continue
  end
  
  %try 
    % Find the subject info in the sinfo structure
    sidx = strmatch(subid,{params.sinfo.id});
    
    if isempty(sidx)
      error(sprintf('Subject %s not found in sinfo', subid))
    end
        
    % Load the EEG data
    EEG = pop_loadset('filename',[subid srcfstub '.set'],'filepath',set_path);
    
    % See if we are writing to a new file
    if strcmp(srcfstub,destfstub)
      savemode = 'resave';
    else
      savemode = 'twofiles';
    end
    
    EEG.setname = sprintf('%s, spherical spline interpolation', EEG.setname);

    % Perform the interpolation
    EEG = EEG_transform(EEG,EEGinterp.chanlocs);
    
    % Check and save the information to the destination file. Might be same
    % as source file
    EEG = eeg_checkset(EEG);
    destfname = [subid destfstub '.set'];
    fprintf('Saving data to: %s\n', fullfile(set_path, destfname));
    EEG = pop_saveset (EEG,'filepath',set_path, ...
      'filename',destfname, 'savemode', savemode);
    try
    catch
      fprintf('Failed to remove bad channels for subject: %s\n', subid);
    end
end % for isub=

end

function EEG = remove_chans(EEG,chanlist,action)
  if nargin < 3
    action = 'remove';
  end
  
  switch action
    case 'remove'
      retain_chan_idxs = find(~ismember({EEG.chanlocs.labels}, chanlist));
    case 'keep'
      retain_chan_idxs = find(ismember({EEG.chanlocs.labels}, chanlist));
  end
  
  if ~isempty(retain_chan_idxs)
    fprintf('Retaining %d channels ...\n', length(retain_chan_idxs));
    try EEG.data(setdiff(1:size(EEG.data,1),retain_chan_idxs),:) = []; catch EEG.data = []; end
    EEG.chanlocs =EEG.chanlocs(retain_chan_idxs);
    EEG.nbchan = size(EEG.data,1);
  end
end % remove_chans


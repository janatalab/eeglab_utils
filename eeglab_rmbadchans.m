function out_st = eeglab_rmbadchans(data_st,params)
% Ensemble wrapper script for using EEGLAB to remove bad channels from the
% data contained in the input data structure.
%

% 02/22/09 PJ - adapted from eeglab_ica
%
% 02/14/13 BH - added some code to handel differences in parameter naming 
%               conventions between current & previous projects
%             - added parallel computing code for looping over subjects

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

% Fire up parallel compute pool if possible
if exist('matlabpool') && ~matlabpool('size') 
    matlabpool
end

% Loop over subjects in parallel
parfor isub = 1:nsub
  subid = subids{isub};
  fprintf('\n%s: Subject (%d/%d): %s\n', mfilename,isub,nsub, subid);

  % Deal with finding the file to process
  if isfield(params.path,'datapath')
      datapath = params.path.datapath;
  else
      datapath = params.path.project_root;
  end
  subject_path = fullfile(datapath,subid);
  set_path = fullfile(subject_path,'set');
  try srcfstub = params.srcfstub; catch srcfstub = ''; end
  try destfstub = params.destfstub; catch destfstub = ''; end
  
  fname = fullfile(set_path, [subid srcfstub '.set']);
  if ~exist(fname)
    fprintf('Could not find file: %s\n', fname);
    continue
  end
  
  if isfield(params.sinfo,'subject_id')
      sinfo_subids = {params.sinfo.subject_id};
  else 
      sinfo_subids = {params.sinfo.id};
  end
  
  try 
    % Find the subject info in the sinfo structure
    sidx = strmatch(subid,sinfo_subids);
    
    if isempty(sidx)
      error(sprintf('Subject %s not found in sinfo', subid))
    end
    
    bad_chans = params.sinfo(sidx).bad_chans;
    
    if isempty(bad_chans)
      fprintf('No bad channels for subject\n')
      continue
    end
    
    % Just load the file information but not the data so that we can
    % determine if the channels have already been removed
    EEG = pop_loadset('filename',[subid srcfstub '.set'],'filepath',set_path, ...
      'loadmode', 'info');
    
    % Determine the channel indices
    rm_idxs = find(ismember({EEG.chanlocs.labels}, bad_chans));
    if isempty(rm_idxs)
      fprintf('Bad channels already removed\n')
      continue
    end

    % If we are here, then deal with loading the data
    EEG = pop_loadset('filename',[subid srcfstub '.set'],'filepath',set_path);
    
    % remove the channels
    EEG = pop_select(EEG,'nochannel', rm_idxs);

    % See if we are writing to a new file
    if strcmp(srcfstub,destfstub)
      savemode = 'resave';
    else
      savemode = 'twofiles'
    end
    
    EEG.setname = sprintf('%s, -%d bad chans', EEG.setname,length(bad_chans));

    try EEG.etc.badchans = [EEG.etc.badchans bad_chans]; 
    catch EEG.etc.badchans = bad_chans;
    end
    
    % Check and save the information to the destination file. Might be same
    % as source file
    EEG = eeg_checkset(EEG);
    destfname = [subid destfstub '.set'];
    fprintf('Saving data to: %s\n', fullfile(set_path, destfname));
    EEG = pop_saveset (EEG,'filepath',set_path, ...
      'filename',destfname, 'savemode', savemode);
  catch
    fprintf('Failed to remove bad channels for subject: %s\n', subid);
  end
end % parfor isub=

matlabpool close
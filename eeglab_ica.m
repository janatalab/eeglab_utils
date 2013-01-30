function out_st = eeglab_ica(data_st,params)
% Ensemble wrapper script for using EEGLAB to perform an ICA analysis on the
% data contained in the input data structure.
%
% The input data structure should contain an EEG variable which either contains
% the filenames of .set files that contain EEG structures that serve as input,
% or contains the EEG structures themselves. 

% 12/03/07 Petr Janata
% 02/22/09 PJ - Improved subsampling of original data in the event that the
% number of points in the data exceed maximum allowed for ICA.

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
	if isfield(params.path, 'eegpath')
		eegpath = params.path.eegpath;
	else
		eegpath = params.path.project_root;
	end
  subject_path = fullfile(eegpath,subid);
  set_path = fullfile(subject_path,'set');
  try fstub = params.ica.fstub; catch fstub = ''; end
  
  fname = fullfile(set_path, [subid fstub '.set']);
  if ~exist(fname)
    fprintf('Could not find file: %s\n', fname);
    continue
  end
  
  try 
    EEG = pop_loadset('filename',[subid fstub '.set'],'filepath',set_path);
    
    % Make sure channel data is appropriately attached
    EEG = eeglab_check_chaninfo(EEG,params);
    
    % Deal with the situation in which we need/want to subdivide the data
    % that we train on
    try max_samps = params.ica.max_samps; catch max_samps = EEG.pnts; end
    try num_segments = params.ica.num_segments; catch num_segments = 1; end
    
     if max_samps < EEG.pnts
       subsample = true;
     else
       subsample = false;
     end
    
    if subsample
      fprintf('Number of points in data exceeds limit for ICA algorithm\n')
      fprintf('Sampling data from %d segments\n', num_segments);
      
      tmpdata = EEG.data; % make a copy of the data
      EEG.data = []; % clear out the existing data
      segsize_in_orig = fix(EEG.pnts/num_segments); % size of subdivisions in original data
      max_samps_per_seg = fix(max_samps/num_segments);
      fprintf('%d samples per segment\n', max_samps_per_seg);
      
      fprintf('Segment\tStartSamp\tStopSamp\n');
      for iseg = 1:num_segments
        % Figure out range of possible start samples for this segment
        min_seg_start = (iseg-1)*segsize_in_orig+1;
        max_seg_start = iseg*segsize_in_orig-max_samps_per_seg+1;
       
        % Choose a random starting location within the current segment
        start_samp = min_seg_start+fix(rand*(max_seg_start-min_seg_start));
        stop_samp = start_samp+max_samps_per_seg-1;

        fprintf('%d\t%d\t%d\n', iseg, start_samp, stop_samp);
        
        % Add the current chunk to the larger chunk, making sure to remove
        % the local baseline
        EEG.data = [EEG.data rmbase(tmpdata(:,start_samp:stop_samp))];
        
      end % for iseg
      EEG.pnts = size(EEG.data,2);      
    end % if max_samps > EEG.pnts && num_segments > 1
  
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
    if subsample
      EEG.data = tmpdata;
      EEG.pnts = size(tmpdata,2);
    end
    
    % Check and save the ICA information to the same file
    try destfstub = params.destfstub; catch destfstub = ''; end

    destfname = [subid destfstub '.set'];
    EEG = eeg_checkset(EEG);
    EEG = pop_saveset (EEG,'filepath',set_path,'filename',destfname,'savemode','twofiles');
  catch
    fprintf('Failed to complete ICA analysis for subject: %s\n', subid);
  end
end % for isub=

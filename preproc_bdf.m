function out_st = preproc_bdf(data_st,params)
% Performs preprocessing of data stored in Biosemi BDF files using EEGLAB
% routines. Specifically, this routine will:
%
% 1) Select or remove channels specified in params.bdf.use_chans and
%    params.bdf.remove_chans, respectively, using the remove_chans()
%    subfunction.
%
% 2) Attach electrode location information specified in
%    params.eeglab.chanlocs.locfile. Uses pop_chanedit()
% 
% 3) Rereference data to be average reference if desired. Default is to not
%    rereference. Uses pop_reref()
%
% 4) Removes DC offset for each channel. Default is to remove it.
%
% 5) Low-pass and high-pass filtering of each channel. Filtering is applied
%    in two separate calls to eegfilt().
%
% 6) If data are SCR data, convert to microsiemens.
%
% 7) Merges datasets if multiple BDF files were encountered and they were
%    designated to be merged.
%
% The input data structure (data_st) should contain a variable that contains
% all of the subject IDs to be processed.
%
% If preproc_bdf is called with 'get_default_params', a params structure is
% returned that lists all of the fields that can be set in order to control the
% behavior of this function.
%
% NOTE: If you encounter strange problems reading event data correctly from
% the BDF files, it might be due to using a problematic version of sopen()
% which is part of the Biosig routines. To solve this problem, add
% '~/svn/public/matlab/eeg/fixes/' to your path

% 05/20/07 Petr Janata - started script
% 07/29/08 Petr Janata - added support for finding bdf directories buried
%                        within session directories 
% 05/12/09 Fred Barrett - attaching bdf file header (returned by sopen) to
% the EEG data structure, in EEG.etc.bdfhdr. This contains information that
% may be useful or vital for later processing and interpretation of SCR
% (and other) signals (for instance, bdfhdr.PhysDim)
% 05/13/09 Fred Barrett - added mechanism to convert SCR signal from Ohms or
% nanosiemens to microsiemens. Older versions of BioSemi hardware save SCR
% data in Ohms, newer versions save SCR data in nanosiemens, but the standard
% dimension with which SCR is discussed is microsiemens (or micromhos)
% 12/11/09 PJ - added dynamic handling of adding biosig directories to path
% 12/13/09 PJ - rearranged code to accommodate recent changes in pop_reref()
% 24Jan2013 PJ - added option to detrend the EEG data.
%

% Initialize EEGLAB paths
eeglab('nogui')

% Need to add paths for eeglab/external/biosig subdirectories. These include
% critical biosig functions such as sopen(). Ideally calling
% eeglab('nogui') would add these paths
p = which('eeglab.m');
p = p(1:findstr(p,'eeglab.m')-1);
extlist = dir(fullfile(p,'external')); % list external directory
biosig_idx = strmatch('biosig',{extlist.name});
biosig_dir = fullfile(p,'external',extlist(biosig_idx).name);
path(genpath(biosig_dir),path);

if ischar(data_st) && any(ismember(data_st,{'get_default_params','getDefaultParams'}))
  out_st = get_default_params;
  return
end

out_st = ensemble_init_data_struct;
out_st.vars = {'subject_id','session_id','filenum','EEG','origfname'};
out_st.data = cell(1,length(out_st.vars));
outcols = set_var_col_const(out_st.vars);

% Perform any filtering that the calling function wants to have performed
if isfield(params, 'filt')
  data_st = ensemble_filter(data_st, params.filt);
end

% Get the column ids
cols = set_var_col_const(data_st.vars);

% Get a list of unique subjects
subids = unique(data_st.data{cols.subject_id});
nsub = length(subids);

% Probably need some handling here for situations in which we don't have a
% separate sinfo data structure that has subject specific information about bad
% channels, etc.
try save_set = params.eeglab.save.set; catch save_set = true; end
try merge_sets = params.eeglab.merge_sets; catch merge_sets = false; end
try return_eeg = params.eeglab.return.EEG; catch return_eeg = false; end

% Loop over subjects
fprintf('%s: Pre-processing data for %d subjects\n', mfilename,nsub);
for isub = 1:nsub
  clear EEG ALLEEG
  
  subid = subids{isub};
  fprintf('\nProcessing subject (%d/%d): %s\n', isub, nsub, subid);
	
	sinfoIdx = find(strcmp(subid, {params.sinfo.subject_id}));
	curr_sinfo = params.sinfo(sinfoIdx);
  % Initialize the list of bdf_paths for this subject that we are going to process
  bdf_paths = {};

  % Initialize the list of session directories containing each of the
  % corresponding bdf files in bdf_paths
  bdf_paths_session_map = {};
	
	% Construct a default subject path
	if isfield(params.paths,'eegpath')
		eegpath = params.paths.eegpath;
	else
		eegpath = params.paths.project_root;
	end
  subject_path = fullfile(eegpath,subid);
 	
	% Make sure we have a valid subject_path
	if ~exist(subject_path,'dir')
		fprintf('Cannot find default subject directory: %s\n', subject_path)
		
		% See if a sinfo structure was passed in
		if ~isfield(params,'sinfo')
			error('Do not have sinfo struct that might contain alternative subject ID mappings')
		end
		
		% Get list of fields
		sinfoFields = fieldnames(curr_sinfo);
		
		% Find an alternate path
		while ~exist(subject_path,'dir')
			% Check a mapping from ensemble_id to subject_id
			match_flds = {'ensemble_id', 'subject_id'};
			haveFieldsMask = ismember(sinfoFields, match_flds);
			if sum(haveFieldsMask) == length(match_flds);
				subid = curr_sinfo.subject_id;
				fprintf('Re-mapped subject ID. New subid: %s\n', subid);
				subject_path = fullfile(eegpath,subid);
				continue
			end
			
			error('Found no alternative subject ID mappings')
		end % while ~exist(subject_path,'dir')
	end % ~exist(subject_path,'dir')

	fprintf('Using subject_path: %s\n', subject_path);
	
  % Read the BDF file
  bdf_prefix = sprintf('%s%s', params.bdf.file_prefix, subid);
  if exist(fullfile(subject_path,'bdf'))
    bdf_paths = {fullfile(subject_path,'bdf')};
  end
  
  % Check to see if bdf directory exists at this level or if it is buried in
  % session directories

  % Check the sinfo session variable to figure out which sessions to use
  if isfield(curr_sinfo, 'use_session')
      use_sessions = curr_sinfo.use_session;
  else
      use_sessions = []; 
  end

  if isempty(bdf_paths)
    fprintf(['Could not locate bdf directory in subject directory: %s\n' ...
      'Checking for session directories ...\n'], subject_path);
    sesslist = dir(fullfile(subject_path,'session*'));
    if ~isempty(sesslist)
      fprintf('Found %d session directories\n', length(sesslist));
      
      
      % If no info about which session to use was found, check to see if
      % there is only one session and use that by default
      if isempty(use_sessions)
        fprintf(['Found no session information in sinfo structure\n' ...
          'Will use session1 if that is the only session directory\n']);
        if length(sesslist) == 1
          use_sessions = 1;
        end
      end % if isempty(use_sessions)
		else
			use_sessions = [];
    end
    
    % Construct the list of bdf_paths

    for isess = 1:length(use_sessions)
        session_stub = sprintf('session%d', use_sessions(isess));
      tmppath = fullfile(subject_path, session_stub, 'bdf');
      if exist(tmppath)
        bdf_paths{isess} = tmppath;
        bdf_paths_session_map{isess} = session_stub;
      end
    end % for isess
    
    % Make sure we have at least one directory that is a bdf directory
    if isempty(bdf_paths)
      error('Could not locate any bdf directories for this subject\n');
    end
    
  end % if isempty(bdf_paths)
  
  %
  % Loop over bdf directories
  %
  nbdf_dir = length(bdf_paths);
  fprintf('Processing data in %d BDF directories\n', nbdf_dir);
  for idir = 1:nbdf_dir
    bdf_path = bdf_paths{idir};
    
    % Figure out how many BDF files there are in this directory. There could be
    % multiple files due to multiple runs.
    bdfname = sprintf('*%s*.bdf',params.srcfstub) ;
    if params.sufix
      bdflist = dir(fullfile(bdf_path,bdfname));%'*%s.bdf');
    else
      bdflist = dir(fullfile(bdf_path,'*.bdf'));
    end
    
    num_bdf_files = length(bdflist);
    fprintf('Found %d BDF files in %s\n', num_bdf_files, bdf_path);
    
    if ~num_bdf_files
      continue
    end
    
    for ifile = 1:num_bdf_files
      bdffname = fullfile(bdf_path,bdflist(ifile).name);
      [dummy,orig_fstub] = fileparts(bdffname);
      
      % Set up some variables for saving the resulting .set file
      if ~merge_sets
        set_stub = sprintf('%s', orig_fstub);
      else
        set_stub = sprintf('%s_allbdf', subid);
      end

      % Determine the directory into which we'll put the resulting .set
      % file
      if ~isempty(use_sessions) && params.eeglab.use_session_dirs
        set_path = fullfile(subject_path, bdf_paths_session_map{isess}, 'set');
      else
        set_path = fullfile(subject_path,'set'); 
      end
      check_dir(set_path); % make sure the output directory exists
      %curr_fstub = set_stub;
      
      % Read the BDF header
      fprintf('Reading BDF header from file: %s\n', bdffname);
      bdfhdr = sopen(bdffname);
      sclose(bdfhdr);
      
      
      fprintf('Loading BDF data using pop_biosig from file: %s\n', bdffname);
      
      EEG = pop_biosig(bdffname, ...
        'rmeventchan', params.bdf.rmeventchan, ...
        'blockrange', params.bdf.blockrange, ...
        'ref', params.bdf.refchan);
      
      % attach bdf header structure to EEG.etc
      EEG.etc.bdfhdr = bdfhdr;
      
			%
			% Check whether we want to re-read the event information correctly
			% (bits 9-16)
			%
			try reread_events = params.eeglab.fix_eeg_event_info; catch reread_events = 0; end
			if reread_events
				EEG = fix_eeg_event_info(EEG);
			end
			
      %
      % Deal with removing (un)wanted channels.
      % The use_chans field lists channels to keep.  This pretty much needs
      % to include the full set of EEG channels listed in the electrode
      % location information for the electrode location information to load
      % properly
			origChans = EEG.chanlocs;
      try use_chans = params.bdf.use_chans; catch use_chans = []; end
      if ~isempty(use_chans)
        EEG = remove_chans(EEG,use_chans,'keep');
			end
			
      %
      % Attach the electrode location information
      %
      % NOTE: We need to make this part of it a bit smarter to handle sets of
      % channels that don't necessarily match the chanlocs file exactly
      try attach_locs = params.eeglab.chanlocs.locfile; catch attach_locs = ...
          false; end
      if attach_locs
        locfile = params.eeglab.chanlocs.locfile;
        fprintf('Attaching electrode locations from file: %s\n', locfile);
				
        EEG = pop_chanedit(EEG,'load',locfile);
        nchanlocs = length(EEG.chanlocs);

        % Make sure that the number of channels for which we have labels matches the
        % actual number of channels
        if EEG.nbchan ~= nchanlocs
          warnmsg = sprintf(['Mismatch between number of channel locations ' ...
            'in data (%d) and location spec file(%d)\n'], EEG.nbchan, nchanlocs);
          warning(warnmsg)
				end
				
				% Deblank channel labels
				for ichan = 1:length(EEG.chanlocs)
					EEG.chanlocs(ichan).labels = strrep(EEG.chanlocs(ichan).labels, ' ', '');
				end
      end % if attach_locs
      
      % Remove channels that were inadvertently recorded but are empty.  The list
      % of these is stored in defs.eeg.extra_chans.  Match them up by comparing
      % chan labels read from BDF header
      try toss_chans = params.bdf.remove_chans; catch toss_chans = []; end
      if ~isempty(toss_chans)
        EEG = remove_chans(EEG, toss_chans,'remove');
      end
      
      % Remove channels specified as bad for this individual subject
      try bad_chans = curr_sinfo.bad_chans; catch bad_chans = []; end

       % Check whether any session-specific bad channels have been
       % specified
       if isempty(bad_chans) 

       end


      if ~isempty(bad_chans)
           EEG = remove_chans(EEG, bad_chans,'remove');
			end
          


      % Re-reference the data to an average reference
      % overwrite unreferenced dataset?
      try reref = params.eeglab.reref; catch reref = false; end
      if reref
        fprintf('Re-referencing the data to be average reference\n');
        
        % The refchan that was used when the BDF data were read in now
        % needs to be ignored, otherwise it is duplicated in the data. In
        % 2009, pop_reref() in eeglab was changed which prompted the change
        % in the code here.  The old code is simply commented out for
        % reference.
%         if ~isempty(params.bdf.refchan)
%             refloc{1} = EEG.chanlocs(params.bdf.refchan).labels;
%             reflEEG = pop_chanedit(EEG,'load',locfile,oc{2} = EEG.chanlocs(params.bdf.refchan).theta;
%             refloc{3} = EEG.chanlocs(params.bdf.refchan).radius;
%             EEG = pop_reref( EEG, [], 'refloc', refloc);
%         else
%           EEG = pop_reref( EEG, []);
%         end
%         
        EEG = pop_reref( EEG, []);
       
      end
      
      % See if we want to detrend the data
      if isfield(params.eeglab, 'detrend')
          detrend_eeg = params.eeglab.detrend;
      else
          detrend_eeg = false;
      end
      
      % Remove the offset of each channel
      try remove_dc = params.eeglab.remove_dc; catch remove_dc = true; end
      if remove_dc || detrend_eeg
          if detrend_eeg
              str = '';
              dtype = 'linear trend';
          else
              str = 'constant';
              dtype = 'offset';
          end
          
        fprintf('Removing the %s from each channel ...\n', dtype);
        EEG.data = single(detrend(double(EEG.data'),str))';        
      end
      
      %
      % Perform channel-type specific filtering if specified.
      % The filtering is done with eegfilt. Note that because we are filtering
      % different channel types with different filter settings, we have to call
      % eegfilt() directly, rather than pop_eegfilt().  The latter would send all
      % of the data channels.  The only issue here is that the data MUST BE
      % CONTINUOUS, i.e. it musn't have been segmented so as to create
      % boundaries/discontinuities.
      %
      try overwrite_filtered = params.eeglab.save.clobber; catch ...
          overwrite_filtered = 0; end
      typelist = fieldnames(params.eeglab.filter);
      ntypes = length(typelist);

			% Construct the putative name for the output file and check to see whether
      % it exists.
      filt_stub = sprintf('_filt');
      curr_fname = sprintf('%s%s.set', set_stub,filt_stub);
      filtfname = fullfile(set_path, curr_fname);
      
      % Check to see if the file already exists
      if ~exist(filtfname,'file') || overwrite_filtered
        if exist(filtfname,'file') && overwrite_filtered
          fprintf('Overwriting existing file: %s\n', filtfname);
        end
        
        for itype = 1:ntypes
          chantype = typelist{itype};
          
          % Make sure we have channels that fit the criteria
          chan_idxs = ...
            find(ismember({EEG.chanlocs.labels},params.eeglab.filter.(chantype).chan_labels));
          
          try low_cutoff = params.eeglab.filter.(chantype).low_cutoff;
          catch low_cutoff = []; end
          
          try high_cutoff = params.eeglab.filter.(chantype).high_cutoff;
          catch high_cutoff = []; end
          
          if (~isempty(low_cutoff) || ~isempty(high_cutoff)) && ~isempty(chan_idxs)
            
            fprintf('Filtering %d channels of type: %s\n', length(chan_idxs), upper(chantype));
            try filt_order = params.eeglab.filter.(chantype).filt_order; catch filt_order = []; ...
            end
          
          % Filter data  - do highpass and low-pass in separate stages to avoid
          % numerical problems
          if low_cutoff
            if isempty(filt_order)
              try filt_order = params.eeglab.filter.(chantype).filt_order_multiplier*fix(EEG.srate/low_cutoff);
              catch
                filt_order = [];
              end
            end
            
            EEG.data(chan_idxs,:) = eegfilt(EEG.data(chan_idxs,:), EEG.srate, ...
              low_cutoff, 0, 0, filt_order);
          end
          if high_cutoff
            EEG.data(chan_idxs,:) = eegfilt(EEG.data(chan_idxs,:), EEG.srate, ...
              0, high_cutoff, 0, []);
          end
          end % if low_cutoff | high_cutoff
          
          % convert data to microsiemens? useful for SCR data
          try c2mS = params.eeglab.filter.(chantype).convert2microsiemens;
          catch c2mS = 0; end
          
          if c2mS
            %%% convert to microsiemens??
            %%% the new BioSemi SCR system saves data in nanosiemens
            %%% the old BioSemi SCR system saves data in Ohms
            
            % iterate over channels in chan_idxs
            for ilci = 1:length(chan_idxs)
              % Check physical dimension of the current channel
              lchan = chan_idxs(ilci);
              etcidx = strmatch(EEG.chanlocs(lchan).labels,...
                EEG.etc.bdfhdr.Label);
              if isempty(etcidx), continue, end
              if ~isempty(strmatch('nS',EEG.etc.bdfhdr.PhysDim(etcidx),...
                  'exact'))
                % convert from nS to mS
                % mS = nS./1000
                EEG.data(lchan,:) = EEG.data(lchan,:)./10^3;
              elseif ~isempty(strmatch('Ohm',...
                  EEG.etc.bdfhdr.PhysDim(etcidx),'exact'))
                % convert from Ohm to micromho (=mS)
                % mS = (1./Ohm)*1,000,000
                EEG.data(lchan,:) = 1./EEG.data(lchan,:);
                EEG.data(lchan,:) = EEG.data(lchan,:).*10^6;
              end
            end % for ilci
          end % if nS2mS
        end % for itype=
			end
				
      % Update output structure if necessary
      if merge_sets
        ALLEEG(ifile) = EEG; %%
        if save_set
          pop_saveset(EEG, 'filename', curr_fname, 'filepath', set_path, ...
            'check', 'on', ...
            'savemode', params.eeglab.savemode);
        end
      else
        out_st.data{outcols.subject_id}{end+1,1} = subid;
        out_st.data{outcols.filenum}(end+1,1) = ifile;
        
        if save_set
          pop_saveset(EEG, 'filename', curr_fname, 'filepath', set_path, ...
            'check', 'on', ...
            'savemode', params.eeglab.savemode);
        end
      end % if merge_sets
      
      if return_eeg
        out_st.data{outcols.EEG}{end+1,1} = EEG;
      else
        out_st.data{outcols.EEG}{end+1,1} = fullfile(set_path, curr_fname);
      end
    end % for ifile=
    
    if merge_sets && length(ALLEEG)>1
      EEG = pop_mergeset(ALLEEG,1:length(ALLEEG));
      pop_saveset(EEG, 'filename', curr_fname, 'filepath', set_path, ...
        'check', 'on', ...
        'savemode', params.eeglab.savemode);
      out_st.data{outcols.subject_id}{end+1,1} = subid;
      out_st.data{outcols.filenum}(end+1,1) = 1;
      if return_eeg
        out_st.data{outcols.EEG}{end+1,1} = EEG;
      else
        out_st.data{outcols.EEG}{end+1,1} = fullfile(set_path, curr_fname);
      end
    end
  end % for idir=

end % for isub

end % preproc_bdf

function [EEG, remove_chan_idxs] = remove_chans(EEG,chanlist,action)
  if nargin < 3
    action = 'remove';
  end
  
  switch action
    case 'remove'
      retain_chan_idxs = find(~ismember({EEG.chanlocs.labels}, chanlist));
    case 'keep'
      retain_chan_idxs = find(ismember({EEG.chanlocs.labels}, chanlist));
	end
  
	remove_chan_idxs = setdiff(1:size(EEG.data,1),retain_chan_idxs);
  if ~isempty(retain_chan_idxs)
    fprintf('Retaining %d channels ...\n', length(retain_chan_idxs));
    EEG.data(remove_chan_idxs,:) = [];
    EEG.chanlocs(remove_chan_idxs) = [];
    EEG.nbchan = size(EEG.data,1);
  end
end % remove_chans

function params = get_default_params
  params.bdf.file_prefix = '';
  params.srcfstub = '';
    params.sufix = 0;
  
  % Options for pop_biosig
  params.bdf.rmeventchan = 'on'; % should the event channel be removed
  params.bdf.blockrange = []; % range of data blocks to read in
  
  params.bdf.remove_chans = [];
  params.bdf.keep_chans = [];
  params.bdf.refchan = [];
  params.eeglab.merge_sets = false;  % merge EEG sets
  params.eeglab.use_session_dirs = false;
  
  params.eeglab.chanlocs.attach = false;
  params.eeglab.chanlocs.locfile = '';
  
  params.eeglab.reref = true;
  params.eeglab.remove_dc = true;
  
  % Parameters related to filtering. Different types of data, e.g. EEG, GSR,
  % Pleth, require different filtering parameters.  Thus, they are specified in
  % separate sub-structures of the filter structure.  If both high and low cutoff
  % are empty, no filtering is performed using pop_eegfilt
  
  params.eeglab.filter.eeg.chan_labels = biosemi_chan_labels('eeg');
  params.eeglab.filter.eeg.low_cutoff = [0.3];
  params.eeglab.filter.eeg.high_cutoff = [128];
  params.eeglab.filter.eeg.filt_order = [];
  params.eeglab.filter.eeg.filt_order_multiplier = 1;
  
  params.eeglab.filter.gsr.chan_labels = biosemi_chan_labels('gsr');
  params.eeglab.filter.gsr.low_cutoff = [];
  params.eeglab.filter.gsr.high_cutoff = [];
  params.eeglab.filter.gsr.filt_order = [];
  
  params.eeglab.filter.pleth.chan_labels = biosemi_chan_labels('pleth');
  params.eeglab.filter.pleth.low_cutoff = [1];
  params.eeglab.filter.pleth.high_cutoff = [];
  params.eeglab.filter.pleth.filt_order = [];
  
  params.eeglab.filter.resp.chan_labels = biosemi_chan_labels('resp');
  params.eeglab.filter.resp.low_cutoff = [];
  params.eeglab.filter.resp.high_cutoff = [5];
  params.eeglab.filter.resp.filt_order = [];

  params.eeglab.filter.other.chan_labels = {};
  params.eeglab.filter.other.low_cutoff = [];
  params.eeglab.filter.other.high_cutoff = [];
  params.eeglab.filter.other.filt_order = [];
  
  % Parameters related to saving various intermediates
  params.eeglab.save.set = true;
  params.eeglab.save.filtered = true;
  params.eeglab.save.clobber = 1;
  
  params.eeglab.savemode = 'twofiles';  % twofiles creates .set and .dat, whereas onefile creates .setp
  
  % Parameters influencing what data should be returned in the data struct. Note
  % that returning the entire EEG structure, as compared with only the filename,
  % is very memory intensive.
  params.eeglab.return.EEG = false;  % returns filename if false
  
end % get_default_params

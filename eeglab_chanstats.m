function out = eeglab_chanstats(data_st,params)

% Wrapper for using EEGLAB to get channel statistics for multiple subjects
% in parallel. Prints channel stats table to each subject's dataset directory. 
% This is helpful as a cross reference when deciding which channels should be rejected.
%
% Set for now to only calculate SD & kurtosis. The choice of statistics to
% calculate should eventually dynamically as a paramter.
%
% 02/07/13  Brian Hurley & Devin Platt  

eeglab('nogui')

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
  catch
    fprintf('Failed to complete channel statistic calculation: %s\n', subid);
  end

    numchans = size(EEG.data,1);

    SD = zeros(numchans,1);
    k = zeros(numchans,1);

    % Fire up parallel compute pool if possible
    if exist('matlabpool') && ~matlabpool('size') 
        matlabpool
    end
    
    % calculate desired channel stats
    parfor i=1 : numchans
        [~,SD(i),~,k(i),~,~,~,~,~,~,~]=signalstat(EEG.data(i,:), 0);
    end

    % print subject's channel stats to file
    table.fname = fullfile(set_path,sprintf('%s_chanstats.txt',subid));
    table.write2file = 1;   
    channum_vect = (1:numchans)';
    chanstat_data{1} = channum_vect;
    chanstat_data{2} = SD;
    chanstat_data{3} = k;
    chanstat = ensemble_init_data_struct;
    chanstat.vars = {'data','column_labels','column_formats'};
    chanstat.data{1} = chanstat_data;
    chanstat.data{2} = {'chan#','SD','kurtosis'};
    chanstat.data{3} = {'%d','%1.2f','%1.2f'};
    ensemble_display_table(chanstat,table);
  
end % for isub=

matlabpool close

out = [];
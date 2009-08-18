function result = eeglab_loadsets(indata,params)
% Loads one or more EEGLAB .set files based on condition information in
% params.condinfo and params.use_conditions
%
% 2009.08.15 PJ

result = ensemble_init_data_struct;
result.type = 'eeglab_setfile';

sessinfo = indata;
sesscols = set_var_col_const(sessinfo.vars);

%
% Perform any initial filtering
%
sessinfo = ensemble_filter(sessinfo,params.filt);
nsess = length(sessinfo.data{sesscols.session_id});

cond_names = params.use_conditions;
ncond = length(cond_names);

result.vars = {'EEG','condition','subject_id'};
resultcols = set_var_col_const(result.vars);
result.data = cell(size(result.vars));

ALLEEG = [];
for isess = 1:nsess
  subid = sessinfo.data{sesscols.subject_id}{isess};
  setpath = fullfile(params.path.project_root, subid, 'set');
  base_fstub = [subid params.fstub];
  
  for icond = 1:ncond
    curr_cond = cond_names{icond};
    setfile = sprintf('%s_%s.set',base_fstub, curr_cond);
    
    % load the dataset
    EEG = pop_loadset('filename', setfile, 'filepath', setpath);
    EEG.filename = setfile;
    EEG.filepath = setpath;
    EEG.subject = subid;
    EEG.condition = curr_cond;
    
    result.data{resultcols.EEG}{end+1,1} = EEG;
    result.data{resultcols.condition}{end+1,1} = curr_cond;
    result.data{resultcols.subject_id}{end+1,1} = subid;
  end % for icond=
end % for isess=


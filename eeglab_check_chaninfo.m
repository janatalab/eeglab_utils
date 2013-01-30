function EEG = eeglab_check_chaninfo(EEG,params)
% Checks to see whether channel information is attached to the EEG
% structure and if not, it trys to insert it from the file specified in
% params.preproc.eeglab.chanlocs.locfile
%
% EEG = eeglab_check_chaninfo(EEG.params);
%

% 29Jan2013 Petr Janata

try 
  locfile = params.preproc.eeglab.chanlocs.locfile;
catch
  locfile = '';
end

if isempty(locfile)
  error('No location file specified')
end

% Load the channel information from the location file. If not channel
% information is attached yet, run it so that pop_chanedit attaches the
% data, otherwise we are just checking for consistency
if ~isfield(EEG, 'chanlocs') || isempty(EEG.chanlocs)
  EEG = pop_chanedit(EEG,'load',params.preproc.eeglab.chanlocs.locfile);
else
  [chans,chaninfo] = pop_chanedit([],[],'load',params.preproc.eeglab.chanlocs.locfile);
  
  % Determine which of the channels are in the data
  [have_chan, srcIdxs] = ismember({EEG.chanlocs.labels}, {chans.labels});
  
  % Copy the data
  EEG.chanlocs = chans(srcIdxs);
  EEG.chaninfo = chaninfo;

  % Run the consistency check
  EEG = eeg_checkchanlocs(EEG);
end

return
end
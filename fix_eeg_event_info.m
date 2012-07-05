function EEG = fix_eeg_event_info(EEG)
% Corrects information in EEG event structure by re-reading the event info
% in the BDF file and correctly processing marker channels.
%
% EEG = fix_eeg_event_info(EEG);
%
% input EEG can either be the name of a .set file or an EEG struct

% Check to see if we are dealing with an EEG file or struct
if ischar(EEG)
	if exist(EEG,'file')
		load(EEG,'-mat');
	else
		error('Could not locate file: %s', EEG)
	end
end

% Get the name of the original BDF file
bdffname = EEG.etc.bdfhdr.FileName;

fprintf('Reading BDF info from %s\n', bdffname);

bdfEvents = read_bdf_marker_chans(bdffname);

% Match against urevents since this is what we need to fix first
oldEvents = EEG.urevent;
newEvents = oldEvents;
oldLatencies = [oldEvents.latency];

numNew = length(bdfEvents);
for iev = 1:numNew
	% Get the sample at which the event occurred 
	currLatency = bdfEvents(iev).latency;
	
	% See which of the old events match this one in latency
	latencyMask = oldLatencies == currLatency;
	
	% Handle case in which there is no matching event
	if ~any(latencyMask)
		fprintf('Found no match for event %d/%d\n', iev, numNew);
		continue 
	end
	
	% Handle case in which there are multiple events with the same latency.
	% This still needs to be fully implemented. For now, just throw an error
	if sum(latencyMask) > 1
		error('Found (%d) events with the same latency\n', sum(latencyMask));
		
		% Disambiguate based on event codes and their pattern of onsets
		
	else
		newEvents(latencyMask).type = bdfEvents(iev).type;
	end
end % for iev

% Copy modified event structure to urevent structure
EEG.urevent = newEvents;

% Now go through and update the event structure information based on the
% new urevent structure
numEvents = length(EEG.event);
for iev = 1:numEvents
	if ~isempty(EEG.event(iev).urevent)
		urtype = EEG.urevent(EEG.event(iev).urevent).type;
		if isnumeric(EEG.event(iev).type)
			EEG.event(iev).type = urtype;
		else
			EEG.event(iev).type = num2str(urtype);
		end
	end
end % for iev

% Save the modified .set file back out
EEG = eeg_checkset(EEG);
outfname = fullfile(EEG.filepath, EEG.filename);
fprintf('Saving modified set file: %s\n', outfname);
save(outfname, 'EEG');

return

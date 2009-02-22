function eeglab_eventmap(eventdata, params)
% eeglab_eventmap(eventdata, params)
%
% Takes EEGLAB event structure (EEG.event) and generates information about
% event distribution densities. Generate an image of when what events
% occur.

% 02/22/09 PJ

% Make sure we've set necessary params
try binsize_sec = params.binsize_sec; catch binsize_sec = 10; end

%
% Convert the structure into a cell array for easier manipulation
%
ev_vars = fieldnames(eventdata); % get variable names
ev_cols = set_var_col_const(ev_vars); % get column indices

eventdata = eventdata(:);  % make sure we are dealing with a column
eventdata = struct2cell(eventdata)';  % convert to cell array

% Pull out specific lists
eventlist = eventdata(:,ev_cols.type);

latencies = cell2mat(eventdata(:,ev_cols.latency));
latencies = latencies/params.srate;

% Determine the number of unique event types we have
unique_events = unique(eventlist);
nevent_types = length(unique_events);

fprintf('Found %d unique event types\n', nevent_types);

% Tabulate the frequencies
tabledata = tabulate(eventlist);

% Sort by count
[count,idx] = sort(cell2mat(tabledata(:,2)),'descend');
tabledata = tabledata(idx,:);

% Print the summary statistics
fprintf('Event\tCount\tPercent\tFirst\tLast\n')
for iev = 1:nevent_types
  event_mask = ismember(eventlist,tabledata(iev,1));
  fprintf('%s\t%d\t%2.2f\t%1.1f\t%1.1f\n', tabledata{iev,:},min(latencies(event_mask)), max(latencies(event_mask)));
end

fprintf('\n');
fprintf('Earliest event time (s): %1.3f \n', min(latencies));
fprintf('Latest event time (s): %1.3f \n', max(latencies));

% Create a histogram
timescale = fix(min(latencies)/binsize_sec)*binsize_sec:binsize_sec:ceil(max(latencies)/binsize_sec)*binsize_sec;
event_mtx = zeros(nevent_types,length(timescale));
for iev = 1:nevent_types
  event_mask = ismember(eventlist,unique_events{iev});
  event_mtx(iev,:) = hist(latencies(event_mask),timescale);
end

figure
imagesc(timescale,1:nevent_types,event_mtx), colorbar
set(gca,'ytick',1:nevent_types,'yticklabel', unique_events)
ylabel('Event')
xlabel('Time (s)')

try dataset = params.dataset; catch dataset = '(unspecified)'; end
title(sprintf('Map of event densities for dataset %s', strrep(dataset,'_','\_')))

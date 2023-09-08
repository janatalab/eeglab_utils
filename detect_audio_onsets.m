function onset_indices = detect_audio_onsets(data, params)
% Processes data and extracts onsets of audio events.
%
% onsets = detect_audio_onsets(data, params)
%
% data - the vector containing the audio channel data
% params - 
%   .Fs - the sampling rate of the data

if nargin < 2
  error('A params structure is required')
end

% Check our parameters
if ~isfield(params, 'Fs')
  error('sampling rate of data must be specified in params.Fs')
else
  Fs = params.Fs;
end

if ~isfield(params, 'min_ioi')
  error('Must specify minimum IOI, in sec, in params.min_ioi')
else
  min_ioi = params.min_ioi;
end

threshold = 50;
if isfield(params, 'threshold')
  threshold = params.threshold;
end

ramp_duration_sec = 4;
if isfield(params, 'ramp_duration_sec')
  ramp_duration_sec = params.ramp_duration_sec;
end

% Calculate the Nyquist frequency
nyquist = Fs/2;

% Apply a ramp at the beginning to eliminate high-frequency transient

nsamps_ramp = round(ramp_duration_sec*Fs);
ramp = (0:(nsamps_ramp-1))/nsamps_ramp;
data(1:nsamps_ramp) = data(1:nsamps_ramp) .* ramp;

% High-pass filter the data
Fc_high = 5;
[b,a] = butter(1, Fc_high/nyquist, 'high');
highpass = filter(b,a,data);


% Full-wave rectify the high-pass filtered signal
rectified = abs(highpass);

% Low-pass filter the signal with a zero-phase filter with a 2nd-order
% Butterworth filter

% Build the filter
Fc = 0.5; % default low-pass cutoff frequency
order = 1;
if isfield(params, 'butter')
  if isfield(params.butter, 'Fc')
    Fc = params.butter.Fc;
  end

  if isfield(params.butter, 'order')
    order = params.butter.order;
  end
end

[b,a] = butter(order, Fc/nyquist);

% Filter the data
filtered = filter(b,a,double(rectified));

% % Remove offset from filtered signal
% filtered_baseline_mean = mean(filtered(baseline_samps));
% filtered = filtered - filtered_baseline_mean;

% Identify periods of time in which the filtered signal exceeds a threshold

% Detect onsets:mat
above_threshold = filtered >= threshold;

% Identify onsets (transitions from zero to one)
onset_indices = find(diff(above_threshold)>0)+1;

onset_times_sec = onset_indices / Fs;

% Get the differences in onset times
iois = diff(onset_times_sec);

% Check whether all our inter-onset-intervals (iois) are valid
all_iois_valid = false;
last_good_ioi = 0;

while ~all_iois_valid
    curr_ioi = last_good_ioi + 1;

    if iois(curr_ioi) < min_ioi
        % Remove the onset
        onset_indices(curr_ioi+1) = [];
        onset_times_sec = onset_indices / Fs;

        % Recalculate the iois
        iois = diff(onset_times_sec);
    else
        last_good_ioi = curr_ioi;
    end

    if last_good_ioi == length(iois)
        all_iois_valid = true;
    end
end

return

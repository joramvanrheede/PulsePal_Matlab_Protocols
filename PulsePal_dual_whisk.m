%% PulsePal stimulus control script


%% Program details: adjust these to generate desired trial characteristics

% The following five properties allow for entering multiple values,
% generating more unique trials:
whisk_delays            = [.5]; % Whisk stimulus delays in s 
whisk_wave_freq         = [50]; % Whisker stimulus velocity (frequency of WAVEFORM in Hz)

%% IMPORTANT CHANGES 2018_01_31: 
multiple_whisks         = true; % if 'true' then it will trigger the whisker multiple times according to the frequency and duration specified
whisk_trig_freq         = [2 5 10 20 40]; % !!! Values must be LOWER THAN every 'whisk_wave_freq' value !!! Whisker stimulus triggering frequency (frequency of successive stimuli in HZ)
total_whisk_duration    = [4.0]; % whisker stimulation duration in seconds
%% END IMPORTANT CHANGES

led_delays              = [.25 5.25]; % LED stimulus delays in s
led_durations           = [4.5]; % LED stimulus durations

% These parameters currently only allow a single value per experiment:
n_repeats               = 10; % How many repeats of each parameter combination?

n_stimulators           = 2; % 1 or 2, depending on whether we are using one or two piezos

trial_length            = 10; % How long is each trial (for trial TTL)
trial_spacing           = 15; % Space between TRIGGERS for successive trials; NOTE - trial spacing incorporates time of trial execution; 


%% Advanced whisker stim parameters

v_max                   = 10; % Maximum voltage for whisking waveform (CAREFUL - it can't deal with negative voltages...)
stim_sample_duration    = [0.0002]; % 
sync_sample_duration    = 0.002;    % duration of each sample / pulse in the trial + whisk sync channel (total can't exceed 5000) 

%% Output channel assignment

whisk_wave_channel      = 1; % Which channel provides the whisking waveform for the amplifier
trial_whisk_ttl_channel = 2; % Which channel provides the TTL signal for when the whisker is on
led_ttl_channel         = 3; % Which channel provides the TTL trigger for the LED module
stim_switch_channel   	= 4; % Which channel provides the signal that determines which stimulator to use


%% Pulsepal startup and initialisation (reconnect to avoid trigger failures)
if exist('PulsePalSystem','var')
    EndPulsePal;
end

disp('Re-initialising PulsePal...')
pause(5)    % give this some time to complete in PulsePal system
PulsePal;   % Initialises PulsePal connection if none exists already; will give notification if PulsePal already connected

% Prerequisite: Pre-assigned PulsePalmatrix with default values
disp('Retrieving default stimulus properties')
load('DefaultPulsePalMatrix.mat');

disp('Sending to PulsePal...')
ProgramPulsePal(DefaultMatrix);
disp('complete.')
pause(0.2)

%% Generate a matrix of one 'block' containing all possible combinations of trials once
block_mat   = [];
counter     = 0;
for a = 1:length(whisk_delays) % for all whisk_delays...
    for b = 1:length(whisk_wave_freq) % for all whisk_wave_freq...
        for c = 1:length(whisk_trig_freq) % for all whisk_trig_freq...
            for d = 1:length(led_delays) % for all led_delays...
                for e = 1:length(led_durations) % for all led_durations...
                    for f = 1:n_stimulators % for each stimulator (currently 1 or 2)
                        counter = counter + 1; % increment trial counter to keep track of number of unique trials...
                        block_mat(counter,:)    = [whisk_delays(a) whisk_wave_freq(b) whisk_trig_freq(c) led_delays(d) led_durations(e) f]; % ...and generate a trial, appending it to the block.
                    end
                end
            end
        end
    end
end

%% Now generate a trial matrix block by block (using n_repeats * blocks), randomising the trial order of each block
trial_mat = [];
for a = 1:n_repeats
    trial_mat = [trial_mat; block_mat(randperm(counter),:)];
end

%% Program execution
n_stims             = size(trial_mat,1); % work out how many trials we are executing
tic     % start timer

for a = 1:n_stims
    disp(['Preparing trial ' num2str(a) '...']);
    
    this_whisk_delay        = trial_mat(a,1);                               % obtain whisking delay for this trial
    this_whisk_wave_freq    = trial_mat(a,2);                               % obtain whisking stimulus waveform frequency (i.e. velocity) for this trial
    this_whisk_trig_freq    = trial_mat(a,3);                               % obtain whisking trigger frequency (i.e. rate) for this trial
    this_led_delay          = trial_mat(a,4);                               % obtain LED onset delay for this trial
    this_led_duration       = trial_mat(a,5);                               % obtain LED duration for this trial
	this_whisk_stim         = trial_mat(a,6);                               % obtain stimulator ID for this trial
    
    % Display stimulus properties for this trial (also useful for debugging)
    disp([...
        'Whisk @ ' num2str(this_whisk_delay) 's; '...
        'Stim speed @ ' num2str(this_whisk_wave_freq) 'Hz; '...
        'Stim freq @ ' num2str(this_whisk_trig_freq) 'Hz; '...
        'LED @ ' num2str(this_led_delay) 's; '...
        'for @ ' num2str(this_led_duration) 's.'...
        'on stimulator ' num2str(this_whisk_stim)...
        ]);
    
    
    % whisker stimulus programming
    stim_duration           = 1 / this_whisk_wave_freq;                     % how long does this wave stimulus last?
    loop_duration           = 1 / this_whisk_trig_freq;
    
    if multiple_whisks
        n_whisk_stims   	= total_whisk_duration / loop_duration;
    else
        n_whisk_stims       = 1;
    end
    
    % how long should 1 stimulus loop (waveform + wait for next trigger) last?
    stim_burst_duration     = loop_duration * n_whisk_stims;                % how long should the entire stimulus burst last?
    n_samples               = stim_duration / stim_sample_duration;         % Calculate how many 1 ms samples the waveform needs to last 
    
    % stimulus waveform generation
    waveform_volts          = ((1 - cos([0:2*pi/n_samples:2*pi])) / 2);     % generates a [n_samples]-sample sinusoidal waveform, from 0 to 1 to 0.
    this_whisk_wave         = waveform_volts * v_max;                       % multiply waveform by target max voltage
    
    % append zero values to whisk_waveform to make it last until the next
    % whisker stimulation is due
    append_duration         = loop_duration - stim_duration;                % how much waiting time should be appended
    append_number           = round(append_duration / stim_sample_duration);% how many zero samples should be appended
    append_zeros            = zeros(1,append_number);                       % create a vector of the appropriate nr of zeros to append
    this_whisk_wave         = [this_whisk_wave append_zeros];               % append the zero values
    
    if length(this_whisk_wave) > 5000
        error(['Number of samples for whisking waveform exceeds PulsePal 2 maximum (5000) - offending number of samples: ' num2str(length(this_whisk_wave)) '; adjust trigger frequency or whisk sample rate'])
    end
    
	prepend_zeros           = zeros(1,this_whisk_delay / sync_sample_duration);                     % part of trial before whisker stim
    whisk_stim_on           = repmat(this_whisk_wave,1,n_whisk_stims) ~= 0; % representation of whenever whisk is happening based on whisker waveform above, but at wrong sample rate
    whisk_stim_inds         = 1:(sync_sample_duration/stim_sample_duration):length(whisk_stim_on);
    whisk_stim_on           = whisk_stim_on(whisk_stim_inds);    % resample to correspond to sync_sample_rate
    append_zeros            = zeros(1,(trial_length / sync_sample_duration) - (length(prepend_zeros) + length(whisk_stim_on))); % find out how many samples need to be appended to complete the trial
    
    sync_pulse_wave         = [prepend_zeros whisk_stim_on append_zeros];           % string all parts of the sync signal together
    sync_pulse_wave         = (sync_pulse_wave + 1) * 2.5;                          % set waveform to go from 2.5V (trial onset) to 5V (whisker stims)
    
    if length(sync_pulse_wave) > 5000
        error(['Number of samples for syncing waveform exceeds PulsePal 2 maximum (5000) - offending number of samples: ' num2str(length(this_whisk_wave)) '; adjust trial length or sync sample rate'])
    end
    
    % Now start populating the PulsePal parameter matrix
    stim_matrix             = DefaultMatrix;                                % copy default matrix for trial-specific editing
    
    % programming whisker waveform channel
    stim_matrix{15,whisk_wave_channel+1}        = 1;                            % 15: 'CustomTrainID'
    stim_matrix{17,whisk_wave_channel+1}        = 1;                            % 17: 'CustomTrainLoop'
    
    stim_matrix{12,whisk_wave_channel+1}        = this_whisk_delay;             % 12: 'PulseTrainDelay'
    stim_matrix{11,whisk_wave_channel+1}        = stim_burst_duration;          % 11: 'PulseTrainDuration'
    stim_matrix{5,whisk_wave_channel+1}         = stim_sample_duration;         % 5: 'Phase1Duration'
    
    % programming trial and whisk sync channel (trial = 2.5V, trial & whisk = 5V)
    stim_matrix{15,trial_whisk_ttl_channel+1}  	= 2;                            % 15: 'CustomTrainID'
    stim_matrix{17,trial_whisk_ttl_channel+1}  	= 1;                            % 17: 'CustomTrainLoop'
    
    stim_matrix{12,trial_whisk_ttl_channel+1} 	= 0;                            % 12: 'PulseTrainDelay' - Delay
    stim_matrix{11,trial_whisk_ttl_channel+1} 	= trial_length;                 % 11: 'PulseTrainDuration' - Duration of total nr of TTLs
    stim_matrix{5,trial_whisk_ttl_channel+1}  	= sync_sample_duration;      	% 5: 'Phase1Duration' - duration of each sample in pulse train (?)
    
    % programming LED TTL channel
    stim_matrix{12,led_ttl_channel+1}           = this_led_delay;               % 12: 'PulseTrainDelay' - Delay
    stim_matrix{5,led_ttl_channel+1}            = this_led_duration;            % 5: 'Phase1Duration' - Duration of TTL up
    stim_matrix{11,led_ttl_channel+1}           = this_led_duration;            % 11: 'PulseTrainDuration' - Duration of total nr of TTLs
    
    % programming stim_switch TTL channel
    stim_matrix{12,stim_switch_channel+1}       = 0;                            % 12: 'PulseTrainDelay' - Delay
    stim_matrix{5,stim_switch_channel+1}        = trial_length;                 % 5: 'Phase1Duration' - Duration of TTL up
    stim_matrix{11,stim_switch_channel+1}       = trial_length;                 % 11: 'PulseTrainDuration' - Duration of total nr of TTLs
    stim_matrix{3,stim_switch_channel+1}        = (this_whisk_stim - 1) * 5;    % 3: Phase1Voltage; 0V (low) for stimulator 1, 5V (high) for stimulator 2.

    % send the edited stimulus matrix to PulsePal
    ProgramPulsePal(stim_matrix);
    
    % Send the custom waveform to custom waveform slot nr. 1
    success1                                    = SendCustomWaveform(1, stim_sample_duration, this_whisk_wave); % send waveform to PulsePal; use 'stim_sample_duration' to set sample rate
    success2                                    = SendCustomWaveform(2, sync_sample_duration, sync_pulse_wave); % send trial synchronisation waveform to PulsePal; use 'stim_sample_duration' to set sample rate
    
    % All stimulus parameters are uploaded; start fast loop to monitor for 
    % elapsed time
    while toc < trial_spacing
    end
    
    % Now trigger all channels at the same time
    TriggerPulsePal('1111');
    tic % start timer until next trigger event
    
    % Report trial number
    disp(['Trial ' num2str(a) ' triggered']);
    pause(trial_length) % IMPORTANT! If you upload next trial parameters while the current trial is still running, it messes with the current trial
    
end

disp('Trial execution completed.')
EndPulsePal;

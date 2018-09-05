%% PulsePal stimulus control script


%% Program details: adjust these to generate desired trial characteristics

% The following five properties allow for entering multiple values,
% generating more unique trials:
whisk_delays            = [1]; % Whisk stimulus delays in s 
whisk_wave_freq         = [20]; % Whisker stimulus velocity (frequency of WAVEFORM in Hz)
whisk_sustain           = [.1]; % how long (in seconds) to sustain the whisker deflection. 

do_whisk                = true;
multiple_whisks         = true; % if 'true' then it will trigger the whisker multiple times according to the frequency and duration specified
whisk_trig_freq         = [4]; % !!! Values must be LOWER THAN every 'whisk_wave_freq' value !!! Whisker stimulus triggering frequency (frequency of successive stimuli in HZ)
total_whisk_duration    = [600]; % whisker stimulation duration in seconds
amplitudes              = [1]; % Amplitude as proportions of v_max_1 & v_max_2 below


%% do sth for this
do_led                  = true;
LED_delays              = [.5]; % LED stimulus delays in s
LED_durations           = [1]; % LED stimulus durations
LED_powers           	= [0]; % Power as PROPORTION of max (max = 1, --> 1000 mA)

LED_stim_type         	= 'step'; % supported: 'step', 'noise', to do: 'hwr noise' (half-wave rectified). To do: 'filtered noise' - Other things to consider: half-wave rectified noise, sine wave, ramped noise, sine wave modulated noise,

%% NOTE: if choosing noise opto stimulus, these need to be specified:
do_noise_filter         = true;     % bandpass filter noise?
noise_min_freq          = [1];      % minimum frequency of bandpass filter
noise_max_freq          = [100];    % maximum frequency of bandpass filter

LED_stim_envelope       = 'ramp up'; % supported: 'none', 'ramp up', 'ramp down', 'sinewave', 'triangle'


multiple_LED_stims     	= false;    % trigger multiple LED stims?
LED_trig_freqs        	= [.5];     % if multiple LED stims, what trigger frequency to use?
total_LED_duration      = [5];      % time window for triggering LED stimuli

stimulators             = 1;        % 1, 2 or [1 2] depending on whether we are using one or two piezos

% These parameters currently only allow a single value per experiment:

n_repeats               = 1;        % How many repeats of each parameter combination?


trial_length            = 8; % How long is each trial (for trial TTL)
trial_spacing           = 10; % Space between TRIGGERS for successive trials; NOTE - trial spacing incorporates time of trial execution; 

%%

debug                   = 1; % debug mode - run code but do not connect to PulsePal or send commands

%% Advanced whisker stim parameters

v_max_1                 = 10; % Maximum voltage for whisking waveform (CAREFUL - it can't deal with negative voltages...)
v_max_2                 = 5;

whisk_stim_sample_duration    = [0.0002]; % 
LED_sample_duration    = [0.001];    % duration of each sample / pulse in the trial + whisk sync channel (total can't exceed 5000) 

%% Output channel assignment

whisk_wave_channel      = 1; % Which channel provides the whisking waveform for the amplifier
trial_whisk_ttl_channel = 2; % Which channel 4provides the TTL signal for when the whisker is on
led_channel             = 3; % Which channel provides the TTL trigger for the LED module
stim_switch_channel   	= 4; % Which channel provides the signal that determines which stimulator to use

%% Checks for stimulus values that could damage the I/O systems:
if any(LED_powers) > 1 || any(LED_powers) < 0
    error('LED_powers values must be between 0 and 1')
end

if v_max_1 < 0 || v_max_2 < 0
    error('v_max values must be positive - negative voltages can damage the piezo system')
end

%% Pulsepal startup and initialisation (reconnect to avoid trigger failures)
if ~debug
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
end
%% Generate a matrix of one 'block' containing all possible combinations of trials once
block_mat   = [];
counter     = 0;
for a = 1:length(whisk_delays) % for all whisk_delays...
    for b = 1:length(whisk_wave_freq) % for all whisk_wave_freq...
        for c = 1:length(whisk_trig_freq) % for all whisk_trig_freq...
            for d = 1:length(LED_delays) % for all LED_delays...
                for e = 1:length(LED_durations) % for all LED_durations...
                    for f = 1:length(stimulators) % for each stimulator (currently 1 or 2)
                        for g = 1:length(amplitudes) % for all whisking amplitudes
                            for h = 1:length(LED_powers) % for all LED powers
                                for i = 1:length(LED_trig_freqs)
                                    counter = counter + 1; % increment trial counter to keep track of number of unique trials...
                                    block_mat(counter,:)    = [whisk_delays(a) whisk_wave_freq(b) whisk_trig_freq(c) LED_delays(d) LED_durations(e) stimulators(f) amplitudes(g) LED_powers(h) LED_trig_freqs(i)]; % ...and generate a trial, appending it to the block.
                                end
                            end
                        end
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

total_seconds       = n_stims * trial_spacing;

disp(' ')
disp(['Starting protocol totaling ' num2str(n_stims) ' trials...'])
disp(['Total execution time will be ' num2str(total_seconds / 60) ' minutes.'])

tic     % start timer

for a = 1:n_stims
    disp(' ')
    disp(['Trial ' num2str(a) ' of ' num2str(n_stims) '...']);
    this_whisk_delay        = trial_mat(a,1);                               % obtain whisking delay for this trial
    whisk_waveform_freq     = trial_mat(a,2);                               % obtain whisking stimulus waveform frequency (i.e. velocity) for this trial
    this_whisk_trig_freq    = trial_mat(a,3);                               % obtain whisking trigger frequency (i.e. rate) for this trial
    this_LED_delay          = trial_mat(a,4);                               % obtain LED onset delay for this trial
    this_LED_duration       = trial_mat(a,5);                               % obtain LED duration for this trial
	this_whisk_stim         = trial_mat(a,6);                               % obtain stimulator ID for this trial
    this_amplitude          = trial_mat(a,7);                               % obtain whisk stim amplitude for this trial
    this_LED_power          = trial_mat(a,8);                               % obtain LED power for this trial
    this_LED_trig_freq      = trial_mat(a,9);                               % obtain LED trigger freq for this trial
    
    % Display stimulus properties for this trial (also useful for debugging)

    disp([...
        'Whisk @ ' num2str(this_whisk_delay) 's; '...
        'Speed ' num2str(whisk_waveform_freq) 'Hz; '...
        'Freq ' num2str(this_whisk_trig_freq) 'Hz; '...
        'LED @ ' num2str(this_LED_delay) 's '...
        'for ' num2str(this_LED_duration) 's; '...
        'stim ' num2str(this_whisk_stim) '; ' ...
        'amplitude ' num2str(this_amplitude) '; ' ...
        'LED ' num2str(this_LED_power*1000) 'mA.' ...
        ]);
    
    % whisker stimulus programming
    whisk_stim_duration         = 1 / whisk_waveform_freq + whisk_sustain;                     % how long does this wave stimulus last?
    whisk_loop_duration      	= 1 / this_whisk_trig_freq;
    
    if whisk_stim_duration > whisk_loop_duration
        error('Whisker stimulus duration is too long for this stimulus trigger frequency')
    end
    
    if multiple_whisks
        n_whisk_stims           = total_whisk_duration / whisk_loop_duration;
    else
        n_whisk_stims           = 1;
    end
    
    % how long should 1 stimulus loop (waveform + wait for next trigger) last?
    whisk_stim_burst_duration 	= whisk_loop_duration * n_whisk_stims;                % how long should the entire stimulus burst last?
    n_samples                   = (1 / whisk_waveform_freq) / whisk_stim_sample_duration;         % Calculate how many 1 ms samples the waveform needs to last 
    
    % stimulus waveform generation
    whisk_waveform_volts    	= [((1 - cos([0:2*pi/n_samples:pi])) / 2) ones(1,whisk_sustain/whisk_stim_sample_duration) ((1 - cos([pi:2*pi/n_samples:2*pi])) / 2)] ;     % generates a [n_samples]-sample sinusoidal waveform, from 0 to 1 to 0.
    
    if this_whisk_stim == 1
        whisk_waveform          = whisk_waveform_volts * v_max_1;                       % multiply waveform by target max voltage
    elseif this_whisk_stim == 2
        whisk_waveform          = whisk_waveform_volts * v_max_2; 
    end
    
    whisk_waveform              = whisk_waveform * this_amplitude;
    
    % append zero values to whisk_waveform to make it last until the next
    % whisker stimulation is due
    whisk_append_duration   	= whisk_loop_duration - whisk_stim_duration;                % how much waiting time should be appended
    whisk_append_number      	= round(whisk_append_duration / whisk_stim_sample_duration);% how many zero samples should be appended
    whisk_append_zeros      	= zeros(1,whisk_append_number);                       % create a vector of the appropriate nr of zeros to append
    whisk_waveform            	= [whisk_waveform whisk_append_zeros];               % append the zero values
    
    if length(whisk_waveform) > 5000
        error(['Number of samples for whisking waveform exceeds PulsePal 2 maximum (5000) - offending number of samples: ' num2str(length(whisk_waveform)) '; adjust trigger frequency or whisk sample rate'])
    end
    

    %% LED stimulus wave generation
    
    LED_stim_duration       = this_LED_duration;
    LED_loop_duration       = 1 / this_LED_trig_freq;
    
    if multiple_LED_stims
        n_LED_stims         = total_LED_duration / LED_loop_duration;
    else
        n_LED_stims     	= 1;
    end
    
    LED_stim_burst_duration = LED_loop_duration * n_LED_stims;
    n_LED_samples        	= this_LED_duration / LED_sample_duration;
    
    switch LED_stim_type
        case 'step' 
            LED_waveform    = ones(1,n_LED_samples)  * this_LED_power;

        case 'noise'
            LED_waveform    = randn(1,n_LED_samples);                 	% sample from gaussian distribution
            
            if do_noise_filter
                [filt_b,filt_a] = butter(2, [noise_min_freq noise_max_freq]/((1/LED_sample_duration)/2));
                LED_waveform    = filter(filt_b, filt_a, LED_waveform);
            end
            
            LED_waveform    = LED_waveform - min(LED_waveform);     % min --> 0
            LED_waveform    = LED_waveform / max(LED_waveform);     % max --> 1
    end
    
    switch LED_stim_envelope
        case 'ramp up'
            this_LED_envelope   = ([1:n_LED_samples] / n_LED_samples);
        case 'ramp down'
            this_LED_envelope   = (1 - [1:n_LED_samples] / n_LED_samples);
        case 'triangle'
            this_LED_envelope   = 1 - abs(1-([1:n_LED_samples] / n_LED_samples)*2);
        case 'sinewave'
            this_LED_envelope   = cos(linspace(1*pi,3*pi,n_LED_samples))/2 + 0.5;
        case 'none'
            this_LED_envelope   = ones(1,n_LED_samples);
        otherwise
            error(['Unsupported LED_stim_envelope: ' LED_stim_envelope])
    end
    
    LED_waveform    = LED_waveform .* this_LED_envelope;
    
    
    LED_append_duration   	= LED_loop_duration - this_LED_duration;
    LED_append_number       = round(LED_append_duration / LED_sample_duration);
    LED_append_zeros        = zeros(1,LED_append_number);
    LED_waveform            = [LED_waveform LED_append_zeros];
    
    if length(LED_waveform) > 5000
        error(['Number of samples for syncing waveform exceeds PulsePal 2 maximum (5000) - offending number of samples: ' num2str(length(whisk_waveform)) '; adjust trial length or sync sample rate'])
    end
    
    if debug
        figure
        plot(whisk_waveform)
        figure
        plot(LED_waveform)
        
        keyboard
        continue
    end
    
    % Now start populating the PulsePal parameter matrix
    stim_matrix                                 = DefaultMatrix;            	% copy default matrix for trial-specific editing
    
    % programming whisker waveform channel
    stim_matrix{15,whisk_wave_channel+1}        = 1;                            % 15: 'CustomTrainID'
    stim_matrix{17,whisk_wave_channel+1}        = 1;                            % 17: 'CustomTrainLoop'
    
    stim_matrix{12,whisk_wave_channel+1}        = this_whisk_delay;             % 12: 'PulseTrainDelay'
    stim_matrix{11,whisk_wave_channel+1}        = whisk_stim_burst_duration;          % 11: 'PulseTrainDuration'
    stim_matrix{5,whisk_wave_channel+1}         = whisk_stim_sample_duration;         % 5: 'Phase1Duration'
    
    % programming whisk TTL channel
    stim_matrix{12,trial_whisk_ttl_channel+1} 	= this_whisk_delay;            	% 12: 'PulseTrainDelay' - Delay
    stim_matrix{5,trial_whisk_ttl_channel+1}   	= whisk_stim_duration;                % 5: 'Phase1Duration' - Duration of TTL up
    stim_matrix{8,trial_whisk_ttl_channel+1}  	= whisk_loop_duration - whisk_stim_duration; % 8: 'InterPhaseInterval' - Duration of TTL up
    
    stim_matrix{11,trial_whisk_ttl_channel+1} 	= whisk_stim_burst_duration;       	% 11: 'PulseTrainDuration' - Duration of total nr of TTLs
    stim_matrix{3,trial_whisk_ttl_channel+1}   	= 5;                            % 3: Phase1Voltage; 5V for 'on'
    % stim_matrix{18,trial_whisk_ttl_channel+1}       = 2.5; ?
    
    % programming LED channel
    stim_matrix{15,led_channel+1}               = 2;                            % 15: 'CustomTrainID'
    stim_matrix{17,led_channel+1}               = 1;                            % 17: 'CustomTrainLoop'
    
    stim_matrix{12,led_channel+1}               = this_LED_delay;               % 12: 'PulseTrainDelay' - Delay
    stim_matrix{11,led_channel+1}               = LED_stim_burst_duration;   	% 11: 'PulseTrainDuration' - Duration of total nr of TTLs
    stim_matrix{5,led_channel+1}                = LED_sample_duration;      	% 5: 'Phase1Duration' - duration of each sample in pulse train (?)
    
    % programming stim_switch TTL channel
    stim_matrix{12,stim_switch_channel+1}       = 0;                            % 12: 'PulseTrainDelay' - Delay
    stim_matrix{5,stim_switch_channel+1}        = trial_length;                 % 5: 'Phase1Duration' - Duration of TTL up
    stim_matrix{11,stim_switch_channel+1}       = trial_length;                 % 11: 'PulseTrainDuration' - Duration of total nr of TTLs
    stim_matrix{3,stim_switch_channel+1}        = (this_whisk_stim - 1) * 5;    % 3: Phase1Voltage; 0V (low) for stimulator 1, 5V (high) for stimulator 2.
    
    % send the edited stimulus matrix to PulsePal
    ProgramPulsePal(stim_matrix);
    
    % Send the custom waveform to custom waveform slot nr. 1
    success1                                    = SendCustomWaveform(1, whisk_stim_sample_duration, whisk_waveform); % send waveform to PulsePal; use 'whisk_stim_sample_duration' to set sample rate
    success2                                    = SendCustomWaveform(2, LED_sample_duration, LED_waveform); % send trial synchronisation waveform to PulsePal; use 'whisk_stim_sample_duration' to set sample rate
    
    % All stimulus parameters are uploaded; start fast loop to monitor for 
    % elapsed time
    while toc < trial_spacing
    end
    
    % Now trigger all channels at the same time
    TriggerPulsePal('1111');
    tic % start timer until next trigger event
    
    % Report trial number
    disp(['Trial ' num2str(a) ' triggered']);
    disp(['Time remaining: ~ ' num2str((n_stims-a)*trial_spacing/60) ' minutes.'])
    disp(' ')
    
    pause(trial_length) % IMPORTANT! If you upload next trial parameters while the current trial is still running, it messes with the current trial
    
end

disp('Trial execution completed.')
EndPulsePal;

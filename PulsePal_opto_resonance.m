%% PulsePal stimulus control script

LED_trig_freq           = [1 2 4 8 16 32 64 100]; % !!! Values must be LOWER THAN every 'whisk_wave_freq' value !!! Whisker stimulus triggering frequency (frequency of successive stimuli in HZ)
LED_delay               = [.5]; % LED stimulus delays in s
total_LED_duration      = [6]; % LED stimulation duration in seconds

% These parameters currently only allow a single value per experiment:
n_repeats               = 10; % How many repeats of each parameter combination?

trial_length            = 12; % How long is each trial (for trial TTL)
trial_spacing           = 20; % Space between TRIGGERS for successive trials; NOTE - trial spacing incorporates time of trial execution; 

%% Output channel assignment

trial_ttl_channel       = 2; % Which channel provides the TTL signal for when the whisker is on
led_ttl_channel     	= 3; % Which channel provides the TTL trigger for the LED module

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

%% Now generate a trial matrix block by block (using n_repeats * blocks), randomising the trial order of each block
trial_mat = [];
for a = 1:n_repeats
    trial_mat = [trial_mat; LED_trig_freq(randperm(counter),:)];
end

%% Program execution
n_stims             = size(trial_mat,1); % work out how many trials we are executing
tic     % start timer

for a = 1:n_stims
    disp(['Preparing trial ' num2str(a) '...']);
    
    this_LED_trig_freq    = trial_mat(a);                               % obtain whisking trigger frequency (i.e. rate) for this trial
    
    % Display stimulus properties for this trial (also useful for debugging)
    disp([...
        'Opto stim freq @ ' num2str(this_LED_trig_freq) 'Hz; '...
        'LED @ ' num2str(LED_delay) 's; '...
        'on stimulator ' num2str(this_whisk_stim)...
        ]);
    
    % whisker stimulus programming
    stim_duration           = .5 / this_LED_trig_freq;                    % how long does this wave stimulus last?
    loop_duration           = 1 / this_LED_trig_freq;
    
    n_opto_stims            = total_whisk_duration / loop_duration;
    total_LED_duration
    
    
    % Now start populating the PulsePal parameter matrix
    stim_matrix             = DefaultMatrix;                                % copy default matrix for trial-specific editing
    
    % programming trial and whisk sync channel (trial = 2.5V, trial & whisk = 5V)

    stim_matrix{12,trial_whisk_ttl_channel+1} 	= 0;                            % 12: 'PulseTrainDelay' - Delay
    stim_matrix{11,trial_whisk_ttl_channel+1} 	= trial_length;                 % 11: 'PulseTrainDuration' - Duration of total nr of TTLs
    stim_matrix{5,trial_whisk_ttl_channel+1}  	= trial_length;                 % 5: 'Phase1Duration' - duration of phase 1 of pulse

    % programming LED TTL channel
    stim_matrix{12,led_ttl_channel+1}           = LED_delay;                    % 12: 'PulseTrainDelay' - Delay
    stim_matrix{5,led_ttl_channel+1}            = stim_duration;                % 5: 'Phase1Duration' - Duration of TTL up
    stim_matrix(8,led_ttl_channel+1)            = loop_duration-stim_duration;  % 8: InterPulseInterval
    stim_matrix{11,led_ttl_channel+1}           = total_LED_duration;           % 11: 'PulseTrainDuration' - Duration of total nr of TTLs
    
    % send the edited stimulus matrix to PulsePal
    ProgramPulsePal(stim_matrix);
    
    % All stimulus parameters are uploaded; start fast loop to monitor for 
    % elapsed time
    while toc < trial_spacing
    end
    
    % Now trigger all channels at the same time
    TriggerPulsePal('0110'); % trigger channels 2 and 3 of pulsepal
    tic % start timer until next trigger event
    
    % Report trial number
    disp(['Trial ' num2str(a) ' triggered']);
    pause(trial_length) % IMPORTANT! If you upload next trial parameters while the current trial is still running, it messes with the current trial
    
end

disp('Trial execution completed.')
EndPulsePal;


function [spindles] = Spindle_Detect(num_CH, spi_CH, varargin )
% SPINDLE_DETECT Detects Spindle events in LFP file
%       SPINDLES = SPINDLE_DETECT(NUM_CH, SPI_CH) uses default parameters
%       SPINDLES = SPINDLE_DETECT(NUM_CH, SPI_CH, passband, state, med_thresh,
%       ied_int, minburstInt, mean_thresh) uses specified parameters 
%
% REQUIREMENTS
%       NUM_CH: total number of channels
%       SPI_CH: channel for spindle analysis
%       Note:   must be in path of .lfp file and -states.mat(from StateEditor)
%                   file
%
% OPTIONAL INPUT ARGUMENTS
%       SPI_THRESHOLDS:     specifies thresholds for spindle beginning/end 
%                           and peak, in multiples of standard deviation 
%       STATE:              select between WAKE, NREM, or REM
%
% OUTPUT
%       SPINDLE.mat
%           named with fbasename+IED_CH where column 1 = start; 
%           column 2 = peak; column 3 = end; column 4 = power;
%       .spi.evt file
%           named with fbasename+IED_CH
%
% EXAMPLES:
% [SPINDLES] = SPINDLE_DETECT(128,94)
% [SPINDLES] = SPINDLE_DETECT(128,94, [2 3], [7 12], 'NREM');
%
% Naureen Ghani and Jennifer Gelinas (2018)

%Parse Inputs
p = inputParser;
defaultSpi_thresholds = [2 4];
defaultState = 'NREM';

addRequired(p, 'num_CH', @isnumeric);
addRequired(p, 'IED_CH', @isnumeric);
addOptional(p, 'spi_thresholds', defaultSpi_thresholds, @isivector);
addOptional(p, 'state', defaultState, @is_string);

parse(p, num_CH, spi_CH , varargin{:});

% Obtain States
state_mat = dir('*-states*');
load (state_mat.name);
StateIntervals = ConvertStatesVectorToIntervalSets(states);                 % 6 Intervalsets representing sleep states
REM = StateIntervals{5};
NREM = or(StateIntervals{2}, StateIntervals{3});
WAKE = StateIntervals{1};

% Define state
state = p.Results.state;

% State parameter
if strcmp(state,'NREM')
    state = NREM;
elseif strcmp(state, 'REM')
    state = REM;
else
    strcmp(state,'WAKE')
    state = WAKE;
end

% Define thresholds
spi_lowThresholdFactor = p.Results.spi_thresholds(1);
spi_highThresholdFactor = p.Results.spi_thresholds(2);

% Load file and restrict to state
filename = dir('*.lfp');
[~, fbasename, ~] = fileparts(filename.name);
lfp = LoadLfp(fbasename,num_CH,spi_CH); 
restricted_lfp = [Range(Restrict(lfp,state),'s') Data(Restrict(lfp,state))];

% Downsample 
dLfp = restricted_lfp(:,2);
dLfp = resample(dLfp, 125, 1250); 
rg = restricted_lfp(:,1);
rg = rg(1:10:end);

% Bandpass Filter
Wn = [9 16]; % spindle freq range [Hz]
rs = 125; % new sampling rate
[b,a] = butter(9, 2*Wn/rs, 'bandpass');
fil_sleep = filtfilt(b,a,dLfp);
fil_sleep = [rg fil_sleep];

% Spindle Detection
[~, sd, ~] = FindRipples(fil_sleep,  'thresholds', [spi_lowThresholdFactor spi_highThresholdFactor], 'durations', [400 3000 300], 'frequency', 125); 
[spindles, ~, ~] = FindRipples(fil_sleep, 'thresholds', [spi_lowThresholdFactor spi_highThresholdFactor], 'durations', [400 3000 300], 'frequency', 125, 'stdev', sd); 

% Remove Cross-Detected IEDs
% IED_mat = dir('*_IED.mat*');
% load (IED_mat.name);
% IED_times = IED(spi_CH,1);

IED_times = [];

if isempty(IED_times)==1
    disp('No IED'); 
    % Save spindle file as is
    spindle_file = strcat(fbasename, '_', num2str(spi_CH), 'spindles'); 
    save(spindle_file, 'spindles');
    spindle_events = strcat(fbasename, '_', num2str(spi_CH), '.spi.evt');
    channelID = spi_CH-1; 
    SaveRippleEvents(spindle_events, spindles, channelID);
elseif isempty(spindles)==1
    disp('No Spindles'); 
    % Save spindle file as is
    spindle_file = strcat(fbasename, '_', num2str(spi_CH), 'spindles'); 
    save(spindle_file, 'spindles');
    spindle_events = strcat(fbasename, '_', num2str(spi_CH), '.spi.evt');
    channelID = spi_CH-1; 
    SaveRippleEvents(spindle_events, spindles, channelID);   
else
    spi = spindles(:,1); 
    diff_val = zeros(1); 
    for i = 1:length(spi)
        ind = knnsearch(IED_times,spi(i,1));
        diff_val(i) = spi(i,1)-IED_times(ind); 
    end
    
    [~,ind2] = find(diff_val<=0.05); % start time within 50 ms
    spindles(ind2,:) = [];

    % Save revised spindle file
    spindle_file = strcat(fbasename, '_', num2str(spi_CH), 'spindles'); 
    save(spindle_file, 'spindles');
    spindle_events = strcat(fbasename,'_', num2str(spi_CH), '.spi.evt');
    channelID = spi_CH-1; 
    SaveRippleEvents(spindle_events, spindles, channelID);
end


end
















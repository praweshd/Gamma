function [gamma] = Gamma_Detect(num_CH, gam_CH, varargin )
% GAMMA_DETECT Detects Gamma events in LFP file
% Prawesh Dahal (last verified Sept 3, 2019)

%       gamma = Gamma_Detect(NUM_CH, SPI_CH) uses default parameters
%       gamma = Gamma_Detect(NUM_CH, SPI_CH, passband, state, med_thresh,
%       ied_int, minburstInt, mean_thresh) uses specified parameters 
%
% REQUIREMENTS
%       NUM_CH: total number of channels
%       GAM_CH: channel for spindle analysis
%       Note:   must be in path of .lfp file and -states.mat(from StateEditor)
%                   file
%
% OPTIONAL INPUT ARGUMENTS
%       GAM_THRESHOLDS:     specifies thresholds for spindle beginning/end 
%                           and peak, in multiples of standard deviation 
%       STATE:              select between WAKE, NREM, or REM
%
% OUTPUT
%       GAMMA.mat
%           named with fbasename+IED_CH where column 1 = start; 
%           column 2 = peak; column 3 = end; column 4 = power;
%       gam.evt file
%           named with fbasename+IED_CH
%
% EXAMPLES:
% [GAMMA] = GAMMA_DETECT(205,10)
% [GAMMA] = GAMMA_DETECT(205,10, [2 3], [7 12], 'NREM');
%
% Adapted from Naureen Ghani and Jennifer Gelinas (2018)
% Prawesh Dahal (July 3, 2019) 

%Parse Inputs
p = inputParser;
defaultGam_thresholds = [1 3]; %low [0.5 3], high [1 5]
defaultState = 'NREM';

addRequired(p, 'num_CH', @isnumeric);
addRequired(p, 'IED_CH', @isnumeric);
addOptional(p, 'spi_thresholds', defaultGam_thresholds, @isivector);
addOptional(p, 'state', defaultState, @is_string);

parse(p, num_CH, gam_CH , varargin{:});

% Obtain States
state_mat = dir('*-states*');
load (state_mat.name);
StateIntervals = ConvertStatesVectorToIntervalSets(states);                 % 6 Intervalsets representing sleep states
REM = StateIntervals{5};
%NREM = or(StateIntervals{2}, StateIntervals{3});
NREM = StateIntervals{3};
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
 
% disp(['Loading LFP for CH ', num2str(gam_CH)])
% Load file and restrict to state
filename = dir('*.lfp');
[~, fbasename, ~] = fileparts(filename.name);
lfp = LoadLfp(fbasename,num_CH,gam_CH); 
restricted_lfp = [Range(Restrict(lfp,state),'s') Data(Restrict(lfp,state))];

interg = 50; %Inter-Gamma Duration - low 100, high 50
maxg = 200; %Max Gamma Duration     - low 300, high 350
ming = 50; %Min Gamma Duration      - low 50, high 20

disp(['Gamma Detection for CH ', num2str(gam_CH)])

     
[gamma] = FindGamma_I(restricted_lfp,'short','thresholds', [spi_lowThresholdFactor spi_highThresholdFactor], 'durations', [interg maxg ming], 'frequency', 1250); 
    
    
    
    
  
    
     
 


end



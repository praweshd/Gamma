function lfp = LoadGamma(fbasename,ext,nchannels,channels,varargin)

% Loads LFP in a tsd object. Default sampling frequency @ 1250 Hz (or state otherwise)
%
% USAGE:
% 
%     lfp = LoadGamma(fbasename,ext,nchannels,channels)
%
% INPUT:
%     fbasename: session file basename
%     ext: extension can either be '.gam' or '.lfp'
%     nchannels: total number of channels in the 
%     channels: vector of channel index
%
% OPTION
% lfp = LoadLfp(fbasename,ext, nchannels,channels,'optionname',optionvalue)
% same options as LoadBinary

% Prawesh Dahal June 19, 2019
% Adapted from Adrien Peyrache 2012

duration = Inf;
start = 0;
precision = 'int16';
duration = Inf;
frequency = 1250;

for i = 1:2:length(varargin),
  if ~isa(varargin{i},'char'),
    error(['Parameter ' num2str(i+3) ' is not a property (type ''help LoadBinary'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'duration',
      duration = varargin{i+1};
      if ~isa(duration,'numeric') | length(duration) ~= 1 | duration < 0,
        error('Incorrect value for property ''duration'' (type ''help LoadBinary'' for details).');
      end
    case 'frequency',
      frequency = varargin{i+1};
      if ~isa(frequency,'numeric') | length(frequency) ~= 1 | frequency <= 0,
        error('Incorrect value for property ''frequency'' (type ''help LoadBinary'' for details).');
      end
    case 'start',
      start = varargin{i+1};
      if ~isa(start,'numeric') | length(start) ~= 1,
        error('Incorrect value for property ''start'' (type ''help LoadBinary'' for details).');
      end
		if start < 0, start = 0; end
      case 'precision',
      precision = varargin{i+1};
      if ~isa(precision,'char'),
        error('Incorrect value for property ''precision'' (type ''help LoadBinary'' for details).');
      end
  end
end
        
%  if frequency == 125
%      ext = '.spi';
%  else
%      ext = '.lfp';
%  end

lfp = LoadBinary([fbasename ext],'nchannels',nchannels,'channels',channels,'duration',duration,'precision',precision,'start',start,'frequency',frequency);
t = 10000*[0:length(lfp)-1]/frequency;

lfp = tsd(t',lfp);



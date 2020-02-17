%% Create bad channels mat file
% Prawesh Dahal 
% 09/11/19
% bad_ch in neuroscope
function createbadchannels(bad_ch)

lfp_file = dir('*.lfp');
[~, fbasename, ~] = fileparts(lfp_file.name);
bad_ch = sort(bad_ch); 
bad_ch = bad_ch + 1;
savename = [fbasename, '__bad_CH.mat']; 
save(savename,'bad_ch')

end


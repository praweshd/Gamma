%% dat_filhil

% data_filhil filters and Hilbert transforms the data from DatTracker
% Prawesh Dahal (June 25, 2019) 

function [data_fil, data_hil] = data_filhil(data,Rs,Wn,trial_N,duration,CH_N,resample_factor)
  
 % Gamma freq range [Hz]

%Filter
n=3;

Rs_resample=Rs/resample_factor;
[b,a]=butter(n,2*Wn/Rs_resample,'bandpass');

if length(size(data)) > 2
    tmp=reshape (data,CH_N,trial_N*duration);
    data = resample(tmp',1,resample_factor)';
end 


%Hilbert
Data_fil_tmp=filtfilt (b,a,data')';
Data_hil_tmp=hilbert (Data_fil_tmp')';

if length(size(data)) > 2
    data_fil=reshape(Data_fil_tmp,CH_N,duration/resample_factor,trial_N);
    data_hil=reshape(Data_hil_tmp,CH_N,duration/resample_factor,trial_N);  
end 
  
data_fil= Data_fil_tmp;
data_hil= abs(Data_hil_tmp);

    
end
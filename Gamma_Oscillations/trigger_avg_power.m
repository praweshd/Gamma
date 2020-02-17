%% trigger_average_power 

% data_filhil filters and Hilbert transforms the data from DatTracker
% Prawesh Dahal (June 25, 2019) 

function [gam_pwr, gam_hil, data_hil_avg] = trigger_avg_power(data_hil, duration, trigger_window) 

Rs = 1250; 
trigger_window = (trigger_window/1000)*Rs;

start_point=round((duration/2) - trigger_window );
end_point= round(start_point   + 2*trigger_window);

data_hil_avg=mean(abs(data_hil),3);
gam_hil=data_hil_avg(:,start_point:end_point);
%gam_pwr=sum(gam_hil,2)./(length (start_point:end_point)); 

% gam_pwr=sum(gam_hil,2)./median(gam_hil,2);

gam_pwr = max(gam_hil,[],2); 
 
    
end
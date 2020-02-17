%% Create LFP files

close all
clear all
clc

% PrePro_Dat_merger_automatic

%day_n=1;
days=dir('2018*');

 %%

 %
for i = 11:20

   day=days(i).name;
   cd(day);
   
   disp(['Running Preprocessing on:', [pwd]])
   PrePro_Dat_merger_automatic
   
   cd .. 
    
   
end 
 
  
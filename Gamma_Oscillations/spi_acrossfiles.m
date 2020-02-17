%% SPI/RP/GAM across files
% This program goes across all the baseline files with a mainfolder (e.g. ORat_14) and calculates 
% desired detection. 

close all
clear all
clc

%Make sure you are inside a folder of one animal (e.g. inside ORat_14)
mainfolder = pwd; 

%Find all the baseline files 
files = dir('*_baseline');

for i = 1 : length(files)

baselinef = strcat(pwd,'/',files(i).name); 

cd(baselinef)   %go into baseline 

%Detection to run

lfp2spi_detect_wavelet

cd(mainfolder) % go back to mainfolder

end 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
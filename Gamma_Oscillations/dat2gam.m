%% dat2gam - LFP to GAM
%Uses Wavelet-based thresholding to detect power increase only in the gamma
%band
%Adapted from dat2spi.m
%Prawesh Dahal - June 10, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%=========================================================================
%========================================================================= 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%Data_dat is the raw_lfp (yval X channel)

function [S] = dat2gam(Data_dat,resample_factor,arburg_n,Rs,freq,smooth_win,gamb,lowb,highb)

%Downsample
Data_gam=   resample(Data_dat,1,resample_factor);
[c,~]= size(Data_gam);
k=0;

for i = 10 %: c
    
%     disp(['Doing for channel ', num2str(i)])
    
    x=Data_gam(i,:);
    
    %Whiten
    b = arburg(x,arburg_n);
    y = Filter0In(b, x);
    
    %Wavelet
    [S_tmp,~,~] = awt_freqlist(y,Rs/resample_factor,freq,'Gabor');
    
    %Smoothen
    S_tmp_smooth = smoothdata(abs(S_tmp),1,'gaussian',smooth_win);
    
    %Ratio
    gam_band = sum(S_tmp_smooth(:, gamb(1):gamb(2) ), 2);
    
    low_band = sum(S_tmp_smooth(:, lowb(1):lowb(2) ), 2);
    
    high_band = sum(S_tmp_smooth(:, highb(1):highb(2) ), 2);
    
    gam_ratio= (gam_band - low_band - high_band) ./ (gam_band + low_band + high_band);

    k = k+1; 
    S(k,:)  = gam_ratio;
    
    
end
 
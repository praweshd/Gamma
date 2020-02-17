%% LFP to Wavelet-Based Gamma Ratio 
% Prawesh Dahal
% June 18, 2019

%This streams with 60s duration over all the LFP to create the 
%wavelet-based increase in gamma band ratio file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=============================================================
%=============================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lfp2gam_detect_wavelet(freq,gamb,lowb,highb)
%%
tic
%Load file
lfp_file=dir('*.lfp');
lfp_filename=lfp_file.name;
filename=lfp_filename;
basename= lfp_filename(1:end-4);

Rs=1250;
CH_N=xml2CH_N(cat(2,basename,'.xml'));
 

resample_factor= 1;
f_out_spi= cat(2,basename,'__SPI_wavelet_ratio.spi');
fid_lfp=fopen(f_out_spi,'w+');

duration=120*Rs;
s_time=60*Rs;
f_duration_sec=dat_length(lfp_filename,CH_N,Rs);

arburg_n=1;
smooth_win= round(Rs*(0.01)/resample_factor);

disp(cat(2,'The duration of LFP is ',num2str(f_duration_sec),'s... working...'));
pad_length = 12 *Rs;

%% ============================================================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% If the file is less than 60s long
if  f_duration_sec<( duration/Rs)

    duration=f_duration_sec*Rs;
    stream=round(duration/2);
    Data_dat= Dat_tracker(filename,stream,duration,CH_N);
    [S] = dat2gam(Data_dat,resample_factor,arburg_n,Rs,freq,smooth_win);
    [cc,rr]=  size(S);
    fwrite(fid_lfp,int16(reshape(S,1,cc*rr).*1000),'int16');

else
    
    for stream=s_time:duration:f_duration_sec*Rs-round(duration/2)
        
        if stream == s_time
            
            Data_dat=   Dat_tracker(filename,stream+pad_length/2,duration+pad_length,CH_N);
            
%             disp(['Data', num2str(size(Data_dat))])

            
            [S] = dat2gam(Data_dat,resample_factor,arburg_n,Rs,freq,smooth_win,gamb,lowb,highb);
            
%             disp([size(S)])

            S= S(:,1:end-pad_length/(2*resample_factor));
            %S=S (:,1:end-pad_length/(2*resample_factor));
            [cc,rr]=  size(S);

         else
             Data_dat=   Dat_tracker(filename,stream,duration+2*pad_length,CH_N);
             
%              disp([size(Data_dat)])

             
             [S] = dat2gam(Data_dat,resample_factor,arburg_n,Rs,freq,smooth_win,gamb,lowb,highb);
             
%              disp([size(S)])

                         
             S =  S (:,  3*pad_length/(2*resample_factor)+1  : end-pad_length/(2*resample_factor)  ) ;
             %S =  S (:,   pad_length/(2*resample_factor)  : end-pad_length/(2*resample_factor)  ) ;
             [cc,rr]=  size(S);
             
         end
            
            
            fwrite(fid_lfp,int16(reshape(S,1,cc*rr).*1000),'int16');
            disp (num2str(stream/Rs))
            
        
    end
    %% the remaining of the file 
    
    remaining_duration= f_duration_sec*Rs -(stream + round(duration/2));
    remaining_center=round(remaining_duration/2)+stream + s_time;

    if remaining_duration>0
        Data_dat=   Dat_tracker(filename,remaining_center-pad_length/2,remaining_duration+2*pad_length,CH_N);
        [S] = dat2gam(Data_dat,resample_factor,arburg_n,Rs,freq,smooth_win,gamb,lowb,highb);
       % S =  S (:,   pad_length/(2*resample_factor)  : end ) ;
        S =  S (:, 2*pad_length/(resample_factor) +1 : end ) ;
        [cc,rr]=  size(S);
        fwrite(fid_lfp,int16(reshape(S,1,cc*rr).*1000),'int16');
    end
    
end

fclose all;

toc

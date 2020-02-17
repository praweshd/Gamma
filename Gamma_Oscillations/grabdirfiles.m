%% Grab a certain file from all folders and do computation 

%% DETECT GAMMA IN ALL
close all
clear all
clc

%OR13/OR15
%files = [dir('*_0um_baseline*') ; dir('*_0um_posttrain*') ; dir('*_0um_pretrain*'); dir('*_0um_intertrain*')]; 

%OR6
%files = [dir('*_baseline') ; dir('*_posttrain') ; dir('*_pretrain'); dir('*_intertrain')]; 

%OR17
files = [dir('*_0um_baseline*')];
%%
for j = 1 : length(files)
    
    sesfolder = strcat(files(j).folder, '/', files(j).name);
    cd(sesfolder)
    
    if isfile(strcat(files(j).name,'_highGAM_RES.mat')) == 1
        disp(['Gamma already detected for ', files(j).name]) 
        cd ..
    else 
        disp(['Detecting Gamma now for ', files(j).name]) 
                            
        
        lfp_file = dir('*.lfp');
        CH_N=xml2CH_N(cat(2,lfp_file.name(1:end-4),'.xml'));
        [~, fbasename, ~] = fileparts(lfp_file.name);
        gamma_file = strcat(fbasename, '_highGAM_RES'); 
        
        mkdir EventFiles
        
        tic

        parfor i = 1 : 128

            [gamma] = Gamma_Detect(CH_N, i); 

            gammas(i).res = gamma; 

            % gamma_file = strcat(fbasename, '_', num2str(i), 'gamma_wav'); 
            % save(gamma_file, 'gamma');
            gamma_events = strcat(fbasename,'_', num2str(i), '.gam.evt');
            channelID = i-1; 
            SaveRippleEvents(gamma_events, gamma, channelID);
            movefile([strcat(fbasename,'_', num2str(i), '.gam.evt')], strcat(pwd,'/','EventFiles')) 

            
        end 

        toc
        save(gamma_file, 'gammas'); 
 
        clear gammas
        cd ..        
        
    end     
end 

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do CCG in all 

close all
clear all 
clc

%Or13
% files = [dir('*_0um_baseline*') ; dir('*_0um_posttrain*') ; dir('*_0um_pretrain*'); dir('*_0um_intertrain*')]; 


%OR6
files = [dir('*_baseline') ; dir('*_posttrain') ; dir('*_pretrain'); dir('*_intertrain')]; 

%%
for j = 1:length(files)
    
    close all    
    sesfolder = strcat(files(j).folder, '/', files(j).name);
    cd(sesfolder)
    
    if isfile(strcat(files(j).name,'__bad_CH.mat')) == 0
        
        disp(['Creating bad_ch file for ', files(j).name])
        good_file = dir('*_goodCH.mat');
        load(good_file.name); 
        createbadchannels([setdiff(1:128,good_CH+1)]-1); 
%         createbadchannels([]);         
        
    end 
    
    pp = pwd; 
    ani = pp(43:49);
    %Load CH
    lfp_file = dir('*.lfp');
    [~, fbasename, ~] = fileparts(lfp_file.name);
    CH_Nall=xml2CH_N(cat(2,lfp_file.name(1:end-4),'.xml'));

    
    if isfile([ani,' ',fbasename, '_GlobalNORM.png']) == 1
        disp(['CCG already detected for ', files(j).name]) 
        cd ..
    else 
        disp(['Plotting CCG for ', files(j).name]) 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Load the states
        state_mat = dir('*-states*');
        load (state_mat.name);
        StateIntervals = ConvertStatesVectorToIntervalSets(states);                 % 6 Intervalsets representing sleep states
        REM = StateIntervals{5}; 
        NREM = StateIntervals{3};
        WAKE = StateIntervals{1};
        state = NREM;

        %Load detected event 
        detec_file = dir('*_highGAM_RES.mat');
        %detec_file = dir('*_highGAM_res_II.mat');
        load(detec_file.name); 

        %Load Bad channels
        badch = dir('*bad_CH*');
        load(badch.name);
        CH_noise = bad_ch'; 

        NREM_dur = length(find(states==3)); 
        Rs = 1250;

%         %Load grid info - OLD 
%         NG = load('/home/prawesh/Documents/Insync/MATLAB_codes/Gamma_Oscillations/NG2019_3_info.mat'); 
%         NGmap = NG.NG_2019_v3.map; 
%         NGdim = NG.NG_2019_v3.dim;  

        %New Grid for OR15
        NG = load('/home/prawesh/Documents/Insync/MATLAB_codes/Gamma_Oscillations/NG2019_4_info.mat'); 
        NGmap = NG.NG_2019_v4.map; 
        NGdim = NG.NG_2019_v4.dim;
        
 
        for ref = 1 :  length(gammas)

            res1= ((gammas(ref).res(:,1)) + (gammas(ref).res(:,3)))./2;
            disp(['CCG for CH ', num2str(ref)])

            for chs = 1:length(gammas)

                res2= ((gammas(chs).res(:,1)) + (gammas(chs).res(:,3)))./2;

                window=500;
                bin=10;

                [H, B, hiB, loB]=CCG_res(res1,res2,window,bin,Rs); 
                GammaCCG(ref).zerobin(chs) =  H(B==0);

            end
        end 

        gamma_file = strcat(fbasename, '_CCG'); 
        save(gamma_file, 'GammaCCG');

        % Plot CCG, create global CCG

        ccgfile = dir('*_CCG.mat');
        load(ccgfile.name);
        
        mkdir CCG 

        for i = 1 : length(GammaCCG)
            if ismember(i,CH_noise) || ~ismember(i,NGmap)
                disp(['No CCG for bad ch ', num2str(i)])
                i = i +1;                              
            else 

                close all
                ccgval = GammaCCG(i).zerobin';
                ccgval(bad_ch) = NaN; 
                ccgval = ccgval(NGmap,1);                                                    
                grid_ccgX=reshape(ccgval,NGdim(1),NGdim(2));        

                fig=figure_ctrl('grid ccg',600,800);
                heatmap(grid_ccgX','CellLabelColor','none');   
        %         imagesc(grid_ccgX')
                colormap jet;  colorbar; 
                caxis([min(ccgval)  max(ccgval)])
                title(['High GammaCCG ', num2str(i)])
                
                %savethis = strcat(pwd,'/','CCG','/',fbasename, '_CCG_', num2str(j),'.png');

                print([fbasename, '_CCG_', num2str(i)], '-dpng')
                movefile([fbasename, '_CCG_', num2str(i),'.png'], strcat(pwd,'/','CCG')) 

            end
        end

        %Global
        close all
        
        ccgfile = dir('*_CCG.mat');
        load(ccgfile.name);

        globaltest = zeros(length(GammaCCG),1);

        for i = 1:length(GammaCCG) 
            indccg = GammaCCG(i).zerobin';
            indccg(i) = 0;
            globaltest = globaltest + indccg;   
        end 

        globaltest(bad_ch) = NaN; 
        gridtestt = globaltest(NGmap,1);
        grid_gob = reshape(gridtestt,NGdim(1),NGdim(2));

        gridtestnorm = normalize(gridtestt,'range');
        grid_gobnorm = reshape(gridtestnorm,NGdim(1),NGdim(2));

        ccgnorm.val = [gridtestnorm]; 
        ccgnorm.sum = nansum(gridtestnorm);
        ccgnorm_file = strcat(fbasename, '_CCG_VAL'); 
        save(ccgnorm_file, 'ccgnorm');

        fig=figure_ctrl('glob ccg',600,800); 
        heatmap(grid_gob','CellLabelColor','none');  
        colormap jet;  colorbar; 
        caxis([0 max(gridtestt)])
        title([ani,' ',fbasename, ' Global'])
        print([ani,' ',fbasename, '_Global'], '-dpng')

        fig=figure_ctrl('glob ccg norm',600,800); 
        heatmap(grid_gobnorm','CellLabelColor','none');  
        colormap jet;  colorbar; 
        caxis([0 1])
        title([ani,' ',fbasename, ' Global'])
        print([ani,' ',fbasename, '_GlobalNORM'], '-dpng')
        
        cd ..
    
    end

end 

%% Grab and move CCG file to one folder
close all
clear all
clc

%OR13/15
% files = [dir('*_0um_baseline*') ; dir('*_0um_posttrain*') ; dir('*_0um_pretrain*'); dir('*_0um_intertrain*')]; 

%%
%OR6
files = [dir('*_baseline') ; dir('*_posttrain') ; dir('*_pretrain'); dir('*_intertrain')]; 
%%

pdest = strcat(pwd,'/BehaviorSPI');

for j = 1:length(files) 
    
    close all    
    sesfolder = strcat(files(j).folder, '/', files(j).name);
    cd(sesfolder)
    
    globcc = dir('*_SPIGlobalNORM*');
    
    if length(globcc) > 0
        disp([globcc.name])
        source = fullfile(pwd,globcc.name);
        destination = fullfile(pdest,globcc.name);
        copyfile(source,destination)
        disp(['Copied CCG fig for ', files(j).name]) 
    end
    
    cd ..

end  


%% Plot PNG in one file
close all
clear all
clc
fig=figure_ctrl('Gamma_PGD',400,2000);
f= dir('*_GlobalNORM*'); 


for i = 1: length(f)
    
   subaxis(11,3,i,'SpacingVert',0.01,'SpacingHoriz', 0.01, 'ML', 0.01, 'MR',0.01, 'MT', 0.01, 'MB', 0.01)
   img = imread(f(i).name);
   imagesc(img); axis off
   gname = f(i).name(13:20);
   sname = f(i).name(13:24); 
   title(sname,'FontSize', 5, 'Interpreter', 'None')    
    
end 

print('OR15_Gammall', '-depsc', '-tiff','-r300', '-painters')

%%

k=0;
days = [14:29];
fig=figure_ctrl('Gamma_PGD',400,2000);

for i = 1:length(days)
   
  test = ['*08',num2str(days(i)),'*'];
     
  filename = dir(test);
  k=k+1;
  subaxis(16,3,k,'SpacingVert',0.01,'SpacingHoriz', 0.01, 'ML', 0.01, 'MR',0.01, 'MT', 0.01, 'MB', 0.01)
  
  img = imread(filename(end).name);
  imagesc(img); axis off
  title(test, 'FontSize', 5)
  
  k=k+1;
  subaxis(16,3,k,'SpacingVert',0.01,'SpacingHoriz', 0.01, 'ML', 0.01, 'MR',0.01, 'MT', 0.01, 'MB', 0.01)
  img = imread(filename(end-1).name);
  imagesc(img); axis off 
  
  
  k=k+1;
     
  if length(filename)==3
      
    subaxis(16,3,k,'SpacingVert',0.01,'SpacingHoriz', 0.01, 'ML', 0.01, 'MR',0.01, 'MT', 0.01, 'MB', 0.01)
    img = imread(filename(end-2).name);
    imagesc(img); axis off     
      
  end 
     
    
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAKE ONE BIG STRUCT FILE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do CCG in all 

close all
clear all 
clc
%%

%OR13/15
files = [dir('*_0um_baseline*') ; dir('*_0um_posttrain*') ; dir('*_0um_pretrain*'); dir('*_0um_intertrain*')]; 
%%
% %OR6
%files = [dir('*_baseline') ; dir('*_posttrain') ; dir('*_pretrain'); dir('*_intertrain')]; 

for j = 1 : length(files)

    close all    
    sesfolder = strcat(files(j).folder, '/', files(j).name);
    cd(sesfolder)
    disp(['Saving in struct for ', files(j).name]) 
    
    load((strcat(files(j).name,'_highGAM_RES.mat')))
    higamma = gammas;
    
    %Load the states
    state_mat = dir('*-states*');
    load (state_mat.name);
    StateIntervals = ConvertStatesVectorToIntervalSets(states);                
    REM = StateIntervals{5}; 
    NREM = StateIntervals{3};
    WAKE = StateIntervals{1};
    state = NREM;
    pp = pwd; 
    ani = pp(43:49);

    %Load CH
    lfp_file = dir('*.lfp');
    [~, fbasename, ~] = fileparts(lfp_file.name);
    CH_Nall=xml2CH_N(cat(2,lfp_file.name(1:end-4),'.xml'));

    %Load Bad channels
    badch = dir('*bad_CH*');
    load(badch.name);
    CH_noise = bad_ch'; 

    NREM_dur = length(find(states==3)); 
    Rs = 1250;

%         %Load grid info - OLD 
%         NG = load('/home/prawesh/Documents/Insync/MATLAB_codes/Gamma_Oscillations/NG2019_3_info.mat'); 
%         NGmap = NG.NG_2019_v3.map; 
%         NGdim = NG.NG_2019_v3.dim;  

    %New Grid for OR15
    NG = load('/home/prawesh/Documents/Insync/MATLAB_codes/Gamma_Oscillations/NG2019_4_info.mat'); 
    NGmap = NG.NG_2019_v4.map; 
    NGdim = NG.NG_2019_v4.dim;  

    
    for ref = 1 :  length(gammas)

        res1= ((gammas(ref).res(:,1)) + (gammas(ref).res(:,3)))./2;
%         disp(['CCG for CH ', num2str(ref)])

        for chs = 1:length(gammas)

            res2= ((gammas(chs).res(:,1)) + (gammas(chs).res(:,3)))./2;

            window=500;
            bin=10;

            [H, B, hiB, loB]=CCG_res(res1,res2,window,bin,Rs); 

            ccg(chs).H = [H];
            ccg(chs).B = [B];
            zerobin(chs) =  H(B==0);

        end

        higamma(ref).ccgH = [ccg(:).H];
        higamma(ref).ccgB = [ccg(:).B];
        higamma(ref).zerobin = [zerobin]; 

    end 
    
    ccgfile = dir('*_CCG.mat');
    load(ccgfile.name);

    for i = 1 : length(GammaCCG)
        if ismember(i,CH_noise) || ~ismember(i,NGmap)
%             disp(['No CCG for bad ch ', num2str(i)])
            i = i +1;                              
        else 

            ccgval = GammaCCG(i).zerobin';
            ccgval(bad_ch) = NaN; 
            ccgval = ccgval(NGmap,1);                                                    
            grid_ccgX=reshape(ccgval,NGdim(1),NGdim(2));  
            higamma(i).ccgheat = [ccgval];        
        end
    end

    save([strcat(files(j).name,'_summary')],'higamma')

    ORbehavior(j).rat = ani;
    ORbehavior(j).ses = fbasename; 
    ORbehavior(j).gam = [higamma];    

    load((strcat(files(j).name,'_CCG_VAL.mat')))

    ORbehavior(j).glob = [ccgnorm];

    cd ..
end

save([strcat(ani,'_summary')],'ORbehavior')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SPINDLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DETECT SPI IN ALL
close all
clear all
clc

%OR13/15
files = [dir('*_0um_baseline*') ; dir('*_0um_posttrain*') ; dir('*_0um_pretrain*'); dir('*_0um_intertrain*')]; 

% %OR6
% files = [dir('*_baseline') ; dir('*_posttrain') ; dir('*_pretrain'); dir('*_intertrain')]; 
%% Spindle Wavelet Detect

for j = 1 : length(files)
    
    sesfolder = strcat(files(j).folder, '/', files(j).name);
    cd(sesfolder)
    
    if isfile(strcat(files(j).name,'_SPI_RES.mat')) == 1
        disp(['Spindle already detected for ', files(j).name]) 
        cd ..
    else 
        disp(['Detecting Spi now for ', files(j).name])                            
        
        lfp_file = dir('*.lfp');
        CH_N=xml2CH_N(cat(2,lfp_file.name(1:end-4),'.xml'));
        [~, fbasename, ~] = fileparts(lfp_file.name);
        spi_file = strcat(fbasename, '_SPI_RES'); 
        
        tic
        parfor i = 1 : 128

            [spindles] = Spindle_Detect_wavelet(CH_N, i)        

            spindle(i).res = spindles; 

            spindle_events = strcat(fbasename,'_', num2str(i), '.spi.evt');
            channelID = i-1; 
            SaveRippleEvents(spindle_events, spindles, channelID);

            mkdir SPI_Events
            movefile([strcat(fbasename,'_', num2str(i), '.spi.evt')], strcat(pwd,'/','SPI_Events')) 
        end 

        toc
        save(spi_file, 'spindle');
        clear spindle
        cd ..          
    end     
end 

%% SpindleCCG in all
for j = 1 : length(files)
    
    close all    
    sesfolder = strcat(files(j).folder, '/', files(j).name);
    cd(sesfolder)
    
    if isfile(strcat(files(j).name,'__bad_CH.mat')) == 0
        
        disp(['Creating bad_ch file for ', files(j).name])
        good_file = dir('*_goodCH.mat');
        load(good_file.name); 
        createbadchannels([setdiff(1:128,good_CH+1)]-1); 
%         createbadchannels([]);       
    end 
    
    if isfile(strcat(files(j).name,'_spiCCG.mat')) == 1
        disp(['Spi CCG already detected for ', files(j).name]) 
        cd ..
    else 
        disp(['Detect and plot spi CCG for ', files(j).name]) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Load the states
        state_mat = dir('*-states*');
        load (state_mat.name);
        StateIntervals = ConvertStatesVectorToIntervalSets(states);                 % 6 Intervalsets representing sleep states
        REM = StateIntervals{5}; 
        NREM = StateIntervals{3};
        WAKE = StateIntervals{1};
        state = NREM;
        pp = pwd; 
        ani = pp(43:49);

        %Load CH
        lfp_file = dir('*.lfp');
        [~, fbasename, ~] = fileparts(lfp_file.name);
        CH_Nall=xml2CH_N(cat(2,lfp_file.name(1:end-4),'.xml'));
        
        %Load detected event 
        detec_file = dir('*_SPI_RES.mat');
        %detec_file = dir('*_highGAM_res_II.mat');
        load(detec_file.name); 

        %Load Bad channels
        badch = dir('*bad_CH*');
        load(badch.name);
        CH_noise = bad_ch'; 

        NREM_dur = length(find(states==3)); 
        Rs = 1250;
        
        %Load OLD grid info 
        NG = load('/home/prawesh/Documents/Insync/MATLAB_codes/Gamma_Oscillations/NG2019_3_info.mat'); 
        NGmap = NG.NG_2019_v3.map; 
        NGdim = NG.NG_2019_v3.dim;

        %Regular CCG        
        for ref = 1 :length(spindle)

            res1= ((spindle(ref).res(:,1)) + (spindle(ref).res(:,3)))./2;

            if length(res1) > 100
                disp(['CCG for CH ', num2str(ref)])

                for chs = 1:length(spindle) 
                    
                    res2= ((spindle(chs).res(:,1)) + (spindle(chs).res(:,3)))./2;
                    if length(res2) > 100 
                        window=30000;
                        bin=200;

                        [H, B, hiB, loB]=CCG_res(res1,res2,window,bin,Rs); 
                        SpiCCG(ref).zerobin(chs) = H(B==0);
                    else
                        SpiCCG(ref).zerobin(chs) = 0;         
                
                    end 

                end    
            else
                SpiCCG(ref).zerobin = [];
                disp(['No CCG for CH ', num2str(ref)])        
            end   
        end 
        
        gamma_file = strcat(fbasename, '_spiCCG'); 
        save(gamma_file, 'SpiCCG');             

        ccgfile = dir('*_spiCCG.mat');
        load(ccgfile.name);
        
        mkdir CCG_SPI
        
        close all
         
        for i = 1 : length(SpiCCG)
            if ismember(i,CH_noise) || ~ismember(i,NGmap)
                disp(['No CCG for bad ch ', num2str(i)])
                i = i +1;                              
            else 

                close all
                
                ccgval = SpiCCG(i).zerobin';
                if length(ccgval) > 0
                    ccgval(bad_ch) = NaN; 
                    ccgval = ccgval(NGmap,1);                                                    
                    grid_ccgX=reshape(ccgval,NGdim(1),NGdim(2));        

                    fig=figure_ctrl('grid ccg',600,800);
                    heatmap(grid_ccgX','CellLabelColor','none');   
                    colormap jet;  colorbar; 
                    caxis([min(ccgval)  max(ccgval)])
                    title(['SpiCCG ', num2str(i)])

                    print([fbasename, '_spiCCG_', num2str(i)], '-dpng')
                    movefile([fbasename, '_spiCCG_', num2str(i),'.png'], strcat(pwd,'/','CCG_SPI'))
                end
            end
        end       
        
        %Global
        close all         
        ccgfile = dir('*_spiCCG.mat');
        load(ccgfile.name);

        globaltest = zeros(length(SpiCCG),1);
        for i = 1:length(SpiCCG)     
            indccg = SpiCCG(i).zerobin';
            if length(indccg) > 0
                indccg(i) = 0;
                globaltest = globaltest + indccg;
            end
        end 
 
        globaltest(bad_ch) = NaN; 
        gridtestt = globaltest(NGmap,1);
        if nansum(gridtestt) > 0
            grid_gob = reshape(gridtestt,NGdim(1),NGdim(2));

            spigridnorm = normalize(gridtestt,'range');
            grid_gobnorm = reshape(spigridnorm,NGdim(1),NGdim(2));

            spiccgnorm.val = [spigridnorm]; 
            spiccgnorm.sum = nansum(spigridnorm);
            ccgnorm_file = strcat(fbasename, '_spiCCG_VAL'); 
            save(ccgnorm_file, 'spiccgnorm');

            fig=figure_ctrl('glob ccg',600,800); 
            heatmap(grid_gob','CellLabelColor','none');  
            colormap jet;  colorbar; 
            caxis([min(gridtestt) max(gridtestt)])
            title([ani,' ',fbasename, ' SPIGlobal'])
            print([ani,' ',fbasename, '_SPIGlobal'], '-dpng')

            fig=figure_ctrl('glob ccg norm',600,800); 
            heatmap(grid_gobnorm','CellLabelColor','none');  
            colormap jet;  colorbar; 
            caxis([0 1])
            title([ani,' ',fbasename, ' SPIGlobalnorm'])
            print([ani,' ',fbasename, '_SPIGlobalNORM'], '-dpng')
        end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear fbasename
        cd ..
    
    end    
end



%% -----------------------------------------------------------------------------------------------------------
%% Move SPI struct to one file 

close all
clear all
clc











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  
% %% OR13
% files = [dir('*_0um_posttrain*')]; 
% ptfiles = [2:7 12:17]; %or13
% 
% %% OR15
% files = [dir('*_0um_baseline*') ; dir('*_0um_posttrain*')]; 
% ptfiles = [1:5 9:14]; %or15
% 
% %% OR6
% files = [dir('*_posttrain*')]; 
% ptfiles = [2:4 8:13]; %or6
% 
% %%
% 
% for i = 1 : length(ptfiles)
%     
%     filenum = ptfiles(i);
%     ptfolder = strcat(files(filenum).folder, '/', files(filenum).name);
%     cd(ptfolder)
%     ccgfile = dir('*_CCG_VAL*');
%     load(ccgfile.name);
%     ccgbv(i) = ccgnorm.sum;
%     clear ccgnorm
%     cd .. 
%     
%     
% end 
% 
% %% Plot OR13
% close all
% test = reshape(ccgbv, 3,4);
% fig=figure_ctrl('Occ Rate',400,600);
% shadedErrorBar([], mean(test,1), std(test),'b');
% ylim([50 100])
% print(['OR13_sumccg'], '-depsc', '-tiff','-r300', '-painters')
% %%
% fig=figure_ctrl('Occ Rate',400,600);
% boxplot(test,'Labels',{'Baseline','Cheeseboard','Object','Object-Space'})
% ylim([65 95])
% print(['OR13_sumccgII'], '-depsc', '-tiff','-r300', '-painters')
% 
% %% Plot OR15
% 
% tests(1,1) = NaN;
% tests(1,2:12) = ccgbv;
% test = reshape(tests, 3,4);
% fig=figure_ctrl('Occ Rate',400,600);
% boxplot(test,'Labels',{'Baseline','Cheeseboard','Object','Object-Space'})
% ylim([34 50])
% print(['OR15_sumccgII'], '-depsc', '-tiff','-r300', '-painters')
% 
% %% Plot OR6
% close all
% test = reshape(ccgbv, 3,3);
% fig=figure_ctrl('Occ Rate',400,600);
% shadedErrorBar([], mean(test,1), std(test),'b');
% ylim([30 100])
% print(['OR6_sumccg'], '-depsc', '-tiff','-r300', '-painters')
% %%
% fig=figure_ctrl('Occ Rate',400,600);
% boxplot(test,'Labels',{'Cheeseboard','Object','Object-Space'})
% ylim([30 100])
% print(['OR6_sumccgII'], '-depsc', '-tiff','-r300', '-painters')
% 
% 
% %%
% close all
% 
% test = reshape(ccgbv, 3, 4)';
% X = categorical();
% X = reordercats(X,{'Baseline','Cheeseboard','Object','Object-Space'});
% 
% fig=figure_ctrl('Occ Rate',400,600);
% bar(X,test)
% ylim([0 60])
% ylabel('Normalized CCG sum')
% 
% 
% 

%% Testing CCG between detected Gamma events
% Prawesh Dahal
%June 24

close all
clear all
clc
 
%%
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

%%
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

%%
%Load grid info 
NG = load('/home/prawesh/Documents/Insync/MATLAB_codes/Gamma_Oscillations/NG2019_3_info.mat'); 
 
NGmap = NG.NG_2019_v3.map; 
NGdim = NG.NG_2019_v3.dim;


%% Regular CCG
 clc
for ref = 1 :length(spindle)

    res1= ((spindle(ref).res(:,1)) + (spindle(ref).res(:,3)))./2;
    
    if length(res1) > 100
        disp(['CCG for CH ', num2str(ref)])

        for chs = 1:length(spindle)
            disp(['CCG for refCH ', num2str(ref),' for resch', num2str(chs)])

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

%%
gamma_file = strcat(fbasename, '_spiCCG'); 
save(gamma_file, 'SpiCCG');

%% Plot CCG, create global CCG
 
ccgfile = dir('*_spiCCG.mat');
load(ccgfile.name);

%CH_noise = []; bad_ch =[];

close all
clc
for i = 1 : length(SpiCCG)
    if ismember(i,CH_noise) || ~ismember(i,NGmap)
        disp(['No CCG for bad ch ', num2str(i)])
        i = i +1;                              
    else 
        
        close all
        ccgval = SpiCCG(i).zerobin';
        ccgval(bad_ch) = NaN; 
        ccgval = ccgval(NGmap,1);                                                    
        grid_ccgX=reshape(ccgval,NGdim(1),NGdim(2));        
        
        fig=figure_ctrl('grid ccg',600,800);
        heatmap(grid_ccgX','CellLabelColor','none');   
%         imagesc(grid_ccgX')
        colormap jet;  colorbar; 
        caxis([min(ccgval)  max(ccgval)])
        title(['SpiCCG ', num2str(i)])
         
        print([fbasename, '_spiCCG_', num2str(i)], '-dpng')
%         print(['CCGheat_', num2str(i)], '-depsc', '-tiff','-r300', '-painters') 


    end
end

 %%
%Global
close all
clc
ccgfile = dir('*_spiCCG.mat');
load(ccgfile.name);

globaltest = zeros(length(SpiCCG),1);


for i = 1:length(SpiCCG)  
    i
    indccg = SpiCCG(i).zerobin';
    if length(indccg) > 0
        indccg(i) = 0;
        globaltest = globaltest + indccg;
    end
end 

%%
globaltest(bad_ch) = NaN; 
gridtestt = globaltest(NGmap,1);
grid_gob = reshape(gridtestt,NGdim(1),NGdim(2));

gridtestnorm = normalize(gridtestt,'range');
grid_gobnorm = reshape(gridtestnorm,NGdim(1),NGdim(2));

ccgnorm.val = [gridtestnorm]; 
ccgnorm.sum = nansum(gridtestnorm)
ccgnorm_file = strcat(fbasename, '_spiCCG_VAL'); 
save(ccgnorm_file, 'ccgnorm');

fig=figure_ctrl('glob ccg',600,800); 
heatmap(grid_gob','CellLabelColor','none');  
colormap jet;  colorbar; 
caxis([min(gridtestt) max(gridtestt)])
title([ani,' ',fbasename, ' Global'])
print([ani,' ',fbasename, '_Global'], '-dpng')

fig=figure_ctrl('glob ccg norm',600,800); 
heatmap(grid_gobnorm','CellLabelColor','none');  
colormap jet;  colorbar; 
caxis([0 1])
title([ani,' ',fbasename, ' Global'])
print([ani,' ',fbasename, '_GlobalNORM'], '-dpng')


%% IED Population Analysis

close all
clear all
clc

NG = load('/home/prawesh/Documents/Insync/MATLAB_codes/Gamma_Oscillations/NG2019_3_info.mat'); 
 
NGmap = NG.NG_2019_v3.map; 
NGdim = NG.NG_2019_v3.dim;

Rs = 1250; 
%% 

days = [6 8:14 16:17];

for day = 1:length(days)
iedstr = dir(['*_K',num2str(days(day)),'_post*']);
load(iedstr.name);
disp(['day',num2str(days(day))])
allCH = [IED_all.channel];

%%

for refchs = 1 : length(NGmap)
    
REFCH = NGmap(refchs);
ref = IED_all(find(allCH == REFCH)).res;
ref_samp = round(ref*Rs);


for i = 1:length(NGmap)
    
    RESCH = NGmap(i);
    if REFCH ~= RESCH 
    
        res  = IED_all(find(allCH == RESCH)).res;
        res_samp = round(res*Rs);
        in = [];

        for k = 1:length(ref)

            for j = 1:length(res)            

                tmp = abs(res_samp(j) - ref_samp(k));
                if tmp < 0.05*Rs
                    if ismember(res(j),in) == 0
                        in = [in res(j)];    
                    end
                end 
            end            
        end  


        matched = in;
        unmatched = setdiff(res,matched);
        IEDpop(i).refCH = REFCH; 
        IEDpop(i).resCH = RESCH;
        IEDpop(i).matched = matched;
        IEDpop(i).unmatched = unmatched; 
    
    else
        IEDpop(i).refCH = REFCH; 
        IEDpop(i).resCH = RESCH;
        IEDpop(i).matched = [];
        IEDpop(i).unmatched = []; 
        
    end
    
    indIED(days(day)).day = days(day); 
    indIED(days(day)).ref(REFCH).chs = [IEDpop];
end
%  disp(['Saving K17_pop_CH_',])
%  save(['K17_pop_CH_', num2str(REFCH)],'IEDpop')
end  

end

save('RatC17_indIED','indIED')

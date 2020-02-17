%% IED spatial distribution in NG
% Prawesh Dahal
% Nov 20, 2019

close all 
clear all
clc

day = 1; 
iedstr = dir(['*_K',num2str(day),'_post*']);
load(iedstr.name);
Rs = 1250;

%%
%Load Rat name Filename
filename = dir('B*_HC_*');
animal = filename.name(1:3); 
chs = {IED_all.reg};

%Find hippocampal CHs
HC_ind = find(strcmp(chs, 'HC'));
mpfc = find(strcmp(chs, 'mPFC'));

NG = load('/home/prawesh/Documents/Insync/MATLAB_codes/Gamma_Oscillations/NG2019_3_info.mat'); 

NGmap = NG.NG_2019_v3.map; 
NGdim = NG.NG_2019_v3.dim;
allCH = [IED_all.channel];

 
%% Remove IEDs < 20ms
if ~isempty(HC_ind) 
     HC_time = []; 
    for i = 1:length(HC_ind)     
        HC_time = [HC_time; IED_all(HC_ind(i)).res];
    end  
     
    HC_time = sort(HC_time);
    HC_samp = round(HC_time*Rs);
    HC_samp = sort(HC_samp); 
 
    if ~isempty(HC_samp) 
        HC_res = []; j=0;         
        for i = 1 : (length(HC_samp)-1)
            A = HC_samp(i+1); 
            B = HC_samp(i);
            tmp = A - B;
            if tmp >= 0.02*Rs
                j = j+1;
                HC_res(j) = HC_time(i);             
            end 
        end 
        HC_res = HC_res'; 
    else
     HC_res = [];    
    end     
else 
    disp('No hippocampal CHs in this struct >.<')
    HC_res = [];
end     
  

% Find number of IEDs in each channel over day

load('IEDnum_days.mat');

IEDnum(1,day) = length(HC_res); 
for i = 1 : length(NGmap)
    
    REFCH = NGmap(i);
    whereCH = find(allCH == REFCH);
    if ~isempty(whereCH)
        IEDnum(i+1,day) = length(IED_all(whereCH).res);
    else
        IEDnum(i+1,day) = 0;
    end
    
end 

save(['IEDnum_days'],'IEDnum')

%% Plot
close all
clc
load('IEDnum_days.mat');

daysnum = 1:19;
hipp = IEDnum(1,:);
parulac = parula(19); 
grayc = gray(10);
figure_ctrl('grid ccg travel',800,800);
for i = 1:19
scatter(ones(117,1)*i, IEDnum(2:end,i), 'MarkerFaceColor',grayc(7,:),'MarkerEdgeColor',grayc(7,:),'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.6)
hold on :
scatter(daysnum(i), IEDnum(1,i),'filled', 'MarkerFaceColor',parulac(11,:),'MarkerEdgeColor',parulac(11,:),'MarkerFaceAlpha',.9,'MarkerEdgeAlpha',.6)
line(daysnum, hipp)
errorbar(i, mean(IEDnum(2:end,i)), std(IEDnum(2:end,i))./sqrt(118),'r'); 
scatter(daysnum(i), mean(IEDnum(2:end,i)),'filled', 'MarkerFaceColor',parulac(1,:),'MarkerEdgeColor',parulac(1,:),'MarkerFaceAlpha',.9,'MarkerEdgeAlpha',.6)
ylabel('Number of IEDs detected in channel')
xlabel('Kindling days') 
title('Number of IEDs in HP and NG over kindling days')
end
print(['IEDNUM'], '-dpng') 


%%
savename1= ['RAT_',num2str(animal),'day_',num2str(day),'_per_',num2str(length(HC_res))];
savename2= ['RAT_',num2str(animal),'day_',num2str(day),'_del_',num2str(length(HC_res))];



load('plotmap.mat') 

left = find(plotmap==0);
right = find(plotmap==1);
mpfcloc = find(plotmap==2); 



figure_ctrl('grid ccg travel',2000,2000);

for i = 1 : length(NGmap)

    REFCH = NGmap(i);
    whereCH = find(allCH == REFCH);

    subplot(13,15,left(i))
    if isempty(whereCH)
        B=0;
        H=0;
        b = bar(B,H);
        netdelay(i) = NaN; sigper(i) = NaN; sigval(i) = NaN;    

    elseif ~isempty(whereCH)
        NGres = IED_all(whereCH).res;
        window=500;
        bin=10;        
        if isempty(NGres)
            B=0;
            H=0;
            b = bar(B,H);
            netdelay(i) = NaN; sigper(i) = NaN; sigval(i) = NaN;            

        else
            [H, B, hiB, loB]=CCG_res(HC_res,NGres,window,bin,Rs);
            b = bar(B,H);
            hold on; plot(B,hiB,'--r'); 
            hold on; plot(B,loB,'--r');
            %set(gca,'YTickLabel',[]);
            xlim([- 100 100]); ylim([0 5])

            cchEvt = (H.*length(HC_res)*bin/1000)';

            if hiB(find(B==0)) == 0
                sigper(i) = 0 ;
            elseif hiB(find(B==0)) ~= 0
                if H(find(B==0)) > hiB(find(B==0))

                    sigper(i) = cchEvt(find(B==0))/length(NGres);
                else
                    sigper(i) = 0 ;
                end 
            end

            %New net delay
            tmp = intersect(find(B>0), find(B<30));
            if sum(hiB) ~=0
                tmp = intersect(tmp, find(H>hiB));

                sigval(i) = sum(H(tmp));  
            elseif sum(hiB) == 0
                sigval(i) = 0; 
            end
        end  
    end 
end 

sigpergridm = reshape(sigper',9,13)';
sigvalgridm = reshape(sigval',9,13)';


subplot(13,15,right)
heatmap(sigvalgridm); colormap jet;  colorbar; caxis([0 4]) 
title([savename2])

for i = 1:length(mpfc)    

    whereCH = mpfc(i);
    mpfcres = IED_all(whereCH).res;

    subplot(13,15,mpfcloc(i))
    [H, B, hiB, loB]=CCG_res(HC_res,mpfcres,window,bin,Rs);
    b = bar(B,H);
    hold on; plot(B,hiB,'--r'); 
    hold on; plot(B,loB,'--r');
    %set(gca,'YTickLabel',[]);
    xlim([- 100 100]);
    ylim([0 5])


    %Zero-bin percentage
    cchEvt = (H.*length(HC_res)*bin/1000)';
    if H(find(B==0)) > hiB(find(B==0))
        sigperm(i) = cchEvt(find(B==0))/length(mpfcres);
    else 
        sigperm(i) = 0; 
    end

    %New net delay
    tmp = intersect(find(B>0), find(B<30));
    tmp = intersect(tmp, find(H>hiB)); 

    sigvalm(i) = sum(H(tmp));             

end

subplot(13,15,mpfcloc(5:8))
% netdelaygridm = reshape(netdelaym',4,1)';
sigvalgridm = reshape(sigvalm',4,1)';

% heatmap(netdelaygridm','CellLabelColor','none'); colormap jet;  colorbar; caxis([-2  2]) 
heatmap(sigvalgridm','CellLabelColor','none'); colormap jet;  colorbar; caxis([0 4]) %caxis([-2  2]) 
print([savename2], '-dpng') 

figure_ctrl('grid ccg',1000,1000);
heatmap(sigpergridm); colormap jet;  colorbar; caxis([0 0.8])  
print([savename1], '-dpng') 
title([savename1])



 
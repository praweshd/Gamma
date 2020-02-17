%% YOU MUST ENTRE THE CORRECT FILM THICKNESS !!!!
clc
clear
close all
IV_file=dir('*.csv');
%%
% IV file strut: Repeat	VAR2 Point VD_Voltage VD_Current VD_Time VG_Voltage	VG_Current VG_Time
% 5 12 30 110 250 500
size_vector = [5 12 30 110 250 500].*1e-6;

for i=1:length (IV_file)
   filename=IV_file(i).name;
   data=csvread(filename,1,0);
   Vg_n=max(data(:,2)) ;
   Vd_n=max(data(:,3)) ;
   
   IV(i).VD= data(1:Vd_n,4);
   IV(i).VG= round ( data(1:Vd_n:end,7) ,2);
   IV(i).ID = reshape ( data(:,5) ,Vd_n, Vg_n);
   IV(i).IG=  reshape ( data(:,8) ,Vd_n, Vg_n);
   
   IV(i).transfer = IV(i).ID (end,:);
   IV(i).Gm  = abs ( diff(IV(i).transfer) ./ abs(IV(i).VG(2) - IV(i).VG(1)) );
   IV(i).onoff = abs (min (IV(i).transfer) / max (IV(i).transfer));
   
   IV(i).corrd = [str2num(filename(2)) str2num(filename(4))];
   
   IV(i).set =   filename(1);
   IV(i).L= size_vector (IV(i).corrd (1));
   IV(i).W= size_vector (IV(i).corrd (2));
   IV(i).d= 2e-6;

   IV(i).volume = IV(i).L * IV(i).W.* IV(i).d ; 
   IV(i).wd_L = (IV(i).W.* IV(i).d)/ IV(i).L  ; 
   
end

%%
close all

fig=figure_ctrl (cat(2,'Summary IV curves of ',IV(1).set, 'in µA'), 1700,1000);
for i =1:length (IV)
    plot_location = ( (IV(i).corrd(1) -1)*6 + IV(i).corrd(2))*2 -1  ;
    
    subaxis(6,12,plot_location,'Spacing',0.01,'Padding',0.01,'Margin',0.01)

    color_vector=[ (1:Vg_n*5)./(Vg_n)];
    
    IID =IV(i).ID.*1e6;
    
    for ii =1:Vg_n
        plot ( IV(i).VD ,IID(:,ii),'Color',shade([0 1 1],color_vector(ii)) ,'LineWidth', 2); hold on    ;
    end
    xlim([-0.6 0]);
    ylim ([min(min (IID))*1.1   5 ]);
    set (gca,'Ydir','reverse');
    set (gca,'Xdir','reverse');
    gm_max=  (max (IV(i).Gm*1e3));
    subaxis(6,12,plot_location+1,'Spacing',0.01,'Padding',0.01,'Margin',0.01)
    
    yyaxis left
    plot (IV(i).VG,IV(i).transfer*1e6,'-o','MarkerSize',5,'LineWidth',2); hold on ; 
    yyaxis right
    plot (IV(i).VG(2:end),IV(i).Gm*1e3,'-o','MarkerSize',5,'LineWidth',2); 
    title (cat(2,'Gmax=', num2str(round (max(IV(i).Gm*1e3),2) ), ' OnOff=', num2str( IV(i).onoff)   ));
    
end

printeps (fig);
printjpg (fig);

%% ========================================================================================
fig = figure_ctrl ('Gm vs geometry ', 500,500);
for j=1:length (IV)
    loglog (IV(j).wd_L* 1e6,max (IV(j).Gm*1e3),'or', 'MarkerSize',8, 'LineWidth',3); hold on 
end

xlim ([0.15e-2 5e2])
printeps (fig);
printjpg (fig);

%% ========================================================================================
fig = figure_ctrl ('OnOFF vs geometry ', 500,500);
for j=1:length (IV)
    loglog (IV(j).volume*1e18,(IV(j).onoff),'or', 'MarkerSize',8, 'LineWidth',3); hold on 
end

xlim ([9 1e6])
printeps (fig);
printjpg (fig);

%% ========================================================================================
tmp=pwd;
slash_index=find (tmp=='/');
IV_filename=cat(2, tmp(slash_index(end)+1:end),'__IV');
save (IV_filename,'IV')
% Set mPFC and HC CHs with the IED detection as 'mpfc_ch_key' and 'hc_ch_key'  
% Check that 'folders_days' detect all the folders to analyze
% Run inside a rat folder with this subfolder struct rat/session_day/IED_detect (subfolder with IED detection per ch)
% Create an structure for each recording day with: 
%    1)ch; 2)region of recording ; 3)IEDs number:  4)IEDs/h of NREM
% Copy the created structure to the rat folder 



clear;
clc;

%Chs to analyze
mpfc_ch_key = [129 130 132 133]+1;
hc_ch_key = [164 165 166 188]+1;
%sessionInfo = LoadParameters (pwd);
ng1_ch_key = [70,50,83,57,76,44,69,51,82,58,75,45,68,52,81]+1;
ng2_ch_key = [37,96,7,89,0,102,38,95,6,88,32,101,39,94,5]+1;
ng3_ch_key = [109,22,122,29,115,15,108,23,121,30,114,17,107,24,120]+1;
ng4_ch_key = [59,74,46,67,53,80,60,73,47,66,54,79]+1;
ng5_ch_key = [87,33,100,40,93,4,86,34,99,41,92,3]+1;
ng6_ch_key = [31,106,10,105,113,25,112,26,18,119,19,118]+1;
ng7_ch_key = [61,72,48,65,55,78,62,71,49,64,56,77]+1;
ng8_ch_key = [85,35,98,42,91,2,84,36,97,43,90,1]+1;
ng9_ch_key = [9,111,20,104,27,117,8,110,21,103,28,116]+1;

folders_days = dir('*post*');
%folders_days ([22]) = []; 
    
 for jj = 1:length(folders_days);
     cd(folders_days(jj).name);
     disp(folders_days(jj).name);
        
        detect_bad_folders = dir ('*no_ACC_correct');
        if ~isempty (detect_bad_folders) == 1;
            for ii = 1: length (detect_bad_folders);
             movefile(detect_bad_folders.name, strcat('_', detect_bad_folders.name));
            end
        end
        
        state_mat = dir('*-states*');
        if isempty (state_mat) == 1;
        state_mat = dir('*-autostates*');
        end 
 
        detect_sub = dir('IED_detect*');
       
    if ~isempty (state_mat) == 1 &&  ~isempty (detect_sub) == 1;

        load (state_mat.name);
        DW_dur = find(states==2);
        NR_dur = find(states==3);
        NREM_duration = length(DW_dur) + length(NR_dur); 
       
        cd (detect_sub.name);
        
        prev_struc = dir ('*_IED_struct.mat');
        if ~isempty (prev_struc) == 1;
            delete (prev_struc.name);
        end
        
         IED_files = dir('*_IED.mat');
         k = 1;
         for ii = 1:length(IED_files);
             load (IED_files(ii).name);
             filename = IED_files(ii).name;
             file_index1 = strfind(filename,'Ch');
             underline_index =find(filename == '_', 4);
             underline_index = underline_index (4);
             filename = filename (1:underline_index-1);
             CH= str2num (filename(file_index1+2 : underline_index-1));
             if ismember(CH,mpfc_ch_key) == 1
                 IED_all(k).channel = CH;
                 IED_all(k).res = IED;
                 IED_all(k).reg = 'mPFC';
                 IED_all(k).IED_h = length(IED)./NREM_duration.*3600;
                 k = k+1;
             elseif ismember(CH,hc_ch_key) == 1
                 IED_all(k).channel = CH;
                 IED_all(k).res = IED;
                 IED_all(k).reg = 'HC';
                 IED_all(k).IED_h = length(IED)./NREM_duration.*3600;
                 k = k+1;
              elseif ismember(CH,ng1_ch_key) == 1
                 IED_all(k).channel = CH;
                 IED_all(k).res = IED;
                 IED_all(k).reg = 'NG1';
                 IED_all(k).IED_h = length(IED)./NREM_duration.*3600;
                 k = k+1; 
              elseif ismember(CH,ng2_ch_key) == 1
                 IED_all(k).channel = CH;
                 IED_all(k).res = IED;
                 IED_all(k).reg = 'NG2';
                 IED_all(k).IED_h = length(IED)./NREM_duration.*3600;
                 k = k+1; 
              elseif ismember(CH,ng3_ch_key) == 1
                 IED_all(k).channel = CH;
                 IED_all(k).res = IED;
                 IED_all(k).reg = 'NG3';
                 IED_all(k).IED_h = length(IED)./NREM_duration.*3600;
                 k = k+1; 
              elseif ismember(CH,ng4_ch_key) == 1
                 IED_all(k).channel = CH;
                 IED_all(k).res = IED;
                 IED_all(k).reg = 'NG4';
                 IED_all(k).IED_h = length(IED)./NREM_duration.*3600;
                 k = k+1; 
              elseif ismember(CH,ng5_ch_key) == 1
                 IED_all(k).channel = CH;
                 IED_all(k).res = IED;
                 IED_all(k).reg = 'NG5';
                 IED_all(k).IED_h = length(IED)./NREM_duration.*3600;
                 k = k+1; 
              elseif ismember(CH,ng6_ch_key) == 1
                 IED_all(k).channel = CH;
                 IED_all(k).res = IED;
                 IED_all(k).reg = 'NG6';
                 IED_all(k).IED_h = length(IED)./NREM_duration.*3600;
                 k = k+1; 
              elseif ismember(CH,ng7_ch_key) == 1
                 IED_all(k).channel = CH;
                 IED_all(k).res = IED;
                 IED_all(k).reg = 'NG7';
                 IED_all(k).IED_h = length(IED)./NREM_duration.*3600;
                 k = k+1; 
              elseif ismember(CH,ng8_ch_key) == 1
                 IED_all(k).channel = CH;
                 ED_all(k).res = IED;
                 IED_all(k).reg = 'NG8';
                 IED_all(k).IED_h = length(IED)./NREM_duration.*3600;
                 k = k+1;    
              elseif ismember(CH,ng9_ch_key) == 1
                 IED_all(k).channel = CH;
                 IED_all(k).res = IED;
                 IED_all(k).reg = 'NG9';
                 IED_all(k).IED_h = length(IED)./NREM_duration.*3600;
                 k = k+1;    
            end % if
               
               if IED_all(k-1).res == 0;
                  IED_all(k-1).IED_h = 0;
               end
               
%                if isempty (IED_all(k-1).res) == 1;
%                   IED_all(k-1).IED_h = 0;
%                end
     
          end % for
         
          %Sort by region
          IED_all_table = struct2table(IED_all);
          IED_all_sort = sortrows(IED_all_table,'reg');
          IED_all = table2struct(IED_all_sort);
          
          %Save file
          file_index2 = strfind(filename,'post');
          fbasename= strcat(filename(1:file_index2+3),'_IED_struct');
          save(fbasename,'IED_all');
         
          %Copy to a folder in rat directory
          path_sub = pwd;
          file_name = dir('*_struct.mat');
          path_file = strcat (path_sub,'/',file_name.name);
          
          cd ../..;
          path_rat = pwd;
          index_rat = strfind (path_rat,'/');
          index_rat = index_rat(end);
          rat_ID = path_rat (index_rat+1:end);
          path_copy = strcat (path_rat, '/','IED_rate_analysis','_',rat_ID);
          folder_name_rate = strcat ('IED_rate_analysis','_',rat_ID);
          prev_rat_dir = dir ('IED_rate_analysis*');
          if isempty (prev_rat_dir) == 1;
             mkdir (folder_name_rate);
          end
          copyfile(path_file,path_copy);
       
     else %if states
          disp ('NO_states_OR_detect_IED_files');  
          cd ..;
    end % if states
     
 end
 


close all
clear
clc

addpath('jsonlab')
main_path = '/home/lmiguelgato/Documents/DAP_project/corpus';

files = dir();
dirFlags = [files.isdir] & ~strcmp({files.name},'.') & ~strcmp({files.name},'..') & ~strcmp({files.name},'jsonlab') & ~strcmp({files.name},'config');
subFolders = files(dirFlags);

all_data = [];
no_data  = [];
no_10 = 0;
no_5 = 0;

for k = 1 : length(subFolders)
  cd(subFolders(k).name)
  files = dir();
  dirFlags = [files.isdir] & ~strcmp({files.name},'.') & ~strcmp({files.name},'..');
  subFolders1 = files(dirFlags);
  for j = 1 : length(subFolders1)
    cd(subFolders1(j).name)
    files = dir();
    dirFlags = [files.isdir] & ~strcmp({files.name},'.') & ~strcmp({files.name},'..');
    subFolders2 = files(dirFlags);
    
    nsources = str2double(subFolders1(j).name(1));
    str_regexp = '';
    for isources = 1:nsources
        str_regexp = strcat(str_regexp, '%d_');
    end
    str_regexp = strcat(str_regexp, '_');
    for isources = 1:nsources
        str_regexp = strcat(str_regexp, '%f_');
    end
    str_regexp = strcat(str_regexp, '_%u');
    
    for i = 1 : length(subFolders2)
      
      cd(subFolders2(i).name)
      files = dir();
      dirFlags = [files.isdir] & ~strcmp({files.name},'.') & ~strcmp({files.name},'..');
      subFolders3 = files(dirFlags);
      
      
      
      for h = 1 : length(subFolders3)
          thisPath = strcat(subFolders(k).name, '/', subFolders1(j).name, ...
            '/', subFolders2(i).name, '/', subFolders3(h).name);
        
          cd ([main_path '/' thisPath])
          
          mae15 = 0;
          mae10 = 0;
          mae5 = 0;
          rms15 = 0;
          rms10 = 0;
          rms5 = 0;
          e15 = 0;
          e10 = 0;
          e5 = 0;
          some_source = false;
      
          for isources = 1:nsources
            fileID = fopen(['source_' num2str(isources) '.txt'],'r');
            if fileID ~= -1
                some_source = true;
                C = textscan(fileID,'%s',1,'delimiter','\n', 'headerlines',3);
                mae15 = mae15 + sscanf(C{1}{1},'mean absolute error up to 15 degrees in error: %f');
                frewind(fileID)
                C = textscan(fileID,'%s',1,'delimiter','\n', 'headerlines',4);
                rms15 = rms15 + sscanf(C{1}{1},'RMS error up to 15 degrees in error: %f');
                frewind(fileID)
                C = textscan(fileID,'%s',1,'delimiter','\n', 'headerlines',8);
                tmp = sscanf(C{1}{1},'mean absolute error up to 10 degrees in error: %f');
                if tmp == -1.0
                    no_10 = no_10 + 1;
                    tmp = 0.0;
                end
                mae10 = mae10 + tmp;
                frewind(fileID)
                C = textscan(fileID,'%s',1,'delimiter','\n', 'headerlines',9);
                tmp = sscanf(C{1}{1},'RMS error up to 10 degrees in error: %f');
                if tmp == -1.0
                    tmp = 0.0;
                end
                rms10 = rms10 + tmp;
                frewind(fileID)
                C = textscan(fileID,'%s',1,'delimiter','\n', 'headerlines',13);
                tmp = sscanf(C{1}{1},'mean absolute error up to 5 degrees in error: %f');
                if tmp == -1.0
                    no_5 = no_5 + 1;
                    tmp = 0.0;
                end
                mae5 = mae5 + tmp;
                frewind(fileID)
                C = textscan(fileID,'%s',1,'delimiter','\n', 'headerlines',14);
                tmp = sscanf(C{1}{1},'RMS error up to 5 degrees in error: %f');
                if tmp == -1.0
                    tmp = 0.0;
                end
                rms5 = rms5 + tmp;
                frewind(fileID)
                C = textscan(fileID,'%s',1,'delimiter','\n', 'headerlines',16);
                e15 = e15 + sscanf(C{1}{1},'E15 metric: %f');
                frewind(fileID)
                C = textscan(fileID,'%s',1,'delimiter','\n', 'headerlines',17);
                e10 = e10 + sscanf(C{1}{1},'E10 metric: %f');
                frewind(fileID)
                C = textscan(fileID,'%s',1,'delimiter','\n', 'headerlines',18);
                e5 = e5 + sscanf(C{1}{1},'E5 metric: %f');
                fclose(fileID);
            end
          end
          mm = nsources;
          if some_source
            tmp = [k, j, i, mae15/mm, mae10/mm, mae5/mm, rms15/mm, rms10/mm, rms5/mm, e15/mm, e10/mm, e5/mm];
            all_data = [all_data; tmp];
          else
            no_data = [no_data; [k, j, i, h]];
          end
      end
        
      cd (strcat(main_path, '/', subFolders(k).name, '/', subFolders1(j).name))
    end
    cd ('..')
  end
  cd ('..')
end

mm = size(all_data, 1);

mean_mae15 = sum(all_data(:, 4))/mm;
mean_mae10 = sum(all_data(:, 5))/(mm-no_10);
mean_mae5 = sum(all_data(:, 6))/(mm-no_5);
mean_rms15 = sum(all_data(:, 7))/mm;
mean_rms10 = sum(all_data(:, 8))/(mm-no_10);
mean_rms5 = sum(all_data(:, 9))/(mm-no_5);
mean_e15    = mean(all_data(:, 10));
mean_e10    = mean(all_data(:, 11));
mean_e5     = mean(all_data(:, 12));
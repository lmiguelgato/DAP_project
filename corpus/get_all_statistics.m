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
all_missedSources = [];
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
          
          fileID = fopen('errors.txt','r');
          C = textscan(fileID,'%s',1,'delimiter','\n', 'headerlines',1);
          Pfa = sscanf(C{1}{1},'false alarm rate: %f');
          frewind(fileID)
          C = textscan(fileID,'%s',1,'delimiter','\n', 'headerlines',2);
          missedSources = sscanf(C{1}{1},'missed sources: %u');
          all_missedSources = [all_missedSources; [missedSources, nsources, Pfa]];
          fclose(fileID);
          
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
            tmp = [k, j, i, mae15/mm, mae10/mm, mae5/mm, rms15/mm, rms10/mm, rms5/mm, e15/mm, e10/mm, e5/mm, nsources];
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

mean_mae15  = sum(all_data(:, 4))/mm;
mean_mae10  = sum(all_data(:, 5))/(mm-no_10);
mean_mae5   = sum(all_data(:, 6))/(mm-no_5);
mean_rms15  = sum(all_data(:, 7))/mm;
mean_rms10  = sum(all_data(:, 8))/(mm-no_10);
mean_rms5   = sum(all_data(:, 9))/(mm-no_5);
mean_e15    = mean(all_data(:, 10));
mean_e10    = mean(all_data(:, 11));
mean_e5     = mean(all_data(:, 12));

mean_mae15_1s  = mean(all_data(all_data(:, 13) == 1, 4));
mean_rms15_1s  = mean(all_data(all_data(:, 13) == 1, 7));
mean_e15_1s    = mean(all_data(all_data(:, 13) == 1, 10));
mean_e10_1s    = mean(all_data(all_data(:, 13) == 1, 11));
mean_e5_1s     = mean(all_data(all_data(:, 13) == 1, 12));

mean_mae15_2s  = mean(all_data(all_data(:, 13) == 2, 4));
mean_rms15_2s  = mean(all_data(all_data(:, 13) == 2, 7));
mean_e15_2s    = mean(all_data(all_data(:, 13) == 2, 10));
mean_e10_2s    = mean(all_data(all_data(:, 13) == 2, 11));
mean_e5_2s     = mean(all_data(all_data(:, 13) == 2, 12));


mean_mae15_3s  = mean(all_data(all_data(:, 13) == 3, 4));
mean_rms15_3s  = mean(all_data(all_data(:, 13) == 3, 7));
mean_e15_3s    = mean(all_data(all_data(:, 13) == 3, 10));
mean_e10_3s    = mean(all_data(all_data(:, 13) == 3, 11));
mean_e5_3s     = mean(all_data(all_data(:, 13) == 3, 12));

mean_mae15_4s  = mean(all_data(all_data(:, 13) == 4, 4));
mean_rms15_4s  = mean(all_data(all_data(:, 13) == 4, 7));
mean_e15_4s    = mean(all_data(all_data(:, 13) == 4, 10));
mean_e10_4s    = mean(all_data(all_data(:, 13) == 4, 11));
mean_e5_4s     = mean(all_data(all_data(:, 13) == 4, 12));



Pmd = all_missedSources(:, 1)./all_missedSources(:, 2);

N1 = sum(all_missedSources(:, 2) == 1);
N01 = sum(all_missedSources(:, 1) == 0 & all_missedSources(:, 2) == 1);
N11 = sum(all_missedSources(:, 1) == 1 & all_missedSources(:, 2) == 1);
Pfa1 = mean(all_missedSources(all_missedSources(:, 2) == 1, 3));

N2 = sum(all_missedSources(:, 2) == 2);
N02 = sum(all_missedSources(:, 1) == 0 & all_missedSources(:, 2) == 2);
N12 = sum(all_missedSources(:, 1) == 1 & all_missedSources(:, 2) == 2);
N22 = sum(all_missedSources(:, 1) == 2 & all_missedSources(:, 2) == 2);
Pfa2 = mean(all_missedSources(all_missedSources(:, 2) == 2, 3));

N3 = sum(all_missedSources(:, 2) == 3);
N03 = sum(all_missedSources(:, 1) == 0 & all_missedSources(:, 2) == 3);
N13 = sum(all_missedSources(:, 1) == 1 & all_missedSources(:, 2) == 3);
N23 = sum(all_missedSources(:, 1) == 2 & all_missedSources(:, 2) == 3);
N33 = sum(all_missedSources(:, 1) == 3 & all_missedSources(:, 2) == 3);
Pfa3 = mean(all_missedSources(all_missedSources(:, 2) == 3, 3));

N4 = sum(all_missedSources(:, 2) == 4);
N04 = sum(all_missedSources(:, 1) == 0 & all_missedSources(:, 2) == 4);
N14 = sum(all_missedSources(:, 1) == 1 & all_missedSources(:, 2) == 4);
N24 = sum(all_missedSources(:, 1) == 2 & all_missedSources(:, 2) == 4);
N34 = sum(all_missedSources(:, 1) == 3 & all_missedSources(:, 2) == 4);
N44 = sum(all_missedSources(:, 1) == 3 & all_missedSources(:, 2) == 4);
Pfa4 = mean(all_missedSources(all_missedSources(:, 2) == 4, 3));

disp('When 1 source:')
disp(['* ' num2str(N11/N1*100) ' % of times no sources were detected'])
disp(['* ' num2str(N01/N1*100) ' % of times one source was detected'])
disp('***')
disp(['* ' num2str(Pfa1) ' % of detections were false alarms (outside +/- 15 degrees from any source)'])
disp(['* ' num2str(mean_mae15_1s) ' degrees of average mean absolute error'])
disp(['* ' num2str(mean_rms15_1s) ' degrees of average RMS error'])
disp(['* ' num2str(mean_e5_1s) ' % of detections between 0 and 5 error degrees'])
disp(['* ' num2str(mean_e10_1s) ' % of detections between 5 and 10 error degrees'])
disp(['* ' num2str(mean_e15_1s) ' % of detections between 10 and 15 error degrees'])

disp(' ')

disp('When 2 sources:')
disp(['* ' num2str(N22/N2*100) ' % of times no sources were detected'])
disp(['* ' num2str(N12/N2*100) ' % of times one source was detected'])
disp(['* ' num2str(N02/N2*100) ' % of times two sources were detected'])
disp('***')
disp(['* ' num2str(Pfa2) ' % of detections were false alarms (outside +/- 15 degrees from any source)'])
disp(['* ' num2str(mean_mae15_2s) ' degrees of average mean absolute error'])
disp(['* ' num2str(mean_rms15_2s) ' degrees of average RMS error'])
disp(['* ' num2str(mean_e5_2s) ' % of detections between 0 and 5 error degrees'])
disp(['* ' num2str(mean_e10_2s) ' % of detections between 5 and 10 error degrees'])
disp(['* ' num2str(mean_e15_2s) ' % of detections between 10 and 15 error degrees'])
disp(' ')

disp('When 3 sources:')
disp(['* ' num2str(N33/N3*100) ' % of times no sources were detected'])
disp(['* ' num2str(N23/N3*100) ' % of times one source was detected'])
disp(['* ' num2str(N13/N3*100) ' % of times two sources were detected'])
disp(['* ' num2str(N03/N3*100) ' % of times three sources were detected'])
disp('***')
disp(['* ' num2str(Pfa3) ' % of detections were false alarms (outside +/- 15 degrees from any source)'])
disp(['* ' num2str(mean_mae15_3s) ' degrees of average mean absolute error'])
disp(['* ' num2str(mean_rms15_3s) ' degrees of average RMS error'])
disp(['* ' num2str(mean_e5_3s) ' % of detections between 0 and 5 error degrees'])
disp(['* ' num2str(mean_e10_3s) ' % of detections between 5 and 10 error degrees'])
disp(['* ' num2str(mean_e15_3s) ' % of detections between 10 and 15 error degrees'])
disp(' ')

disp('When 4 sources:')
disp(['* ' num2str(N44/N4*100) ' % of times no sources were detected'])
disp(['* ' num2str(N34/N4*100) ' % of times one source was detected'])
disp(['* ' num2str(N24/N4*100) ' % of times two sources were detected'])
disp(['* ' num2str(N14/N4*100) ' % of times three sources were detected'])
disp(['* ' num2str(N04/N4*100) ' % of times four sources were detected'])
disp('***')
disp(['* ' num2str(Pfa4) ' % of detections were false alarms (outside +/- 15 degrees from any source)'])
disp(['* ' num2str(mean_mae15_4s) ' degrees of average mean absolute error'])
disp(['* ' num2str(mean_rms15_4s) ' degrees of average RMS error'])
disp(['* ' num2str(mean_e5_4s) ' % of detections between 0 and 5 error degrees'])
disp(['* ' num2str(mean_e10_4s) ' % of detections between 5 and 10 error degrees'])
disp(['* ' num2str(mean_e15_4s) ' % of detections between 10 and 15 error degrees'])
disp(' ')


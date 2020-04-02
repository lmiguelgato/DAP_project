close all
clear
clc

clean_corpus_txts

Emin_ssl = 0.050;    % for Pfa = 0.10: 0.050  total, 0.049 reverb, 0.050 anech
Emax_ssl = 1.00;

addpath('jsonlab')
main_path = '/home/lmiguelgato/Documents/DAP_project/corpus';

files = dir();
dirFlags = [files.isdir] & ~strcmp({files.name},'.') & ~strcmp({files.name},'..') & ~strcmp({files.name},'jsonlab') & ~strcmp({files.name},'config');
subFolders = files(dirFlags);

max_abs_error_E5 = 5.0;
max_abs_error_E10 = 10.0;
max_abs_error_E15 = 15.0;

for k = 1 : length(subFolders)
%     if strcmp(subFolders(k).name, 'Anechoic_Chamber_16')
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
              E5 = sscanf(subFolders2(i).name, str_regexp);
              theta_original = E5(1:nsources);

              cd(subFolders2(i).name)
              files = dir();
              dirFlags = [files.isdir] & ~strcmp({files.name},'.') & ~strcmp({files.name},'..');
              subFolders3 = files(dirFlags);

              for h = 1 : length(subFolders3)
                  thisPath = strcat(subFolders(k).name, '/', subFolders1(j).name, ...
                    '/', subFolders2(i).name, '/', subFolders3(h).name);

                  cd ([main_path '/' thisPath])

                  load('ssl0-1.mat')

                  false_alarm = ones(size(theta_ssl));

                    total_detections_E5 = 0;
                    total_detections_E10 = 0;
                    total_detections_E15 = 0;
                    n_missed_detections = 0;

                    nanIdx = find(or(E_ssl < Emin_ssl, E_ssl > Emax_ssl));
                    if ~isempty(nanIdx)
                        theta_ssl(nanIdx) = NaN;
                        false_alarm(nanIdx) = 0;
                    end

                  for isources = 1:nsources
                    E5 = find(abs(theta_ssl*180/pi-theta_original(isources)) < max_abs_error_E5);
                    E10 = find(abs(theta_ssl*180/pi-theta_original(isources)) < max_abs_error_E10);
                    E15 = find(abs(theta_ssl*180/pi-theta_original(isources)) < max_abs_error_E15);

                    num_detections_E5 = 0;
                    num_detections_E10 = 0;
                    num_detections_E15 = 0;

                    if ~isempty(E15)
                        fileID = fopen(['source_' num2str(isources) '.txt'],'w');
                        fprintf(fileID,'best match: %6.2f\n', theta_original(isources));
                        fprintf(fileID,'number of detections up to %u degrees in error: %u\n', max_abs_error_E15, length(E15));
                        fprintf(fileID,'mean detection up to %u degrees in error: %6.2f\n', max_abs_error_E15, mean(theta_ssl(E15)*180/pi));
                        fprintf(fileID,'mean absolute error up to %u degrees in error: %6.2f\n', max_abs_error_E15, mean(abs(theta_original(isources)-theta_ssl(E15)*180/pi)));
                        fprintf(fileID,'RMS error up to %u degrees in error: %6.2f\n', max_abs_error_E15, sqrt(sum((theta_original(isources)-theta_ssl(E15)*180/pi).^2)/length(E15)));
                        fprintf(fileID,'------------------------\n');
                        num_detections_E15 = num_detections_E15 + length(E15);
                        false_alarm(E15) = 0;


                        if ~isempty(E10)
                            num_detections_E10 = num_detections_E10 + length(E10);
                            num_detections_E15 = num_detections_E15 - num_detections_E10;
                            fprintf(fileID,'number of detections up to %u degrees in error: %u\n', max_abs_error_E10, length(E10));
                            fprintf(fileID,'mean detection up to %u degrees in error: %6.2f\n', max_abs_error_E10, mean(theta_ssl(E10)*180/pi));
                            fprintf(fileID,'mean absolute error up to %u degrees in error: %6.2f\n', max_abs_error_E10, mean(abs(theta_original(isources)-theta_ssl(E10)*180/pi)));
                            fprintf(fileID,'RMS error up to %u degrees in error: %6.2f\n', max_abs_error_E10, sqrt(sum((theta_original(isources)-theta_ssl(E10)*180/pi).^2)/length(E10)));
                            fprintf(fileID,'------------------------\n');

                            if ~isempty(E5)
                                num_detections_E5 = num_detections_E5 + length(E5);
                                num_detections_E10 = num_detections_E10 - num_detections_E5;                

                                fprintf(fileID,'number of detections up to %u degrees in error: %u\n', max_abs_error_E5, length(E5));
                                fprintf(fileID,'mean detection up to %u degrees in error: %6.2f\n', max_abs_error_E5, mean(theta_ssl(E5)*180/pi));
                                fprintf(fileID,'mean absolute error up to %u degrees in error: %6.2f\n', max_abs_error_E5, mean(abs(theta_original(isources)-theta_ssl(E5)*180/pi)));
                                fprintf(fileID,'RMS error up to %u degrees in error: %6.2f\n', max_abs_error_E5, sqrt(sum((theta_original(isources)-theta_ssl(E5)*180/pi).^2)/length(E5)));

                                fprintf(fileID,'------------------------\n');
                            else
                                fprintf(fileID,'number of detections up to %u degrees in error: %u\n', max_abs_error_E5, 0);
                                fprintf(fileID,'mean detection up to %u degrees in error: %6.2f\n', max_abs_error_E5, 0);
                                fprintf(fileID,'mean absolute error up to %u degrees in error: %6.2f\n', max_abs_error_E5, -1);
                                fprintf(fileID,'RMS error up to %u degrees in error: %6.2f\n', max_abs_error_E5, -1);
                                fprintf(fileID,'------------------------\n');
                            end
                        else
                            fprintf(fileID,'number of detections up to %u degrees in error: %u\n', max_abs_error_E10, 0);
                            fprintf(fileID,'mean detection up to %u degrees in error: %6.2f\n', max_abs_error_E10, 0);
                            fprintf(fileID,'mean absolute error up to %u degrees in error: %6.2f\n', max_abs_error_E10, -1);
                            fprintf(fileID,'RMS error up to %u degrees in error: %6.2f\n', max_abs_error_E10, -1);
                            fprintf(fileID,'------------------------\n');
                            fprintf(fileID,'number of detections up to %u degrees in error: %u\n', max_abs_error_E5, 0);
                            fprintf(fileID,'mean detection up to %u degrees in error: %6.2f\n', max_abs_error_E5, 0);
                            fprintf(fileID,'mean absolute error up to %u degrees in error: %6.2f\n', max_abs_error_E5, -1);
                            fprintf(fileID,'RMS error up to %u degrees in error: %6.2f\n', max_abs_error_E5, -1);
                            fprintf(fileID,'------------------------\n');
                        end

                        fprintf(fileID,'E15 metric: %6.2f %%\n', num_detections_E15/(num_detections_E15 + num_detections_E10 + num_detections_E5)*100);
                        fprintf(fileID,'E10 metric: %6.2f %%\n', num_detections_E10/(num_detections_E15 + num_detections_E10 + num_detections_E5)*100);
                        fprintf(fileID,'E5 metric: %6.2f %%\n', num_detections_E5/(num_detections_E15 + num_detections_E10 + num_detections_E5)*100);
                        fprintf(fileID,'------------------------\n');
                        fprintf(fileID,'%u  %6.2f\r\n', [E15 theta_ssl(E15)*180/pi]');
                        fclose(fileID);  

                        total_detections_E5 = total_detections_E5 + num_detections_E5;
                        total_detections_E10 = total_detections_E10 + num_detections_E10;
                        total_detections_E15 = total_detections_E15 + num_detections_E15;
                    else
                        n_missed_detections = n_missed_detections+1;
                    end
                  end

                  false_alarm(isnan(theta_ssl)) = 0;
                  false_alarm = sum(sum(false_alarm));

                  fileID = fopen('errors.txt','w');
                  fprintf(fileID,'number of false alarms: %u\r\n', false_alarm);
                  if (false_alarm+total_detections_E15+total_detections_E10+total_detections_E5 ~= 0)
                    fprintf(fileID,'false alarm rate: %6.2f %%\r\n', abs(false_alarm/(false_alarm+total_detections_E5+total_detections_E10+total_detections_E15)*100));
                  else
                    fprintf(fileID,'false alarm rate: -1\r\n');
                  end
                  fprintf(fileID,'missed sources: %u out of %u\n', n_missed_detections, nsources);
                  fclose(fileID);
              end

              cd (strcat(main_path, '/', subFolders(k).name, '/', subFolders1(j).name))
            end
            cd ('..')
        end
        cd ('..')
%     end
end
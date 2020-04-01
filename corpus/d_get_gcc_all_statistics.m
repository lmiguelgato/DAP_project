close all
clear
clc

addpath('jsonlab')
main_path = '/home/lmiguelgato/Documents/DAP_project/corpus';
enableML = 0;
max_num_sources = 4;
enable_plots = 0;

files = dir();
dirFlags = [files.isdir] & ~strcmp({files.name},'.') & ~strcmp({files.name},'..') & ~strcmp({files.name},'jsonlab') & ~strcmp({files.name},'config');
subFolders = files(dirFlags);

all_data = [];
no_data  = [];
all_missedSources = [];
no_10 = 0;
no_5 = 0;

Emin_ssl = 0.10;
Emax_ssl = 1.00;
max_abs_error_E5 = 5.0;
max_abs_error_E10 = 10.0;
max_abs_error_E15 = 15.0;

for k = 1 : length(subFolders)
    if ~strcmp(subFolders(k).name, 'Anechoic_Chamber_16')
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

              clc
              disp(['Processing: ' thisPath])        

              %for isources = 1:nsources
                fileID = fopen('gcc_ssl.txt','r');
                if fileID ~= -1

                    format = '%d:';
                    for ii = 1:max_num_sources
                        format = strcat(format, ' %f');
                    end

                    A = textscan(fileID, format, 'Delimiter',',','EmptyValue',NaN);
                    fclose(fileID);

                    false_alarm = ones(length(A{1}), length(A)-1);

                    theta_ssl = [A{2}, A{3}, A{4}, A{5}];

                    theta_ssl(theta_ssl == 181) = NaN;

                    total_detections_E5     = 0;
                    total_detections_E10    = 0;
                    total_detections_E15    = 0;
                    n_missed_detections     = 0;

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
                        E5 = find(abs(theta_ssl-theta_original(isources)) < max_abs_error_E5);
                        E10 = find(abs(theta_ssl-theta_original(isources)) < max_abs_error_E10);
                        E15 = find(abs(theta_ssl-theta_original(isources)) < max_abs_error_E15);

                        num_detections_E5 = 0;
                        num_detections_E10 = 0;
                        num_detections_E15 = 0;

                        if ~isempty(E15)
                            mae15 = mae15 + mean(abs(theta_original(isources)-theta_ssl(E15)));
                            rms15 = rms15 + sqrt(sum((theta_original(isources)-theta_ssl(E15)).^2)/length(E15));
                            num_detections_E15 = num_detections_E15 + length(E15);
                            false_alarm(E15) = 0;

                            if ~isempty(E10)
                                num_detections_E10 = num_detections_E10 + length(E10);
                                num_detections_E15 = num_detections_E15 - num_detections_E10;
                                mae10 = mae10 + mean(abs(theta_original(isources)-theta_ssl(E10)));
                                rms10 = rms10 + sqrt(sum((theta_original(isources)-theta_ssl(E10)).^2)/length(E10));

                                if ~isempty(E5)
                                    num_detections_E5 = num_detections_E5 + length(E5);
                                    num_detections_E10 = num_detections_E10 - num_detections_E5;
                                    mae5 = mae5 + mean(abs(theta_original(isources)-theta_ssl(E5)));
                                    rms5 = rms5 + sqrt(sum((theta_original(isources)-theta_ssl(E5)).^2)/length(E5));
                                end
                            end

                            e15 = e15 + num_detections_E15/(num_detections_E15 + num_detections_E10 + num_detections_E5)*100;
                            e10 = e10 + num_detections_E10/(num_detections_E15 + num_detections_E10 + num_detections_E5)*100;
                            e5 = e5 + num_detections_E5/(num_detections_E15 + num_detections_E10 + num_detections_E5)*100;

                            total_detections_E5 = total_detections_E5 + num_detections_E5;
                            total_detections_E10 = total_detections_E10 + num_detections_E10;
                            total_detections_E15 = total_detections_E15 + num_detections_E15;
                        else
                            n_missed_detections = n_missed_detections+1;
                        end
                    end

                    if sum(sum(~isnan(theta_ssl))) ~= 0
                        some_source = true;
                    end

                    mm = nsources;
                    if some_source
                    tmp = [k, j, i, mae15/mm, mae10/mm, mae5/mm, rms15/mm, rms10/mm, rms5/mm, e15/mm, e10/mm, e5/mm, nsources];
                    all_data = [all_data; tmp];
                    else
                    no_data = [no_data; [k, j, i, h]];
                    end

                    false_alarm(isnan(theta_ssl)) = 0;
                    false_alarm = sum(sum(false_alarm));

                    if (false_alarm+total_detections_E15+total_detections_E10+total_detections_E5 ~= 0)
                        Pfa = abs(false_alarm/(false_alarm+total_detections_E5+total_detections_E10+total_detections_E15)*100);
                    else
                        Pfa = 0;
                    end

                    all_missedSources = [all_missedSources; [n_missed_detections, nsources, Pfa]];

                    fileID  = fopen('gcc_sst.txt','r');
                    format = '%d %f %f';
                    K = textscan(fileID, format,'EmptyValue',NaN);
                    fclose(fileID);

                    instantaneousSource = K{1};
                    instantaneousDOA = K{2};
                    kalmanDOA = K{3};

                    numDOAs   = zeros(max_num_sources,1);

                    for ii = 1:max_num_sources
                        temp = A{1+ii};
                        for jj = 1:length(temp)
                            if temp(jj) == 181
                                temp(jj) = NaN;
                            else
                                numDOAs(ii) = numDOAs(ii) + 1;
                            end
                        end
                        A{1+ii} = temp;
                    end

                    numKalmanDOAs   = zeros(max_num_sources,1);

                    temp = K{1};
                    for ii = 1:max_num_sources
                        for jj = 1:length(temp)
                            if temp(jj) == ii-1
                                numKalmanDOAs(ii) = numKalmanDOAs(ii) + 1;
                            end
                        end
                    end

                    meanDOAs  = zeros(max_num_sources,1);
                    kalmanDOAs  = zeros(max_num_sources,1);
                    stdevDOAs = meanDOAs;
                    kalmanStdevDOAs = kalmanDOAs;
                    percentDOAs   = numDOAs/sum(numDOAs)*100;
                    percentKalmanDOAs   = numKalmanDOAs/sum(numKalmanDOAs)*100;

                    for ii = 1:max_num_sources
                        idx = find(instantaneousSource == ii-1);
                        meanDOAs(ii)  = nanmean(A{ii+1});
                        kalmanDOAs(ii) = nanmean(kalmanDOA(idx));
                        stdevDOAs(ii) = sqrt(nanvar(A{ii+1}));
                        kalmanStdevDOAs(ii) = sqrt(nanvar(kalmanDOA(idx)));
                    end

                    firstDetection = sum(numDOAs);
                    for ii = 1:max_num_sources
                        tmp = find(~isnan(A{1+ii}), 1);
                        if tmp < firstDetection
                            firstDetection = tmp;
                        end
                    end

                    if enable_plots
                        if (~enableML) 
                            figure
                            axis([1 sum(numDOAs)+firstDetection-1 -180.5 180.5])
                            grid on

                            labels = {};

                            for ii = 1:max_num_sources    
                                idx = find(instantaneousSource == ii-1);
                                hold on; plot(instantaneousDOA(idx),'LineWidth',2)
                                labels{end+1} = strcat('Instantaneous_', num2str(ii));
                                hold on; plot(kalmanDOA(idx),'LineWidth',2)
                                labels{end+1} = strcat('Kalman_', num2str(ii));
                            end
                            legend(labels)
                        end

                        %% ---------- kalman polar
                        if (~enableML) 
                            figure
                            for jj = 1:max_num_sources
                                if ~isnan(kalmanDOAs(jj)) && ~isnan(stdevDOAs(jj))
                                    tmp = exp(1i*[kalmanDOAs(jj) - kalmanStdevDOAs(jj) kalmanDOAs(jj) kalmanDOAs(jj) + kalmanStdevDOAs(jj)]/180*pi);
                                    polarplot(tmp,':','LineWidth',ceil(8*percentKalmanDOAs(jj)/100));
                                    hold on;
                                    text(imag(tmp(2)),real(tmp(2)),[num2str(round(percentKalmanDOAs(jj))) ' %'])
                                end
                            end
                            polarplot(0.1*exp(1i*[-60 60 180 -60]/180*pi),'k','LineWidth',2);

                            for jj = 1:max_num_sources
                                hold on;
                                polarplot([0.9*exp(1i*(kalmanDOAs(jj) - kalmanStdevDOAs(jj))/180*pi) 1.1*exp(1i*(kalmanDOAs(jj) - kalmanStdevDOAs(jj))/180*pi)],'k');
                                hold on;
                                polarplot([0.9*exp(1i*(kalmanDOAs(jj) + kalmanStdevDOAs(jj))/180*pi) 1.1*exp(1i*(kalmanDOAs(jj) + kalmanStdevDOAs(jj))/180*pi)],'k');
                            end

                            ax = gca;
                            ax.ThetaLim = [-180 180];
                            ax.RTickLabel = {''};
                            ax.ThetaZeroLocation = 'top';
                            ax.ThetaDir = 'clockwise';
                        end
                    end

                    numDetectedSources = sum(percentKalmanDOAs > max(percentKalmanDOAs/2));

                    if (~enableML) 
                        disp(['Total number of sources detected: ' num2str(sum(percentKalmanDOAs>0)) ', found at: '])
                        disp(num2str(kalmanDOAs(find(percentKalmanDOAs))))
                        disp(' ')
                    end
                    disp(['Maximum-likelihood number of soures: ' num2str(numDetectedSources) ', found at: '])
                    disp(num2str(kalmanDOAs(find(percentKalmanDOAs > max(percentKalmanDOAs)/2))))

                    %% ---------- ML tracking

                    if enable_plots
                        labels = {};
                        figure
                        MLindexes = find(percentKalmanDOAs > max(percentKalmanDOAs)/2);
                        for ii = 1:length(MLindexes)
                            idx = find(instantaneousSource == MLindexes(ii)-1);
                            hold on; plot(instantaneousDOA(idx),'LineWidth',2)
                            labels{end+1} = strcat('Instantaneous_', num2str(ii));
                            hold on; plot(kalmanDOA(idx),'LineWidth',2)
                            labels{end+1} = strcat('Kalman_', num2str(ii));
                        end
                        legend(labels)

                        axis([1 sum(numDOAs)+firstDetection-1 -180.5 180.5])
                        grid on
                    end

                    %% ---------- ML kalman polar

                    if enable_plots
                        figure
                        for idx = 1:length(MLindexes)
                            jj = MLindexes(idx);
                            if ~isnan(kalmanDOAs(jj)) && ~isnan(stdevDOAs(jj))
                                tmp = exp(1i*[kalmanDOAs(jj) - kalmanStdevDOAs(jj) kalmanDOAs(jj) kalmanDOAs(jj) + kalmanStdevDOAs(jj)]/180*pi);
                                polarplot(tmp,'LineWidth',ceil(8*percentKalmanDOAs(jj)/100));
                                hold on;
                                text(imag(tmp(2)),real(tmp(2)),[num2str(round(percentKalmanDOAs(jj))) ' %'])
                            end
                        end
                        polarplot(0.1*exp(1i*[-60 60 180 -60]/180*pi),'k','LineWidth',2);

                        for idx = 1:length(MLindexes)
                            jj = MLindexes(idx);
                            hold on;
                            polarplot([0.9*exp(1i*(kalmanDOAs(jj) - kalmanStdevDOAs(jj))/180*pi) 1.1*exp(1i*(kalmanDOAs(jj) - kalmanStdevDOAs(jj))/180*pi)],'k');
                            hold on;
                            polarplot([0.9*exp(1i*(kalmanDOAs(jj) + kalmanStdevDOAs(jj))/180*pi) 1.1*exp(1i*(kalmanDOAs(jj) + kalmanStdevDOAs(jj))/180*pi)],'k');
                        end

                        ax = gca;
                        ax.ThetaLim = [-180 180];
                        ax.RTickLabel = {''};
                        ax.ThetaZeroLocation = 'top';
                        ax.ThetaDir = 'clockwise';
                    end
                end
          end

          cd (strcat(main_path, '/', subFolders(k).name, '/', subFolders1(j).name))
        end
        cd ('..')
      end
      cd ('..')
    end
end

clc

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
    mean_e15_1s    = 100*sum(all_data(all_data(:, 13) == 1, 10))/(sum(all_data(all_data(:, 13) == 1, 10))+sum(all_data(all_data(:, 13) == 1, 11))+sum(all_data(all_data(:, 13) == 1, 12)));
    mean_e10_1s    = 100*sum(all_data(all_data(:, 13) == 1, 11))/(sum(all_data(all_data(:, 13) == 1, 10))+sum(all_data(all_data(:, 13) == 1, 11))+sum(all_data(all_data(:, 13) == 1, 12)));
    mean_e5_1s     = 100*sum(all_data(all_data(:, 13) == 1, 12))/(sum(all_data(all_data(:, 13) == 1, 10))+sum(all_data(all_data(:, 13) == 1, 11))+sum(all_data(all_data(:, 13) == 1, 12)));
%     mean_e15_1s    = mean(all_data(all_data(:, 13) == 1, 10));
%     mean_e10_1s    = mean(all_data(all_data(:, 13) == 1, 11));
%     mean_e5_1s     = mean(all_data(all_data(:, 13) == 1, 12));

    mean_mae15_2s  = mean(all_data(all_data(:, 13) == 2, 4));
    mean_rms15_2s  = mean(all_data(all_data(:, 13) == 2, 7));
    mean_e15_2s    = 100*sum(all_data(all_data(:, 13) == 2, 10))/(sum(all_data(all_data(:, 13) == 2, 10))+sum(all_data(all_data(:, 13) == 2, 11))+sum(all_data(all_data(:, 13) == 2, 12)));
    mean_e10_2s    = 100*sum(all_data(all_data(:, 13) == 2, 11))/(sum(all_data(all_data(:, 13) == 2, 10))+sum(all_data(all_data(:, 13) == 2, 11))+sum(all_data(all_data(:, 13) == 2, 12)));
    mean_e5_2s     = 100*sum(all_data(all_data(:, 13) == 2, 12))/(sum(all_data(all_data(:, 13) == 2, 10))+sum(all_data(all_data(:, 13) == 2, 11))+sum(all_data(all_data(:, 13) == 2, 12)));


    mean_mae15_3s  = mean(all_data(all_data(:, 13) == 3, 4));
    mean_rms15_3s  = mean(all_data(all_data(:, 13) == 3, 7));
    mean_e15_3s    = 100*sum(all_data(all_data(:, 13) == 3, 10))/(sum(all_data(all_data(:, 13) == 3, 10))+sum(all_data(all_data(:, 13) == 3, 11))+sum(all_data(all_data(:, 13) == 3, 12)));
    mean_e10_3s    = 100*sum(all_data(all_data(:, 13) == 3, 11))/(sum(all_data(all_data(:, 13) == 3, 10))+sum(all_data(all_data(:, 13) == 3, 11))+sum(all_data(all_data(:, 13) == 3, 12)));
    mean_e5_3s     = 100*sum(all_data(all_data(:, 13) == 3, 12))/(sum(all_data(all_data(:, 13) == 3, 10))+sum(all_data(all_data(:, 13) == 3, 11))+sum(all_data(all_data(:, 13) == 3, 12)));

    mean_mae15_4s  = mean(all_data(all_data(:, 13) == 4, 4));
    mean_rms15_4s  = mean(all_data(all_data(:, 13) == 4, 7));
    mean_e15_4s    = 100*sum(all_data(all_data(:, 13) == 4, 10))/(sum(all_data(all_data(:, 13) == 4, 10))+sum(all_data(all_data(:, 13) == 4, 11))+sum(all_data(all_data(:, 13) == 4, 12)));
    mean_e10_4s    = 100*sum(all_data(all_data(:, 13) == 4, 11))/(sum(all_data(all_data(:, 13) == 4, 10))+sum(all_data(all_data(:, 13) == 4, 11))+sum(all_data(all_data(:, 13) == 4, 12)));
    mean_e5_4s     = 100*sum(all_data(all_data(:, 13) == 4, 12))/(sum(all_data(all_data(:, 13) == 4, 10))+sum(all_data(all_data(:, 13) == 4, 11))+sum(all_data(all_data(:, 13) == 4, 12)));



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
N44 = sum(all_missedSources(:, 1) == 4 & all_missedSources(:, 2) == 4);
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

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
    
    for i = 1 : length(subFolders2)      
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
                for i = 1:max_num_sources
                    format = strcat(format, ' %f');
                end

                A = textscan(fileID, format, 'Delimiter',',','EmptyValue',NaN);
                fclose(fileID);

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

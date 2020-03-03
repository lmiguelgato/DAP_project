close all
clear
clc

addpath('jsonlab')
main_path = '/home/lmiguelgato/Documents/DAP_project/corpus';

files = dir();
dirFlags = [files.isdir] & ~strcmp({files.name},'.') & ~strcmp({files.name},'..') & ~strcmp({files.name},'jsonlab') & ~strcmp({files.name},'config');
subFolders = files(dirFlags);

progress = 0;

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
    for i = 1 : length(subFolders2)
      cd(subFolders2(i).name)
      files = dir();
      dirFlags = [files.isdir] & ~strcmp({files.name},'.') & ~strcmp({files.name},'..');
      subFolders3 = files(dirFlags);
      cd (main_path)
      for h = 1 : length(subFolders3)
        thisPath = strcat(subFolders(k).name, '/', subFolders1(j).name, ...
            '/', subFolders2(i).name, '/', subFolders3(h).name);
        
        Emin_ssl = 0.10;
        Emax_ssl = 1.00;
        
        if ~isfile([main_path '/' thisPath '/ssl' num2str(Emin_ssl) ...
                '-' num2str(Emax_ssl) '.mat'])
        
            clc
            disp(['Processing: ' thisPath])

            ssl_data = loadjson([main_path '/' thisPath '/ssl.json']);

            n_detections_ssl    = length(ssl_data);
            n_sources_ssl       = length(ssl_data{1}.src);

            x_ssl       = zeros(n_detections_ssl, n_sources_ssl);
            y_ssl       = zeros(n_detections_ssl, n_sources_ssl);
            theta_ssl   = zeros(n_detections_ssl, n_sources_ssl)+NaN;

            for i_detection = 1:n_detections_ssl
                for i_source = 1:n_sources_ssl
                    this_E = ssl_data{i_detection}.src{i_source}.E;

                    if this_E > Emin_ssl && this_E < Emax_ssl
                        this_x = ssl_data{i_detection}.src{i_source}.x;
                        this_y = ssl_data{i_detection}.src{i_source}.y;

                        if this_x ~= 0
                            this_theta = atan(this_y/this_x);
                        else
                            if this_y >= 0
                                this_theta = pi/2;
                            else
                                this_theta = -pi/2;
                            end
                        end

                        x_ssl(i_detection, i_source)        = this_x;
                        y_ssl(i_detection, i_source)        = this_y;
                        theta_ssl(i_detection, i_source)    = this_theta;
                    end
                end
            end

            save([main_path '/' thisPath '/ssl' num2str(Emin_ssl) ...
                '-' num2str(Emax_ssl) '.mat'], 'theta_ssl')
            
            if ~isfile([main_path '/' thisPath '/sst.mat'])
            
                sst_data = loadjson([main_path '/' thisPath '/sst.json']);

                n_detections_sst    = length(sst_data);
                n_sources_sst       = length(sst_data{1}.src);

                x_sst       = zeros(n_detections_sst, n_sources_sst);
                y_sst       = zeros(n_detections_sst, n_sources_sst);
                id_sst      = zeros(n_detections_sst, n_sources_sst)-NaN;
                theta_sst   = zeros(n_detections_sst, n_sources_sst)-NaN;

                for i_detection = 1:n_detections_sst
                    for i_source = 1:n_sources_sst
                        this_id = sst_data{i_detection}.src{i_source}.id;

                        if this_id ~= 0
                            this_x  = sst_data{i_detection}.src{i_source}.x;
                            this_y  = sst_data{i_detection}.src{i_source}.y;

                            if this_x ~= 0
                                this_theta = atan(this_y/this_x);
                            else
                                if this_y >= 0
                                    this_theta = pi/2;
                                else
                                    this_theta = -pi/2;
                                end
                            end

                            x_sst(i_detection, i_source)        = this_x;
                            y_sst(i_detection, i_source)        = this_y;
                            id_sst(i_detection, i_source)       = this_id;
                            theta_sst(i_detection, i_source)    = this_theta;
                        end
                    end
                end

                id_min = min(min(id_sst));
                id_max = max(max(id_sst));

                if isnan(id_min)
                    disp('No tracks found.')
                else
                    sst_theta = zeros(id_max-id_min+1, n_detections_sst);
                    for id = id_min:id_max
                        idx = find(id_sst == id);
                        tmp_theta = zeros(n_detections_sst,n_sources_sst)+181;
                        tmp_theta(idx) = theta_sst(idx);
                        tmp_theta = sum(tmp_theta, 2) - 3*181;
                        tmp_theta(tmp_theta == 181) = NaN;
                        sst_theta(id-id_min+1, :) = tmp_theta;
                    end
                    save([main_path '/' thisPath '/sst.mat'], 'sst_theta')
                end
            end
            
        end
        
      end
      cd (strcat(main_path, '/', subFolders(k).name, '/', subFolders1(j).name))
    end
    cd ('..')
  end
  cd ('..')
end
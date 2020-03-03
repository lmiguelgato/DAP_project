close all
clear
clc

main_path = '/home/lmiguelgato/Documents/DAP_project/corpus';

files = dir();
dirFlags = [files.isdir] & ~strcmp({files.name},'.') & ~strcmp({files.name},'..');
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
        
        config1   = fileread('config/config1.txt');
        config2   = fileread('config/config2.txt');
        config3   = fileread('config/config3.txt');
        config4   = fileread('config/config4.txt');
        
        fid = fopen([main_path '/' thisPath '/fileConfig.cfg'],'wt');
        fprintf(fid, config1);
        fprintf(fid, [main_path '/' thisPath '/']);
        fprintf(fid, config2);
        fprintf(fid, [main_path '/' thisPath '/']);
        fprintf(fid, config3);
        fprintf(fid, [main_path '/' thisPath '/']);
        fprintf(fid, config4);
        fclose(fid);
        
        clc
        %progress = progress+1;
        %disp([num2str(round(progress/(length(subFolders)*length(subFolders1)*length(subFolders2)*length(subFolders3))*100)) ' %'])
        disp(['Processing: ' thisPath])
        
        out_system = system(['./odaslive -c ' main_path '/' thisPath '/fileConfig.cfg']);
      end
      cd (strcat(main_path, '/', subFolders(k).name, '/', subFolders1(j).name))
    end
    cd ('..')
  end
  cd ('..')
end
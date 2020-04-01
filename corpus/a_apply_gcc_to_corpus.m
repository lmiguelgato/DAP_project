close all
clear
clc

main_path = '/home/lmiguelgato/Documents/DAP_project/corpus';

files = dir();
dirFlags = [files.isdir] & ~strcmp({files.name},'.') & ~strcmp({files.name},'..');
subFolders = files(dirFlags);

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
        for i = 1 : length(subFolders2)
          cd(subFolders2(i).name)
          files = dir();
          dirFlags = [files.isdir] & ~strcmp({files.name},'.') & ~strcmp({files.name},'..');
          subFolders3 = files(dirFlags);
          cd (main_path)
          for h = 1 : length(subFolders3)
            thisPath = strcat(subFolders(k).name, '/', subFolders1(j).name, ...
                '/', subFolders2(i).name, '/', subFolders3(h).name);

            clc
            disp(['Processing: ' thisPath])        
            out_system = system(['./gcc_beamformer_offline ' main_path '/' thisPath '/3mics.wav ' main_path '/' thisPath '/']);
          end
          cd (strcat(main_path, '/', subFolders(k).name, '/', subFolders1(j).name))
        end
        cd ('..')
      end
      cd ('..')
    end
end
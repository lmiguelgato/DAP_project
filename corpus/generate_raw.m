function generate_raw (single_wavs_path, N)
%     
  root_wav_name = 'wav_mic';
  multichannel_wav_path = [single_wavs_path '/'];
  
  if exist(multichannel_wav_path, 'dir') ~= 7
    mkdir(multichannel_wav_path);
  end
  
  [tmp, fs] = audioread([single_wavs_path '/' root_wav_name '1.wav']);
  [m, n] = size(tmp);
  
  if n ~= 1
    error('Wrong audio file dimensions. Aborting ...')
  end

  y = zeros(m, N);
  for i = 1:N
    tmp = audioread([single_wavs_path '/' root_wav_name num2str(i) '.wav']);
    if length(tmp) ~= m
      disp('warning: dimensions mismatch')
    end
    y(:, i) = tmp(1:m);
  end

  audiowrite([multichannel_wav_path num2str(N) 'mics.wav'], y, fs)

  out_system = system(['./wav2raw ' multichannel_wav_path num2str(N) 'mics.wav']);
 
  if out_system ~= 0
    error('Unable to run bash. Aborting ...')
  end
end
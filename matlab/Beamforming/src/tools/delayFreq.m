function y = delayFreq(x, t, fs)

    % Apply time-delay using frequency domain processing

    % INPUT
    % x     - signal to be delayed
    % t     - time delay
    % fs    - sampling rate

    % OUTPUT
    % y     - delayed signal

    N = length(x);                          % signal duration (in samples)

    Xf = fft(x);                            % Fourier domain

    w = [0 1:N/2 (-N/2+1):-1]/N*fs;         % frequency vector (wrapped)

    Yf = zeros(1,N);
    for i = 1:N
        freqDelay = exp(-2j*pi*w(i)*t);     % steering vector for this frequency
        Yf(i) = Xf(i)*freqDelay;            % apply delay in frequency domain
    end

%     y = ifft(Yf);
    y = real(ifft(Yf));

end
function y = delayTime( x, t, fs )

    % Apply time-delay using time domain processing

    % INPUT
    % x     - signal to be delayed
    % t     - time delay
    % fs    - sampling rate

    % OUTPUT
    % y     - delayed signal

    y = zeros(size(x));
    
    nDelay = round(t*fs);
    
    if nDelay >= 0
        y(nDelay+1:end) = x(1:end-nDelay);
    else
        y(1:end+nDelay) = x(-nDelay+1:end);
    end
    

end


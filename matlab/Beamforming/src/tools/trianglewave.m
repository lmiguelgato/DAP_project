function triwave = trianglewave(freq, samplerate)
 
x = [1:samplerate];
wavelength = samplerate/freq;
saw = 2*mod(x, wavelength)/wavelength-1; %start with sawtooth wave
triwave = 2*abs(saw)-1;
 
end
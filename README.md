# Robot audition system
## Location and tracking of multiple speakers using an array of 3 microphones
- Estimate the directions of arrival (DOA) of multiple speakers from a set of recordings.
- Spatial filtering to improve the signal-to-interference ratio.

Requirements:
- libsndfile (version 1.0.28 or higher)
- fftw (version 3.3.8 or higher)
- Eigen (version 3.3.7 or higher)

How to use:

1	- For compilation, run: 'make'.

2(a)- If the directions of arrival (azimuth) are known, and just to perform beamforming, run: 'beamformer'; or ...

2(b)- if the directions of arrival (azimuth) are unkown, first perform azimuth estimation and then beamforing by running: 'gcc_beamformer'.

	$ ./gcc_beamformer <separation between microphones (meters)> <maximum number of sources to localize> <which source to filter (set to 0 to filter all)>

	For example: $ ./gcc_beamformer 0.23 2 1

3(a)- Use JACK to manually connect the audio sources to the corresponfing ports of 'jack_doa_beamformer' JACK client; or ...

3(b)- if the audio sources are .wav files, run: 'ReadMicWavs'

	$ ./ReadMicWavs <jack agent's name> <audio file root name> <audio file path> <number of channels>

	For example: $ ./ReadMicWavs jack_doa_beamformer wav_mic corpus/ 3
	
All output files are stored in './output'. MATLAB scripts and functions in './matlab' are meant for analysis of the output files.

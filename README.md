Final project of "Digital Audio Processing" (robotic audition)

- Array signal processing.
- Estimate the directions of arrival (DOA) of multiple speakers from a set of recordings.
- Spatial filtering to improve the signal-to-interference ratio.

How to use:

1	- For compilation, run: 'make'.

2(a)- If the directions of arrival (azimuth) are known, and just to perform beamforming, run: 'beamformer'; or ...

2(b)- if the directions of arrival (azimuth) are unkown, first perform azimuth estimation and then beamforing by running: 'gcc_beamformer'.

3(a)- Use JACK to manually connect the audio sources to the corresponfing ports of 'jack_doa_beamformer' JACK client; or ...

3(b)- if the audio sources are .wav files, run: 'ReadMicWavs'

All output files are stored in './output'. MATLAB scripts and functions in './matlab' are meant for analysis of the output files.
Technical information regarding the localization and separation algorithms implemented can be found in './doc'.
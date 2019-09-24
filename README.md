Final project of "Digital Audio Processing" (robotic audition)

- Array signal processing.
- Estimate the directions of arrival (DOA) of multiple speakers from a set of recordings.
- Spatial filtering to improve the signal-to-interference ratio.

How to compile:
	$ make

How to use:
1.- Run the JACK agent for source localization and separation:
	$ ./gcc_beamformer <separation between microphones (meters)> <maximum number of sources to localize> <which source to filter (set to 0 to filter all)>
2.- Run the JACK agent reading audio files:
	$ ./ReadMicWavs <jack agent's name> <audio file root name> <audio file path> <number of channels>
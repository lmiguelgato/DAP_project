Baudline by SigBlips for Unix
------------------------------------------------------------------------------
Baudline is a time-frequency browser designed for scientific visualization of
the spectral domain.  Signal analysis is performed by Fourier, correlation,
transfer function, impulse response, and raster transforms that create colorful
spectrograms with vibrant detail.  Conduct test and measurement experiments
with the built-in function generator, or play back audio files with a multitude
of effects and filters.

Download the most recent version at http://www.baudline.com/download.html


Installation
------------
Copy the "baudline" executable to a suitable place such as /usr/local/bin,
/usr/bin, /home/user_name/bin, your home directory, or leave it right where it
is.  After a "rehash" run "baudline" or "./baudline" if dot isn't in your path.

For operation with the JACK Audio Connection Kit sound server use the
"baudline_jack" executable version and add the -jack command line option the
first time you run baudline.  Running "baudline_jack -jack" will setup JACK as
baudline's default audio device.

Use the included icons/ directory for creating baudline launch icons on your
desktop.


Instructions
------------
Running "baudline" will either start it up in record, pause, or play mode
depending on how you last quit the program.  Use the third mouse button to
access the main popup menu and control all of baudline's features.

Type "baudline -help" for a list of command line options.  These options are
meant as convenience controls and are not necessary for normal operation of
baudline.

For detailed instructions see the on-line manual at http://baudline.com/manual

Don't be afraid to experiment with different settings because baudline can
always be restored to the initial factory default state by starting it with
the command "baudline -reset"

Remember all the hidden power of this program is in the third mouse button.


Power Users
-----------
Install the file loading helper applications to open mp3, ogg, flac, adpcm,
and gsm audio files. See http://baudline.com/faq.html#input_decoder_helpers

Install more user color palettes by copying the included palettes/ directory
to your ~/.baudline/ directory. For more information about color palettes see
http://baudline.com/manual/color_picker.html#color_picker

Enable backing store for more video performance.  Setup baudline to be your
web browser's audio helper of choice.  Run multiple instances of baudline and
discover how the -session command line option works.  Learn some of the mouse
button and keyboard hot key short cuts. Particularly the the click, shift, and
drag measurements and the X / Y axis zooming with the Alt+<arrow keys>.  For
details and links see http://baudline.com/faq.html#power_user


License Agreement
-----------------
* This software is free and it comes with no warranty. 
* We are not liable for any damage caused by the use of this product. 
* You are not allowed to distribute this software. 
* You are not allowed to reverse engineer this software. 
* By downloading the baudline .tar package you agree to the terms of our
  license agreement. 
* If you desire a warranty on this product and you wish to purchase a support
  contract then please contact us.


Contact
-------
If you find a bug in baudline, it crashes, or behaves oddly then we would be
very interested in hearing about it.  In such a case please send us some email
at contact@baudline.com which describes the problem.  Also the following
information will be very useful to us when it comes to tracking down the
perplexity:

* text output from the command "baudline -sysinfo" 
* name and version number of your distribution
* name of the window manager you are using
* make and model of your soundcard

If you like using baudline then please keep in mind that SigBlips is a contract
engineering design and development firm that can be hired to perform custom
baudline modifications, add features, and for general signal analysis
consultation.  If you require baudline licensing, services, or support then
check out the http://www.sigblips.com website. 

------------------------------------------------------------------------------
Copyright (c) 2010 by SigBlips


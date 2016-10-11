# StructureFromSound
Code for the Structure from Sound project (SFS)
The aim of the project is to develop tools for estimating microphone positions
and sound source movements using a set of audio files as input only.

This is a challenging parameter estimation problem. 
For the project we also develop tools that can be of use for other TDOA type
problems.

## Code

Run main.m

## Data

Example data for the system can be downloaded from the
StructureFromSoundDatabase (sfsdb)
at [http://vision.maths.lth.se/sfsdb/](http://vision.maths.lth.se/sfsdb/)

Briefly the datasets are:

1.  bassh1 - lunch room (3DR, 3DS, Cont, Echo, Single, Distinct)
2.  bassh2 - lunch room (3DR, 3DS, Cont, Echo, Single, Distinct)
3.  bassh3 - lunch room (3DR, 3DS, Cont, Echo, Single, Distinct)
4.  grieg1 - lecture room (3DR, 3DS, Cont, Echo, Single, Distinct)
5.  grieg2 - lecture room (3DR, 3DS, Cont, Echo, Single, Distinct)
6.  grieg3 - lecture room (3DR, 3DS, Cont, Echo, Single, Distinct)
7.  grieg4 - lecture room (3DR, 3DS, Cont, Echo, Single, Distinct)
8.  grieg5 - lecture room (3DR, 3DS, Cont, Echo, Single, Distinct)
9.  spacecure1 - anechoic chamber (3DR, 3DS, Cont, NoEcho, Single, Distinct)
10. whywereyouawayayearroy - anechoic chamber (3DR, 3DS, Cont, NoEcho, Multiple, NotDistinct)
11. axelf - kitchen (2DR, 3DS, Cont, Echo, Single, Distinct)
12. gone - kitchen (2DR, 3DS, Cont, Echo, Single, Distinct)

Here
* 3DR - Dimension of affine hull of receivers is 3. (Microphones not in a plane)
* 2DR - Dimension of affine hull of receivers is 2. (Microphones in a plane)
* 1DR - Dimension of affine hull of receivers is 1. (Microphones on a line)
* (0DR - Dimension of affine hull of receivers is 0. (Microphones all in the same position)
* 3DS- Dimension of affine hull of sound events is 3. (Sound events not in a plane)
* 2DS- Dimension of affine hull of sound events is 2. (Sound events in a plane)
* 1DS - Dimension of affine hull of sound events is 1. (Sound events on a line)
* (0DS Dimension of affine hull of sound events is 0. (Sound events all in the same position)
* Claps - Claps, a few distinct sound events
* Cont - Continuous sound events
* Echo - In a normal (reverberant) room
* NoEcho - In anaechoic chamber
* Single - Single sound source
* Multiple - Multiple sound sources
* Distinct - Spatially distinct sound sources 
* NotDistinct - Spatially spread out sound sources

### Experiments in lunch room  September 18, 2014

The microphone setup is illustrated in the following figure. 

![Lunch Room](/tex/images/IMG_2283.JPG "Lunch Room")

We made 3 experiments (bassh1, bassh2, bassh3). There are several still 
images of the microphone setup and one stationary few film recording of the
experiments.

### Experiments in lecture room on November 27, 2014

The microphone setup is illustrated in the following figure. 

![Lecture Room](/tex/images/IMG_3442.JPG "Lecture Room")

We made 5 experiments (greig1 to grieg5). There are several still 
images of the microphone setup and a few film recordings of the
experiment. During theses film recordings, the camera was moving.

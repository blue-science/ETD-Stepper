# ETD Stepper

ETD Stepper, short for Exponential Time Difference Stepper, is a general
C++ code base for modeling the evolution of reaction-diffusion system
with multiple reactants in 1, 2, and 3 dimensions of space.

## Requirements

 * C++11 compiler
 * FFTW3 http://www.fftw.org/
 * Boost http://www.boost.org/
 * Standard Unix/Linux tools (bash, zcat, etc.)
 * For Image Processing:
    * Ruby http://www.ruby-lang.org/
    * GnuPlot http://www.gnuplot.info/
    * FFmpeg http://www.ffmpeg.org/

Most of these come standard on various Linux distros. For example,
Slackware Linux has all of these installed as part of the normal
installation with the exception of FFmpeg.

## Installation

It is recommended to create a new directory for each "experiment".
Though many settings are saved in an editable text file, a number of
settings, such as the initial conditions, are specified in the code,
and thus require separate compilations for different initial conditions.

    make

## Usage

    ./model_g

This can take anywhere from minutes to days. The first time this runs, it
will create an FFTW "wisdom file" called `wisdom-fftw.txt`. This can take
a variable amount of time, which subsequent runs can skip by reading the
wisdom-fftw.txt file. Next, it will spend some time pre-calculating the
`phi` values which are an effective cache to be using during the main
time-stepping loop. Finally, it will run the time-stepping loop where
it generally spends the majority of its time. Depending on the value of
the frames per unit time (fput) set in the config file (model\_g.txt by
default) it will output a frame file in a frames subdirectory, which is
later used for image processing.

Next run

    ./mkvideo.sh

which will print how to use it. This will read from the frames generated
in the previous step, generate a png image for each frame, and put them
together into an mp4 video.

An example experiment, in which the video was put together using
Mathematica instead of GnuPlot:

<http://www.blue-science.org/articles/2014/06/10/subquantum-kinetics-autogenesis-in-2-dimensions/>

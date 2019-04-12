# moch
Biophysical model of the mechanisms underlying pitch processing in auditory cortex

## Installation
These instructions are written for UNIX, but tranferring them to a windows system should be trivial. 
The libs were tested in both, Linux and Windows.

This libs use both MATLAB and python. MATLAB is the main wrapper and everything can be runned from 
there. Python is used to compute the output of the peripheral system and the cortical input. The
peripheral system runs from a external library contained in the pip package cochlea, which in turn
ses thorns to adjust the stimulus loudness. You will need both. The code should work fine in both,
in python 2.7 and python 3.

First, install cochlea and thorns using pip from the terminal: `pip install thorns cochlea`

Next, set line 54 in tdoch.m to point to your python installation. Just "python" should work in most
systems. I wrote a few extra tips in the comments preceding the line in case that does not work. You
ust be able to run python and import cochlea from matlab using `system(pythonCommandString)`

That's it!

## Usage
The model wrapper is tdoch.m. Thorought examples reproducing the figures shown in the original paper
re available in the folder "figures". Some easy UNIX-to-windows modifications might be necessary to
make the functions work in windows O:).

Have fun!

## Citation
If you use the model, please cite https://doi.org/10.1371/journal.pcbi.1006820

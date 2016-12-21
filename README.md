# WeakRates
Weak Interaction rates in Core-Collapse supernovae

## Installation
 * Install dependencies :

```apt-get install libhdf5-dev libgsl-dev```

 * Edit the Makefiles for the appropriate location of libhdf5 headers and libs
```
Makefile
compose/Makefile
```
 * Run ```./install.sh``` to compile the sources

## Usage

 * Extract compose files (eos.compo, eos.t, etc.) in a directory, preserving their original names.
 * Execute the program, passing compose files directory as the first and only argument (without trailing slash) :
 
```
 ./run path/to/compose/files
```

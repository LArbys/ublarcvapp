# DLLEE Cosmic Tagger

This folder contains the code for the 'classic' tagger.

We've stripped out much of the cruft from the version that lived
in `dllee_unified/larlitecv/app` 

The cosmic tagger consists of the following components:

1) apply PMT precuts
2) apply ROI based on opflash position
3) find end points
4) thrumu tagger connecting cadidate endpoints


# Running the Tagger

An executable which runs the `TaggerProcessor` class can be found
at `bin/run_tagger.cxx`.

The program creates an instance of the class, `TaggerProcessor`,
which deploys the various collection of Tagger algorithms.
It inherits from `larcv::Process`, allowing it to be folded into
a `larcv::ProcessDriver`.

## Inputs

The `larcv` data products needed are
* `Image2D` containing the deconvolved wire response for whole detector view
* `ChStatus` containing the status of the channels

The `larlite` data products needed are
* `ophit` individual optical photon hits for both in and out of the beam readout window
* `opflash` cluster of optical photon hits for both in and out of the beam readout window
* `trigger` event trigger information

## Outputs

## Inspecting the outputs for a given event

To inspect ouputs of each stage use the following notebook ...

## 

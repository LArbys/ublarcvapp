# Tagger Version 2

We aim to use the following inputs:

* mask-rcnn
* infill (within ROI)
* larflow (within ROI)

Mask R-CNN and Infill is our primary input. LArFlow used if the net is ready.

Goal is to match Mask R-CNN clusters (and use the tracker or A-STAR) to give us something we can reject/flash-match.

## Matching
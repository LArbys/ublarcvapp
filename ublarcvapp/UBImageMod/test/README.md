Test/debug scripts for the contents of UBImageMod

## UBSplitDetector

Provides a way to produce regular crops that breaks apart a full UB image in a 3D-consistent manner.
Default settings produce 832x512 images for each plane.
The center of each image are centered around the same 3D position in the detector.
The middle 310 wires on the Y-plane are gaurenteed to have some wire in both the U and V plane images
that it intersects with.

Test script:

    test_ubsplit.py [input larcv1 larcvtruth file] [adc producer name, typically "wiremc"]


## UBCropLArFlow

Takes the cropped images and ROIs defined in the `UBSplitDetector` and crops out the larflow
truth images.
The algorithm needs to adjust the flow values from the full image for the cropped images.

Test script:

    test_croplarflow.py [input larcvtruth image]
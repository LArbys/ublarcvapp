Scripts for DLTagger Development


## `test_mrmcnnmatch.py`

Use to test the MRCNNMatch class which is responsible for matching Mask R-CNN clusters across planes.
It uses mutual overlap in time and detector-z along with
`AStar3DAlgo` to check if a 3D-consistent path through charge can be defined.

Example out:

![](all_combos_pass1.png?raw=true)
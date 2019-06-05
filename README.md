# ublarcvapp: UB LArCV applications

The purpose of this repository is to house various tools and algos
that rely on LArCV (and larlite) data products and are used by
several different applications.

## Tool List

| tool                   | description | dependencies (other than larcv) |
|:----------------------:|:----------- | :------------|
| LArliteHandler         | syncronize larlite data file with larcv file using run,subrun,event indices | larlite |
| ContourTools           | create 2d charge clusters which are simply connected | larlite, geo2d, laropencv  |
| UBImageMod             | ub-specific image manipulation for networks (e.g. splitting the image) | larlite |
| UBWireTool             | code to query wire positions or calculate wire intersections | larlite |
| Filter                 | select events using truth for analysis | larlite |
| LArOpenCVHandle        | algorithm chain primarily dedicated to vertex finding | larlite, geo2d, laropencv |
| Stitchers              | consolidate network output produced on subimages | larlite |
| ubdllee                | tools dedicated to ub DLLEE analysis | larlite |




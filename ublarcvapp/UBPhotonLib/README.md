# UB Photon Library Interface


Where can I get the photon library file?

On the uboone CVMFS:

```
/cvmfs/uboone.opensciencegrid.org/products/uboone_photon_propagation/v01_01_00/PhotonPropagation/LibraryData/
```

This folder contains:
```
README                             uboone_photon_library_v6_70kV_EnhancedExtraTPCVis.root                   uboone_photon_library_variation_Rayleigh.root
uboone_photon_library_v5_0kV.root  uboone_photon_library_v6_70kV.root                                       uboone_photon_library_variation_Suppression75percent.root
uboone_photon_library_v5.root      uboone_photon_library_variation_Attenuation8m.root
uboone_photon_library_v6_0kV.root  uboone_photon_library_variation_Attenuation8m_Suppression75percent.root
```

The README in the folder contains the following info:
```
Available photon libraries
--------------------------

uboone_photon_library_v5_0kV.root  - MCC7 no E-field
uboone_photon_library_v5.root      - MCC7 70kV (this file used for mcc7 g4)
uboone_photon_library_v6_0kV.root  - MCC8 no E-field
uboone_photon_library_v6_70kV.root - MCC8 70kV (standard hv)
uboone_photon_library_v6_70kV_EnhancedExtraTPCVis.root
                                   - Increased visibility outside TPC by 50%.

The zero E-field photon libraries contain the probability that a
photon generated in a particular voxel will be detected by any pmt.

The nonzero E-field are modified from the zero E-field photon libraries
by including the relative effect of the E-field on photon yield due
to recombination.
```


Of the options, we use for now: `uboone_photon_library_v6_70kV_EnhancedExtraTPCVis.root`

This file is 101 MB.


Inside the file is the following ROOT tree `PhotonLibraryData` located in the TDirectory, `pmtresponse`.
```
******************************************************************************
*Tree    :PhotonLibraryData:                                                        *
*Entries : 33615139 : Total =       403773735 bytes  File  Size =  105069173 *
*        :          : Tree compression factor =   3.84                       *
******************************************************************************
*Br    0 :Voxel     : Voxel/I                                                *
*Entries : 33615139 : Total  Size=  134587491 bytes  File Size  =    5521593 *
*Baskets :     1206 : Basket Size=    3264000 bytes  Compression=  24.37     *
*............................................................................*
*Br    1 :OpChannel : OpChannel/I                                            *
*Entries : 33615139 : Total  Size=  134592331 bytes  File Size  =    7064466 *
*Baskets :     1206 : Basket Size=    3264000 bytes  Compression=  19.05     *
*............................................................................*
*Br    2 :Visibility : Visibility/F                                          *
*Entries : 33615139 : Total  Size=  134593541 bytes  File Size  =   92451277 *
*Baskets :     1206 : Basket Size=    3264512 bytes  Compression=   1.46     *
*............................................................................*
```

Voxels are labeled by a simple integer index.
This means we need to know the voxelization defintion and where it comes from.
We also need to know the OpChannel labeling scheme.


Parseing the Table

Number of entries in the tree is: 33615139.

Found the following config in the `standard_g4_uboone.fcl`:

```
   PhotonVisibilityService: {
      DoNotLoadLibrary: false
      Interpolate: false
      LibraryBuildJob: false
      LibraryFile: "PhotonPropagation/LibraryData/uboone_photon_library_v6_70kV.root"
      NX: 75
      NY: 75
      NZ: 400
      UseCryoBoundary: true
      XMax: 120
      XMin: -120
      YMax: 120
      YMin: -120
      ZMax: 1200
      ZMin: 0
      service_type: "PhotonVisibilityService"
   }
```
That implies a total of 75*75*400*32 = 72,000,000 possible entries.
So less than half is filled. I guess that makes sense.

The bounds don't make total sense.
The X width is only 240 cm. The TPC drift length is 256 cm!
The Y width is 240 cm, which is fine. (TPC height is 2*117 cm)
The Z width is 1200 cm, which is also fine. (TPC is 1036 cm)
OK, I see in `photpropservices.fcl` that these bounds are not used when "UseCryoBoundary: true".

So then what are the cryo boundaries?
The photonvisibility service gets it from the geometry surface, by calling GetCryostat.
Great, so how does it do this?

Used dump_geometry to extract cryostat and TPC bounds.

The Key info (dump with all relevant det geom below).


```
Cryostat C:0 (383.22 x 383.22 x 1221.75) cm^3 at (128.175,0,518.5)
TPC C:0 T:0 (260 x 256 x 1045) cm^3 at (128.175,0.97,518.5)
```

We work in the TPC coordinates always. Looks like the center is mostly the same.
There is a y=-0.97 offset (cryostat center is lower in global coordinates).
Dimensions x: 383.22 cm/75 = 5.1096 cm voxels
Origin of the C in global is: (-63.435,-191.61,-92.375)
Origin of TPC in the global is: (-1.825,-127.03,-4)
We use origin_cryo-origin_tpc to get cryo origin in TPC coordinate system.



```
Detector description: '/cvmfs/uboone.opensciencegrid.org/products/ubcore/v08_00_00_40/gdml/microboonev12.gdml'
Detector microboonev12 has 1 cryostats and 1168 auxiliary detectors:
  Detector enclosure: (-613.455,-895.31,-223.13) -- (869.805,895.31,1260.13) cm => ( 1483.26 x 1790.62 x 1483.26 ) cm^3
  Cryostat C:0 (383.22 x 383.22 x 1221.75) cm^3 at (128.175,0,518.5)
    hosts 1 TPCs (largest number of planes: 3, of wires: 3456) and 32 optical detectors
    bounding box: (-63.435,-191.61,-92.375) -- (319.785,191.61,1129.38)
    TPC C:0 T:0 (260 x 256 x 1045) cm^3 at (128.175,0.97,518.5)
      drift direction (-1,0,0) from cathode around (254.8,0.97,518.5) through 255.4 cm toward 3 wire planes
      maximum wires on any plane: 3456
      active volume (256.35 x 233 x 1036.8) cm^3, front face at (126.625,0.97,0.1) cm;
      main directions: width (1,0,0) height (0,1,0) length (0,0,1)
      bounding box: (-1.825,-127.03,-4) -- (258.175,128.97,1041)
      active volume box: (-1.55,-115.53,0.1) -- (254.8,117.47,1036.9)
      plane C:0 T:0 P:0 at (-6.34622e-14,0.97,518.5) cm, theta: 0.523599 rad
        normal to wire: -1.0472 rad, with orientation vertical, has 2400 wires measuring U with a wire pitch of 0.3 cm
        normal to plane: (1,0,0), direction of increasing wire number: (0,-0.866025,0.5) [wire frame normal: (1,0,0)] (increases with z)
        wire direction: (0,0.5,0.866025); width 1037.35 cm in direction: (0,0,-1), depth 233 cm in direction: (0,1,0) [normal: (1,0,0)]
        wires cover width -518.411 to 518.411, depth -116.387 to 116.387 cm
        bounding box: (-0.075,-115.53,-0.175) -- (0.075,117.47,1037.17)
      ... wire stuff ...
      plane C:0 T:0 P:1 at (-0.3,0.97,518.5) cm, theta: 2.61799 rad
        normal to wire: 1.0472 rad, with orientation vertical, has 2400 wires measuring V with a wire pitch of 0.3 cm
        normal to plane: (1,0,0), direction of increasing wire number: (0,0.866025,0.5) [wire frame normal: (1,0,0)] (increases with z)
        wire direction: (0,0.5,-0.866025); width 1037.35 cm in direction: (0,0,1), depth 233 cm in direction: (0,-1,0) [normal: (1,0,0)]
        wires cover width -518.411 to 518.411, depth -116.387 to 116.387 cm
        bounding box: (-0.375,-115.53,-0.175) -- (-0.225,117.47,1037.17)
      plane C:0 T:0 P:2 at (-0.6,0.97,518.5) cm, theta: 1.5708 rad
        normal to wire: 0 rad, with orientation vertical, has 3456 wires measuring Z with a wire pitch of 0.3 cm
        normal to plane: (1,0,0), direction of increasing wire number: (0,0,1) [wire frame normal: (1,0,0)] (increases with z)
        wire direction: (0,1,0); width 1037 cm in direction: (0,0,1), depth 233 cm in direction: (0,-1,0) [normal: (1,0,0)]
        wires cover width -518.4 to 518.4, depth -116.5 to 116.5 cm
        bounding box: (-0.675,-115.53,0) -- (-0.525,117.47,1037)
      ... wire stuff ...
      
    [OpDet #0] centered at (-11.4545,-28.625,990.356) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #1] centered at (-11.4175,27.607,989.712) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #2] centered at (-11.7755,-56.514,951.865) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #3] centered at (-11.6415,55.313,951.861) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #4] centered at (-12.0585,-56.309,911.939) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #5] centered at (-11.8345,55.822,911.065) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #6] centered at (-12.1765,-0.722,865.599) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #7] centered at (-12.3045,-0.502,796.208) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #8] centered at (-12.6045,-56.284,751.905) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #9] centered at (-12.5405,55.625,751.884) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #10] centered at (-12.6125,-56.408,711.274) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #11] centered at (-12.6615,55.8,711.073) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #12] centered at (-12.6245,-0.051,664.203) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #13] centered at (-12.6515,-0.549,585.284) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #14] centered at (-12.8735,55.822,540.929) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #15] centered at (-12.6205,-56.205,540.616) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #16] centered at (-12.5945,-56.323,500.221) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #17] centered at (-12.9835,55.771,500.134) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #18] centered at (-12.6185,-0.875,453.096) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #19] centered at (-13.0855,-0.706,373.839) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #20] centered at (-12.6485,-57.022,328.341) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #21] centered at (-13.1865,54.693,328.212) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #22] centered at (-13.4175,54.646,287.976) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #23] centered at (-13.0075,-56.261,287.639) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #24] centered at (-13.1505,-0.829,242.014) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #25] centered at (-13.4415,-0.303,173.743) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #26] centered at (-13.3965,55.249,128.354) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #27] centered at (-13.2784,-56.203,128.18) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #28] centered at (-13.2375,-56.615,87.8695) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #29] centered at (-13.5415,55.249,87.7605) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #30] centered at (-13.4345,27.431,51.1015) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
    [OpDet #31] centered at (-13.1525,-28.576,50.4745) cm, radius: 15.24 cm, length: 0.2 cm, theta(z): 1.5708 rad
```

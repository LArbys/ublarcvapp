# MC Truth Tools

Module contains functions and classes to more easily parse the truth information from MicroBooNE files.

# Modules

## MCPixelPGraph

Parses the mctrack and mcshower trees (larlite). Makes a graph using the motherid.
Then attempts to associate pixels to the particles in the graph.

Example labeling (1e1p + corsika cosmic image)

![Example image](https://raw.githubusercontent.com/LArbys/ublarcvapp/master/ublarcvapp/MCTools/test/mcpg_example_1e1p_and_cosmic.png)


And portions of the graph.

The image above was made using the example script found in the `test` folder for this module, `run_mcpixelpgraph.py`.

## crossingPointsAnaMethods

Static function, `getFirstStepPosInsideImage`, is used to get the starting (and stopping) point of particles.

Example of end-point labels (corsiak cosmic image)

![Example start and stop points](https://raw.githubusercontent.com/LArbys/ublarcvapp/master/ublarcvapp/MCTools/test/demo_crossingPointAnaMethods.png)

Image above made using example script found in the `test` folder: `run_crossingpointsanamethods.py`

## LArbysMC

Use this to calculate truth variables useful for analysis. Defines target final states.

You can add the truth variables calculate to any TTree.

... before event loop ...

```
TTree* t = new TTree("mytree", "tree we add analysis variables to");

larcv::LArbysMC lmc;

lmc.bindAnaVariables( t );
```

... inside event loop ...

```
lmc.process( iolcv, ioll );
t.Fill();
```

Note, if you dont have access to larcv products or do not want to calculate vertex activity variables,
use `lmc.process( ioll )` instead.
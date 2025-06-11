# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ublarcvapp is a C++ framework with Python bindings for analyzing Liquid Argon Time Projection Chamber (LArTPC) data from the MicroBooNE neutrino physics experiment. It provides tools for:
- Deep learning-based particle tagging and classification
- 3D reconstruction algorithms
- Image manipulation and processing
- Monte Carlo truth analysis
- Event filtering and selection

## Build Commands

### Prerequisites Setup
Ensure these environment variables are set:
- `LARCV_LIBDIR` - LArCV library directory
- `LARLITE_LIBDIR` and `LARLITE_INCDIR` - larlite directories
- `LARCV_INCDIR` - LArCV include directory
- `LAROPENCV_LIBDIR` and `LAROPENCV_INCDIR` - LArOpenCV directories
- Optional: `GEO2D_BASEDIR` - for Geo2D support

### Building
```bash
mkdir build
cd build
cmake ..
make -j4
make install
source ../configure.sh
```

### Running Tests
Tests are located in individual module `test/` directories. Example:
```bash
cd ublarcvapp/Reco3D/test
python test_astar.py

cd ublarcvapp/UBImageMod/test
python test_ubsplit.py [input_larcv_file] [adc_producer_name]
```

## Code Architecture

### Core Components

1. **Processing Framework** (`LLCVProcessor/`)
   - `LLCVProcessBase` - Base class for all processing modules
   - `LLCVProcessDriver` - Manages processing chain execution
   - Configuration via `.cfg` files

2. **Image Processing** (`UBImageMod/`)
   - `UBSplitDetector` - Splits detector images into subregions
   - `UBCropLArFlow` - Crops sparse larflow data
   - `InfillImageStitcher` - Stitches subimages back together

3. **3D Reconstruction** (`Reco3D/`)
   - `AStar3DAlgo` - A* pathfinding for 3D track reconstruction
   - Configuration through `AStar3DAlgoConfig`
   - Supports both regular and proton-specific tracking

4. **Deep Learning Integration** (`DLTagger/`)
   - `DLTagger` - Main interface for DL model predictions
   - `MRCNNMatch` - Mask R-CNN matching algorithms
   - Integrates with external DL model servers

5. **Monte Carlo Tools** (`MCTools/`)
   - Truth matching utilities
   - Space charge effect corrections
   - Flash matching algorithms

### Key Design Patterns

1. **Factory Pattern**: Most algorithms use factory registration:
   ```cpp
   static LArbysImageFactory __global_LArbysImageFactory__;
   ```

2. **Configuration**: Algorithms use dedicated config classes:
   ```cpp
   class AlgoConfig : public larcv::PSet {
     // Configuration parameters
   };
   ```

3. **Process Chain**: Modules inherit from `LLCVProcessBase` and implement:
   - `configure()` - Setup from config file
   - `initialize()` - One-time initialization
   - `process()` - Per-event processing
   - `finalize()` - Cleanup

### Data Flow

1. Input: ROOT files containing larcv::EventImage2D, larlite::event_* objects
2. Processing: Chain of LLCVProcess modules
3. Output: Modified/new data products saved to output ROOT files

## Development Guidelines

### Adding New Modules

1. Create directory under `ublarcvapp/`
2. Add CMakeLists.txt with:
   ```cmake
   set(MODULE_NAME YourModule)
   add_subdirectory(${MODULE_NAME})
   ```
3. Implement classes inheriting from `LLCVProcessBase`
4. Add LinkDef.h for ROOT dictionary generation
5. Register with factory if needed

### Working with Sparse Data

Many algorithms work with sparse representations:
- Use `larcv::EventSparseTensor2D` for sparse image data
- Convert between dense/sparse with utility functions
- Handle coordinate transformations carefully

### Coordinate Systems

Be aware of multiple coordinate systems:
- Wire/tick 2D coordinates per plane
- 3D spatial coordinates (with/without space charge effects)
- Image pixel coordinates
- Detector geometry coordinates
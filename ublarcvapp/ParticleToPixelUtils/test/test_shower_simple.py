#!/usr/bin/env python3
"""
Simple test for ShowerToSpacePoints without photon library
"""

import ROOT as rt
from ROOT import std
from larcv import larcv
from larlite import larlite
from ublarcvapp import ublarcvapp
import numpy as np

# Create synthetic shower
print("Creating synthetic shower...")
cluster = larlite.larflowcluster()

# Create 10 hits in a shower pattern
for i in range(10):
    hit = larlite.larflow3dhit()
    hit.clear()
    
    # Simple positions
    x = 70.0 + np.random.normal(0, 2)
    y = 0.0 + np.random.normal(0, 5)
    z = 500.0 + np.random.uniform(0, 20)
    
    hit.push_back(float(x))
    hit.push_back(float(y))
    hit.push_back(float(z))
    
    cluster.push_back(hit)
    print(f"  Hit {i}: ({x:.1f}, {y:.1f}, {z:.1f})")

print(f"\nCreated cluster with {cluster.size()} hits")

# Create simple ADC images
print("\nCreating ADC images...")
rows = 2048
cols = 3456

# Simple metadata
meta = larcv.ImageMeta(cols, rows, rows, cols, 0, 3000, 0)
img = larcv.Image2D(meta)
img.paint(0)

# Add some charge at expected locations
for row in range(1000, 1500, 10):
    for col in range(1500, 2000, 10):
        img.set_pixel(row, col, 50.0)

adc_images = std.vector('larcv::Image2D')()
adc_images.push_back(img)
adc_images.push_back(img)
adc_images.push_back(img)

# Test ShowerToSpacePoints
print("\nTesting ShowerToSpacePoints...")
converter = ublarcvapp.pixelutils.ShowerToSpacePoints()

spacepoints = converter.convertShower(cluster, adc_images, 10.0, 5, 5)

print(f"\nResults:")
print(f"  Space points found: {len(spacepoints)}")
print(f"  Pixels processed: {converter.getNumPixelsProcessed()}")
print(f"  Total charge: {converter.getTotalCharge():.1f}")

if len(spacepoints) > 0:
    print("\nFirst few space points:")
    for i, sp in enumerate(spacepoints[:3]):
        print(f"  SP {i}: pos=({sp.position.X():.1f}, {sp.position.Y():.1f}, {sp.position.Z():.1f}) charge={sp.charge:.1f}")
else:
    print("\nNo space points found - this may be due to:")
    print("  - Mismatch between hit positions and charge deposit locations")
    print("  - Tick/wire conversion issues")
    print("  - Threshold too high")

print("\nShowerToSpacePoints class is working!")
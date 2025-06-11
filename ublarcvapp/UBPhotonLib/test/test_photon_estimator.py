#!/usr/bin/env python
"""
Test script for PhotonVisibilityEstimator class
Demonstrates how to estimate photon detection from multiple 3D point sources
"""

import ROOT as rt
from ROOT import std
from larcv import larcv
from ublarcvapp import ublarcvapp

# Load the library
rt.gSystem.Load("libLArCVApp_UBPhotonLib.so")

def test_photon_estimator():
    """Test the PhotonVisibilityEstimator functionality"""
    
    # Create estimator instance
    estimator = ublarcvapp.ubphotonlib.PhotonVisibilityEstimator()
    
    print("Created PhotonVisibilityEstimator")
    print("Number of optical channels:", estimator.getNumOpticalChannels())
    
    # Example 1: Add individual point sources
    # These positions should be in detector coordinates (cm)
    # Let's add some points along a track
    print("\n--- Example 1: Individual point sources ---")
    
    # Add points along a hypothetical track
    # Starting near center of detector
    x_start, y_start, z_start = 128.0, 0.0, 500.0  # cm
    x_end, y_end, z_end = 128.0, 50.0, 600.0      # cm
    
    # Create 10 points along this track
    n_points = 10
    photons_per_point = 1000.0  # photons per point
    
    for i in range(n_points):
        t = float(i) / (n_points - 1)
        x = x_start + t * (x_end - x_start)
        y = y_start + t * (y_end - y_start)
        z = z_start + t * (z_end - z_start)
        
        estimator.addPhotonSource(x, y, z, photons_per_point)
    
    print("Added {} point sources".format(estimator.getNumSources()))
    print("Total emitted photons: {}".format(estimator.getTotalEmittedPhotons()))
    
    # Calculate detected photons
    photons_per_opdet = estimator.calculateDetectedPhotons(use_trilinear=True)
    
    # Print results
    print("\nPhotons detected by each optical detector:")
    total_detected = 0.0
    for opch in range(estimator.getNumOpticalChannels()):
        n_photons = photons_per_opdet[opch]
        if n_photons > 0:
            print("  OpDet {:2d}: {:8.2f} photons".format(opch, n_photons))
            total_detected += n_photons
    
    print("\nTotal detected photons: {:.2f}".format(total_detected))
    print("Collection efficiency: {:.4f}".format(estimator.getCollectionEfficiency(True)))
    
    # Example 2: Using PhotonSource struct
    print("\n--- Example 2: Using PhotonSource struct ---")
    
    # Clear previous sources
    estimator.clear()
    
    # Create a vector of sources
    sources = std.vector('ublarcvapp::ubphotonlib::PhotonVisibilityEstimator::PhotonSource')()
    
    # Add sources in a different pattern - a shower-like distribution
    center_x, center_y, center_z = 100.0, 0.0, 400.0
    spread = 20.0  # cm
    
    import random
    random.seed(42)
    
    for i in range(20):
        # Random positions around center
        x = center_x + random.gauss(0, spread)
        y = center_y + random.gauss(0, spread)
        z = center_z + random.gauss(0, spread)
        
        # More photons near center
        distance = ((x-center_x)**2 + (y-center_y)**2 + (z-center_z)**2)**0.5
        n_photons = 2000.0 * (1.0 - distance / (3*spread))
        if n_photons < 0:
            n_photons = 0
            
        source = ublarcvapp.ubphotonlib.PhotonVisibilityEstimator.PhotonSource(x, y, z, n_photons)
        sources.push_back(source)
    
    # Add all sources at once
    estimator.addPhotonSources(sources)
    
    print("Added {} shower-like point sources".format(estimator.getNumSources()))
    print("Total emitted photons: {:.2f}".format(estimator.getTotalEmittedPhotons()))
    
    # Calculate with and without trilinear interpolation
    print("\nWith trilinear interpolation:")
    total_trilinear = estimator.getTotalDetectedPhotons(use_trilinear=True)
    print("  Total detected: {:.2f}".format(total_trilinear))
    print("  Collection efficiency: {:.4f}".format(estimator.getCollectionEfficiency(True)))
    
    print("\nWithout trilinear interpolation:")
    total_no_trilinear = estimator.getTotalDetectedPhotons(use_trilinear=False)
    print("  Total detected: {:.2f}".format(total_no_trilinear))
    print("  Collection efficiency: {:.4f}".format(estimator.getCollectionEfficiency(False)))
    
    # Example 3: Visualize detection pattern
    print("\n--- Example 3: Detection pattern visualization ---")
    
    # Get photons per detector for visualization
    photons_per_opdet = estimator.calculateDetectedPhotons(use_trilinear=True)
    
    # Create a simple histogram
    h_opdet = rt.TH1F("h_opdet", "Photons Detected per Optical Detector;OpDet ID;Photons", 
                      32, 0, 32)
    
    for opch in range(32):
        h_opdet.SetBinContent(opch+1, photons_per_opdet[opch])
    
    # Draw if running interactively
    canvas = rt.TCanvas("c1", "Photon Detection", 800, 600)
    h_opdet.SetFillColor(rt.kBlue-9)
    h_opdet.Draw("HIST")
    canvas.SaveAs("photon_detection_pattern.png")
    print("Saved histogram to photon_detection_pattern.png")

if __name__ == "__main__":
    test_photon_estimator()
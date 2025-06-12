#!/usr/bin/env python3
"""
Complete example demonstrating the full workflow from track to light detection estimation.

This script shows:
1. How to create realistic larlite::track objects from line segments
2. How to create realistic ADC images with proper geometry and charge deposits
3. How to use TrackToSpacePoints to convert tracks to space points with charge
4. How to convert charge to photons using LAr physics
5. How to use PhotonVisibilityEstimator to predict light detection

Usage:
    python3 complete_track_light_example.py [--input-file file.root] [--save-output]
"""

import ROOT as rt
from ROOT import std
import numpy as np
import argparse
import sys

# Load libraries
from larcv import larcv
from larlite import larlite  
from ublarcvapp import ublarcvapp

def create_track_from_segments(segments):
    """
    Create a larlite::track from a list of 3D line segments
    
    Args:
        segments: List of tuples [(start_point, end_point), ...] 
                 where each point is (x, y, z) in cm
    
    Returns:
        larlite::track object
    """
    track = larlite.track()
    
    print(f"Creating track from {len(segments)} line segments")
    
    for i, (start, end) in enumerate(segments):
        # Add start point (skip if it's the same as previous end point)
        if i == 0:
            pos = rt.TVector3(start[0], start[1], start[2])
            track.add_vertex(pos)
        
        # Add end point
        pos = rt.TVector3(end[0], end[1], end[2])
        track.add_vertex(pos)
        
        # Calculate direction vector
        dx = end[0] - start[0]
        dy = end[1] - start[1] 
        dz = end[2] - start[2]
        direction = rt.TVector3(dx, dy, dz)
        if direction.Mag() > 0:
            direction = direction.Unit()
        
        track.add_direction(direction)
        
        print(f"  Segment {i}: ({start[0]:.1f}, {start[1]:.1f}, {start[2]:.1f}) -> "
              f"({end[0]:.1f}, {end[1]:.1f}, {end[2]:.1f})")
    
    return track

def create_realistic_adc_images():
    """
    Create realistic ADC Image2D objects with proper MicroBooNE geometry
    
    Returns:
        std::vector<larcv::Image2D> with U, V, Y plane images
    """
    print("Creating realistic ADC images with MicroBooNE geometry...")
    
    # MicroBooNE detector parameters
    tick_period = 0.5  # microseconds per tick
    rows = 6048        # time ticks
    
    # Wire plane parameters
    cols_u = 2400      # U plane wires
    cols_v = 2400      # V plane wires
    cols_y = 3456      # Y plane wires (collection plane)
    
    # Create image metadata
    # ImageMeta(width, height, row_count, col_count, origin_x, origin_y, plane)
    wire_pitch = 0.3  # cm
    time_pitch = 0.5  # μs per tick
    
    width_u = cols_u * wire_pitch
    width_v = cols_v * wire_pitch  
    width_y = cols_y * wire_pitch
    height = rows * time_pitch
    
    meta_u = larcv.ImageMeta(width_u, height, rows, cols_u, 0.0, 0.0, 0)
    meta_v = larcv.ImageMeta(width_v, height, rows, cols_v, 0.0, 0.0, 1)
    meta_y = larcv.ImageMeta(width_y, height, rows, cols_y, 0.0, 0.0, 2)
    
    # Create empty images
    img_u = larcv.Image2D(meta_u)
    img_v = larcv.Image2D(meta_v) 
    img_y = larcv.Image2D(meta_y)
    
    # Fill with small baseline noise
    for img in [img_u, img_v, img_y]:
        img.paint(0.0)
        # Add realistic noise
        for row in range(0, rows, 50):  # Sample every 50 ticks for speed
            for col in range(0, img.meta().cols(), 20):  # Sample every 20 wires
                noise = np.random.normal(0, 2.0)  # 2 ADC RMS noise
                if noise > 0:
                    img.set_pixel(row, col, noise)
    
    print(f"  U plane: {cols_u} wires x {rows} ticks")
    print(f"  V plane: {cols_v} wires x {rows} ticks") 
    print(f"  Y plane: {cols_y} wires x {rows} ticks")
    
    # Package into vector
    adc_images = std.vector('larcv::Image2D')()
    adc_images.push_back(img_u)
    adc_images.push_back(img_v)
    adc_images.push_back(img_y)
    
    return adc_images

def get_image_and_tracks_from_rootfile(rootfile):
    """
    Get Image2D and track examples for testing TrackToSpacePoints 
    and the resulting photon prediction using ublarcvapp.ubphotonlib.PhotonVisibilityEstimator.
    The test file is a merged file containing TTrees with larcv objects and
    the reconstructed tracks.
    """
    print(f"Loading data from {rootfile}")
    
    # Load larcv images
    ioman = larcv.IOManager(larcv.IOManager.kREAD)
    ioman.add_in_file(rootfile)
    ioman.initialize()

    # Get the first entry
    ioman.read_entry(0)

    ev_img = ioman.get_data(larcv.kProductImage2D, "wire")
    imglist = []
    for p in range(3):
        img = ev_img.as_vector().at(p)
        imglist.append(img)
        print(f"  Plane {p}: {img.meta().cols()} x {img.meta().rows()} pixels")

    # Convert to std::vector for TrackToSpacePoints
    adc_images = std.vector('larcv::Image2D')()
    for img in imglist:
        adc_images.push_back(img)

    # Get the opflash data
    ioll = larlite.storage_manager(larlite.storage_manager.kREAD)
    ioll.add_in_filename(rootfile)
    ioll.open()
    ioll.go_to(0)
    
    # Get opflash data using correct casting
    ev_opflash = ioll.get_data(larlite.data.kOpFlash,"simpleFlashBeam")
    opflashes = std.vector('larlite::opflash')()
    
    # Try to access the opflash data
    print(f"  Found {ev_opflash.size()} optical flashes")
    for iopflash in range(ev_opflash.size()):
        opflashes.push_back(ev_opflash.at(iopflash))

    # Get the larlite::track examples from the reco file
    reco_tfile = rt.TFile(rootfile, 'open')
    recotree = reco_tfile.Get("KPSRecoManagerTree")
    trackout_v = std.vector('larlite::track')()
    
    if not recotree:
        print("Warning: Could not find KPSRecoManagerTree in file")
        return adc_images, trackout_v, opflashes
        
    recotree.GetEntry(0)
    
    # Check if we have the expected structure
    if not hasattr(recotree, 'nuvetoed_v') or recotree.nuvetoed_v.size() == 0:
        print("Warning: No nuvetoed_v found in tree")
        return adc_images, []
        
    nuvtx = recotree.nuvetoed_v.at(0)
    
    print(f"  Found {nuvtx.track_v.size()} tracks")


    for itrack in range(nuvtx.track_v.size()):
        track = nuvtx.track_v.at(itrack)
        if track.NumberTrajectoryPoints()>=2:
            trackout_v.push_back(track)
            print(f"    Track {itrack}: {track.NumberTrajectoryPoints()} points, length {track.Length():.1f} cm")

    ioman.finalize()
    reco_tfile.Close()
    
    return adc_images, trackout_v, opflashes

def add_track_charge_to_images(adc_images, track, charge_per_cm=1000.0):
    """
    Add realistic charge deposits along a track path in the ADC images
    
    Args:
        adc_images: Vector of larcv::Image2D objects
        track: larlite::track object
        charge_per_cm: Average charge deposition per cm of track
    """
    print(f"Adding track charge deposits ({charge_per_cm} ADC/cm)...")
    
    # MicroBooNE parameters
    drift_velocity = 0.1114  # cm/μs (at 273V/cm)
    tick_period = 0.5       # μs/tick
    x_offset = 3200         # ticks (approximate trigger offset)
    
    # Wire geometry (simplified - in reality would use larutil::Geometry)
    wire_pitch = 0.3        # cm
    
    n_deposits = 0
    total_charge = 0
    
    # Step along track and add charge
    n_points = track.NumberTrajectoryPoints()
    if n_points < 2:
        return
        
    for i in range(n_points - 1):
        start_pos = track.LocationAtPoint(i)
        end_pos = track.LocationAtPoint(i + 1)
        
        # Calculate segment length
        dx = end_pos.X() - start_pos.X()
        dy = end_pos.Y() - start_pos.Y()
        dz = end_pos.Z() - start_pos.Z()
        segment_length = np.sqrt(dx*dx + dy*dy + dz*dz)
        
        if segment_length < 0.1:  # Skip very short segments
            continue
            
        # Number of charge deposits along this segment
        n_steps = max(1, int(segment_length / 0.2))  # Every 2mm
        
        for step in range(n_steps):
            t = float(step) / n_steps
            x = start_pos.X() + t * dx
            y = start_pos.Y() + t * dy
            z = start_pos.Z() + t * dz
            
            # Convert 3D position to wire/tick coordinates
            tick = x / drift_velocity / tick_period + x_offset
            
            # Simplified wire calculation (in reality would use proper geometry)
            # U plane: 30° from vertical
            # V plane: -30° from vertical  
            # Y plane: vertical (collection)
            wire_u = int((y * np.cos(np.pi/6) + z * np.sin(np.pi/6)) / wire_pitch + 1200)
            wire_v = int((-y * np.cos(np.pi/6) + z * np.sin(np.pi/6)) / wire_pitch + 1200)
            wire_y = int(z / wire_pitch + 1700)
            
            row = int(tick)
            
            # Add charge if within image bounds
            wires = [wire_u, wire_v, wire_y]
            for plane in range(3):
                img = adc_images[plane]
                col = wires[plane]
                
                if (0 <= row < img.meta().rows() and 
                    0 <= col < img.meta().cols()):
                    
                    # Calculate charge for this step
                    step_length = segment_length / n_steps
                    base_charge = charge_per_cm * step_length
                    
                    # Add some fluctuation (Landau-like)
                    charge = np.random.gamma(2.0, base_charge/2.0)
                    
                    # Spread charge over neighboring pixels (diffusion effect)
                    for dr in range(-1, 2):
                        for dc in range(-2, 3):
                            r = row + dr
                            c = col + dc
                            if (0 <= r < img.meta().rows() and 
                                0 <= c < img.meta().cols()):
                                
                                # Gaussian spread
                                weight = np.exp(-(dr*dr + dc*dc) / 2.0)
                                pixel_charge = charge * weight * 0.1  # 10% of charge per pixel
                                
                                current_val = img.pixel(r, c)
                                img.set_pixel(r, c, current_val + pixel_charge)
                                
                                if plane == 2:  # Count only Y plane for statistics
                                    n_deposits += 1
                                    total_charge += pixel_charge
    
    print(f"  Added {n_deposits} charge deposits")
    print(f"  Total charge: {total_charge:.0f} ADC counts")

def convert_charge_to_photons(adc_charge):
    """
    Convert ADC charge to number of scintillation photons
    
    Args:
        adc_charge: Charge in ADC counts
        
    Returns:
        Number of scintillation photons
    """
    # LAr physics parameters (approximate)
    adc_per_electron = 200.0      # ADC counts per electron
    mev_per_electron = 23.6e-6    # Ionization energy in LAr (MeV)
    photons_per_mev = 24000.0     # Scintillation photons per MeV at 500V/cm
    recombination_factor = 0.7    # Fraction surviving recombination
    
    # Convert: ADC -> electrons -> energy -> photons
    n_electrons = adc_charge / adc_per_electron
    energy_mev = n_electrons * mev_per_electron
    n_photons = energy_mev * photons_per_mev * (1.0 - recombination_factor)
    
    return n_photons

def main():
    parser = argparse.ArgumentParser(description="Complete track to light estimation example")
    parser.add_argument("--input-file", help="ROOT file with track/image data (optional)")
    parser.add_argument("--save-output", action="store_true", help="Save output plots")
    parser.add_argument("--charge-per-cm", type=float, default=2000.0, 
                       help="Charge deposition per cm (ADC)")
    args = parser.parse_args()
    
    print("Complete Track to Light Estimation Example")
    print("="*60)
    
    opflashes = None  # Initialize for synthetic data case
    
    if args.input_file:
        # Use real data from ROOT file
        print(f"\nUsing real data from: {args.input_file}")
        adc_images, tracks, opflashes = get_image_and_tracks_from_rootfile(args.input_file)
        
        if len(tracks) == 0:
            print("No tracks found in input file, exiting")
            return
            
        # Use the first track for demonstration
        track = tracks[0]
        track_start = track.LocationAtPoint(0)
        print(f"\nUsing track 0: {track.NumberTrajectoryPoints()} points, length {track.Length():.1f} cm")
        print(f"Start: ({track_start.X():.1f}, {track_start.Y():.1f}, {track_start.Z():.1f})")
        print(f"End: ({track.End().X():.1f}, {track.End().Y():.1f}, {track.End().Z():.1f})")
        
        # Show information about optical flashes
        print(f"\nFound {len(opflashes)} optical flashes")
        if len(opflashes) > 0:
            flash = opflashes[0]  # Use first flash for comparison
            print(f"Using flash 0: Total PE = {flash.TotalPE():.1f}, Time = {flash.Time():.1f} μs")
        
    else:
        # Use synthetic data
        print("\nUsing synthetic data (use --input-file for real data)")
        
        # Define track as series of line segments
        track_segments = [
            # Start outside detector, enter from top
            ((10.0, 100.0, 200.0), (50.0, 80.0, 300.0)),     # entering
            ((50.0, 80.0, 300.0), (90.0, 50.0, 450.0)),      # through detector  
            ((90.0, 50.0, 450.0), (130.0, 20.0, 600.0)),     # continuing
            ((130.0, 20.0, 600.0), (170.0, -10.0, 750.0)),   # exiting
        ]
        
        # Step 1: Create track from line segments
        print("\n1. Creating track from line segments...")
        track = create_track_from_segments(track_segments)
        
        print(f"Track created with {track.NumberTrajectoryPoints()} trajectory points")
        print(f"Track length: {track.Length():.1f} cm")
        
        # Step 2: Create realistic ADC images
        print("\n2. Creating ADC images...")
        adc_images = create_realistic_adc_images()
        
        # Step 3: Add charge deposits along track
        print("\n3. Adding charge deposits...")
        add_track_charge_to_images(adc_images, track, args.charge_per_cm)
    
    # Convert track to space points
    print("\nConverting track to space points...")
    converter = ublarcvapp.pixelutils.TrackToSpacePoints()
    converter.setUseChargeWeighting(False)
    
    spacepoints = converter.convertTrack(
        track, adc_images,
        threshold=10.0,  # ADC threshold
        dcol=3,         # wire window  
        drow=3,         # tick window
        minstepsize=0.3, # min step (cm)
        maxstepsize=0.5  # max step (cm)
    )
    
    print(f"Generated {len(spacepoints)} space points")
    print(f"Pixels processed: {converter.getNumPixelsProcessed()}")
    print(f"Total charge collected: {converter.getTotalCharge():.0f} ADC")
    
    # Convert to photon sources
    print("\nConverting charge to photons...")
    photon_estimator = ublarcvapp.ubphotonlib.PhotonVisibilityEstimator()
    
    total_photons_added = 0
    for i, sp in enumerate(spacepoints):
        if sp.charge > 5.0:  # Only use points with significant charge
            n_photons = convert_charge_to_photons(sp.charge)
            
            photon_estimator.addPhotonSource(
                sp.position.X(),
                sp.position.Y(), 
                sp.position.Z(),
                n_photons
            )
            total_photons_added += n_photons
            
            if i < 5:  # Print first few
                print(f"  Point {i}: ({sp.position.X():.1f}, {sp.position.Y():.1f}, {sp.position.Z():.1f}) "
                      f"charge={sp.charge:.1f} -> {n_photons:.0f} photons")
    
    print(f"Added {photon_estimator.getNumSources()} photon sources")
    print(f"Total emitted photons: {total_photons_added:.0f}")
    
    # Estimate light detection  
    print("\nEstimating light detection...")
    
    # Calculate with trilinear interpolation
    photons_per_pmt = photon_estimator.calculateDetectedPhotons(use_trilinear=True)
    total_detected = photon_estimator.getTotalDetectedPhotons(True)
    efficiency = photon_estimator.getCollectionEfficiency(True)
    
    print(f"Total detected photons: {total_detected:.1f}")
    print(f"Collection efficiency: {efficiency:.4f} ({efficiency*100:.2f}%)")
    
    # Show detection by PMT
    print("\nPhotons detected by PMT:")
    for pmt in range(32):
        pe = photons_per_pmt[pmt]
        if pe > 0.5:
            print(f"  PMT {pmt:2d}: {pe:7.1f} PE")
    
    # Create visualization
    if args.save_output:
        print("\nCreating visualization...")
        
        # PMT response histogram (predicted)
        h_pmt_pred = rt.TH1F("h_pmt_pred", "PMT Light Response;PMT ID;Normalized Response", 
                             32, 0, 32)
        for pmt in range(32):
            h_pmt_pred.SetBinContent(pmt+1, photons_per_pmt[pmt])
        
        # PMT response histogram (observed) - only if we have opflash data
        h_pmt_obs = None
        if opflashes is not None and len(opflashes) > 0:
            flash = opflashes[0]  # Use first flash
            h_pmt_obs = rt.TH1F("h_pmt_obs", "PMT Light Response;PMT ID;Normalized Response", 
                                32, 0, 32)
            
            # Fill observed PMT response and normalize to 1
            total_obs_pe = 0
            pe_per_pmt = []
            for pmt in range(32):
                pe = flash.PE(pmt)
                pe_per_pmt.append(pe)
                total_obs_pe += pe
            
            # Normalize observed data so total = 1 for shape comparison
            if total_obs_pe > 0:
                for pmt in range(32):
                    normalized_pe = pe_per_pmt[pmt] / total_obs_pe
                    h_pmt_obs.SetBinContent(pmt+1, normalized_pe)
            
            print(f"\nObserved flash: Total PE = {total_obs_pe:.1f}")
            print("Observed PE per PMT (>0.1):")
            for pmt in range(32):
                if pe_per_pmt[pmt] > 0.1:
                    print(f"  PMT {pmt:2d}: {pe_per_pmt[pmt]:6.1f} PE")
        
        # Normalize predicted data to 1 for shape comparison  
        if total_detected > 0:
            for pmt in range(32):
                normalized_pred = photons_per_pmt[pmt] / total_detected
                h_pmt_pred.SetBinContent(pmt+1, normalized_pred)
        
        # Track 3D visualization
        h_track_xz = rt.TH2F("h_track_xz", "Track Projection X-Z;X [cm];Z [cm]", 
                             100, 0, 200, 100, 100, 800)
        h_track_yz = rt.TH2F("h_track_yz", "Track Projection Y-Z;Y [cm];Z [cm]", 
                             100, -50, 150, 100, 100, 800)
        
        # Fill track histograms
        for i in range(track.NumberTrajectoryPoints()):
            pos = track.LocationAtPoint(i)
            h_track_xz.Fill(pos.X(), pos.Z())
            h_track_yz.Fill(pos.Y(), pos.Z())
        
        # Space points with charge
        h_charge_xz = rt.TH2F("h_charge_xz", "Space Points with Charge;X [cm];Z [cm]", 
                              100, 0, 200, 100, 100, 800)
        for sp in spacepoints:
            if sp.charge > 0:
                h_charge_xz.Fill(sp.position.X(), sp.position.Z(), sp.charge)
        
        # Create canvas with space for comparison plot
        canvas = rt.TCanvas("c1", "Track Light Estimation", 1600, 1000)
        canvas.Divide(3, 3)
        
        canvas.cd(1)
        h_track_xz.SetMarkerStyle(20)
        h_track_xz.SetMarkerColor(rt.kBlue)
        h_track_xz.Draw("P")
        rt.gPad.SetTitle("Track Path (X-Z)")
        
        canvas.cd(2)
        h_track_yz.SetMarkerStyle(20)
        h_track_yz.SetMarkerColor(rt.kBlue)
        h_track_yz.Draw("P")
        rt.gPad.SetTitle("Track Path (Y-Z)")
        
        canvas.cd(3)
        h_charge_xz.SetMarkerStyle(20)
        h_charge_xz.SetMarkerColor(rt.kRed)
        h_charge_xz.Draw("COLZ")
        rt.gPad.SetTitle("Charge Deposits")
        
        # PMT comparison plot
        canvas.cd(4)
        h_pmt_pred.SetLineColor(rt.kRed)
        h_pmt_pred.SetLineWidth(2)
        h_pmt_pred.SetTitle("PMT Response Comparison;PMT ID;Normalized Response")
        max_val = h_pmt_pred.GetMaximum()
        
        if h_pmt_obs is not None:
            h_pmt_obs.SetLineColor(rt.kBlue) 
            h_pmt_obs.SetLineWidth(2)
            max_val = max(max_val, h_pmt_obs.GetMaximum())
            h_pmt_pred.GetYaxis().SetRangeUser(0, max_val * 1.1)
            h_pmt_pred.Draw("HIST")
            h_pmt_obs.Draw("HIST SAME")
            
            # Add legend
            legend = rt.TLegend(0.6, 0.7, 0.89, 0.89)
            legend.AddEntry(h_pmt_pred, "Predicted", "l")
            legend.AddEntry(h_pmt_obs, "Observed", "l")
            legend.Draw()
        else:
            h_pmt_pred.SetFillColor(rt.kCyan-9)
            h_pmt_pred.Draw("HIST")
        
        canvas.cd(5)
        # Individual predicted PMT response (unnormalized)
        h_pmt_pred_raw = rt.TH1F("h_pmt_pred_raw", "Predicted PMT Response;PMT ID;Photoelectrons", 
                                 32, 0, 32)
        for pmt in range(32):
            h_pmt_pred_raw.SetBinContent(pmt+1, photons_per_pmt[pmt])
        h_pmt_pred_raw.SetFillColor(rt.kGreen-9)
        h_pmt_pred_raw.Draw("HIST")
        
        canvas.cd(6)
        if h_pmt_obs is not None:
            # Individual observed PMT response (unnormalized)  
            h_pmt_obs_raw = rt.TH1F("h_pmt_obs_raw", "Observed PMT Response;PMT ID;Photoelectrons", 
                                    32, 0, 32)
            for pmt in range(32):
                h_pmt_obs_raw.SetBinContent(pmt+1, pe_per_pmt[pmt])
            h_pmt_obs_raw.SetFillColor(rt.kYellow-9)
            h_pmt_obs_raw.Draw("HIST")
        
        # Summary text
        canvas.cd(7)
        text = rt.TText()
        text.SetTextSize(0.08)
        text.DrawText(0.1, 0.9, f"Track Length: {track.Length():.1f} cm")
        text.DrawText(0.1, 0.8, f"Space Points: {len(spacepoints)}")
        text.DrawText(0.1, 0.7, f"Total Charge: {converter.getTotalCharge():.0f} ADC")
        text.DrawText(0.1, 0.6, f"Emitted Photons: {total_photons_added:.0f}")
        text.DrawText(0.1, 0.5, f"Detected Photons: {total_detected:.0f}")
        text.DrawText(0.1, 0.4, f"Efficiency: {efficiency*100:.2f}%")
        if opflashes is not None and len(opflashes) > 0:
            text.DrawText(0.1, 0.3, f"Observed PE: {total_obs_pe:.0f}")
            text.DrawText(0.1, 0.2, f"Flash Time: {flash.Time():.1f} us")
        
        canvas.SaveAs("complete_track_light_example.png")
        print("Saved visualization to complete_track_light_example.png")
    
    print("\nExample completed successfully!")

if __name__ == "__main__":
    main()

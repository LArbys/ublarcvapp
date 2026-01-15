#!/usr/bin/env python3
"""
Complete example demonstrating the full workflow from shower to light detection estimation.

This script shows:
1. How to create realistic larlite::shower objects from line segments
2. How to create realistic ADC images with proper geometry and charge deposits
3. How to use showerToSpacePoints to convert showers to space points with charge
4. How to convert charge to photons using LAr physics
5. How to use PhotonVisibilityEstimator to predict light detection

Usage:
    python3 complete_shower_light_example.py [--input-file file.root] [--save-output]
"""

import ROOT as rt
from ROOT import std
import numpy as np
import argparse
import sys
rt.gStyle.SetOptStat(0)

# Load libraries
from larcv import larcv
from larlite import larlite  
from ublarcvapp import ublarcvapp


def get_data_from_rootfile(rootfile):
    """
    Get Image2D and shower examples for testing showerToSpacePoints 
    and the resulting photon prediction using ublarcvapp.ubphotonlib.PhotonVisibilityEstimator.
    The test file is a merged file containing TTrees with larcv objects and
    the reconstructed showers.
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

    # Convert to std::vector for showerToSpacePoints
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

    # Get the larlite::shower examples from the reco file
    reco_tfile = rt.TFile(rootfile, 'open')
    recotree = reco_tfile.Get("KPSRecoManagerTree")
    trackout_v  = std.vector('larlite::track')()
    showerout_v = std.vector('larlite::larflowcluster')()
    
    if not recotree:
        print("Warning: Could not find KPSRecoManagerTree in file")
        return adc_images, showerout_v, opflashes
        
    recotree.GetEntry(0)
    
    # Check if we have the expected structure
    if not hasattr(recotree, 'nuvetoed_v') or recotree.nuvetoed_v.size() == 0:
        print("Warning: No nuvetoed_v found in tree")
        return adc_images, []
        
    nuvtx = recotree.nuvetoed_v.at(0)
    
    print(f"  Found {nuvtx.shower_v.size()} showers")

    # Get showers
    for ishower in range(nuvtx.shower_v.size()):
        shower = nuvtx.shower_v.at(ishower)
        if shower.size()>=2:
            showerout_v.push_back(shower)
            print(f"    Shower {ishower}: {shower.size()} points")

    # Get track
    for itrack in range(nuvtx.track_v.size()):
        track = nuvtx.track_v.at(itrack)
        if track.NumberTrajectoryPoints()>=2:
            trackout_v.push_back(track)
            print(f"    Track {itrack}: {track.NumberTrajectoryPoints()} points, length {track.Length():.1f} cm")

    ioman.finalize()
    reco_tfile.Close()
    
    return adc_images, trackout_v, showerout_v, opflashes


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

def get_shower_opdet_predicts( showers, adc_images ):

    # Loop over shower container
    total_photons_added = 0
    shower_predictions = []
    for ishower in range( showers.size() ):
        shower = showers.at(ishower)

        # make a spacepoint converter for this shower
        converter = ublarcvapp.pixelutils.ShowerToSpacePoints()
        converter.setUseChargeWeighting(False)
        spacepoints = converter.convertShower(
            shower, adc_images,
            threshold=10.0,  # ADC threshold
            dcol=3,         # wire window  
            drow=3          # tick window
        )

        # Make photon estimator
        photon_estimator = ublarcvapp.ubphotonlib.PhotonVisibilityEstimator()

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
        print(f"shower[{ishower}] num photons predicted = {n_photons}")
        photons_per_pmt = photon_estimator.calculateDetectedPhotons(use_trilinear=True)
        shower_predictions.append( photons_per_pmt )

    print(f'Number of total photons: {total_photons_added}')

    return shower_predictions

def get_track_opdet_predicts( tracks, adc_images ):

    # Loop over track container
    total_photons_added = 0
    track_predictions = []
    for itrack in range( tracks.size() ):
        track = tracks.at(itrack)

        # make a spacepoint converter for this track
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

        # Make photon estimator
        photon_estimator = ublarcvapp.ubphotonlib.PhotonVisibilityEstimator()

        for i, sp in enumerate(spacepoints):
            #if sp.charge > 5.0:  # Only use points with significant charge
            n_photons = convert_charge_to_photons(sp.charge)
            
            photon_estimator.addPhotonSource(
                sp.position.X(),
                sp.position.Y(), 
                sp.position.Z(),
                n_photons
            )
            total_photons_added += n_photons
        print(f"track[{itrack}] num photons predicted = {n_photons}")
        photons_per_pmt = photon_estimator.calculateDetectedPhotons(use_trilinear=True)
        track_predictions.append( photons_per_pmt )

    print(f'Number of total photons: {total_photons_added}')

    return track_predictions

def make_showers_visualization( tracks, showers, adc_images, offset=1 ):
    # shower 3D visualization
    h_shower_xz = rt.TH2F("h_shower_xz", "shower Projection X-Z;X [cm];Z [cm]", 
                         100, 0, 256.0, 100, 0, 1036.0)
    h_shower_yz = rt.TH2F("h_shower_yz", "shower Projection Y-Z;Y [cm];Z [cm]", 
                          100, -120, 120.0, 100, 0, 1036.0)
        
    # Fill shower histogram
    for ishower in range( showers.size()):
        shower = showers.at(ishower)
        for i in range(shower.size()):
            hit = shower.at(i)
            h_shower_xz.Fill(hit[0],hit[2], offset+ishower)
            h_shower_yz.Fill(hit[1],hit[2], offset+ishower)
        
    return h_shower_xz, h_shower_yz

def main():
    parser = argparse.ArgumentParser(description="Complete shower to light estimation example")
    parser.add_argument("--input-file", required=True, help="ROOT file with shower/image data (optional)")
    parser.add_argument("--save-output", action="store_true", help="Save output plots")
    parser.add_argument("--charge-per-cm", type=float, default=2000.0, 
                       help="Charge deposition per cm (ADC)")
    args = parser.parse_args()
    
    print("Complete Shower to Light Estimation Example")
    print("="*60)
    
    opflashes = None  # Initialize for synthetic data case
    
    
    # Use real data from ROOT file
    print(f"\nUsing real data from: {args.input_file}")
    adc_images, tracks, showers, opflashes = get_data_from_rootfile( args.input_file )
    
    # Show information about optical flashes
    print(f"\nFound {len(opflashes)} optical flashes")
    if len(opflashes) > 0:
        flash = opflashes[0]  # Use first flash for comparison
        print(f"Using flash 0: Total PE = {flash.TotalPE():.1f}, Time = {flash.Time():.1f} μs")
        
    
    # Convert shower to space points
    shower_pmt_preds = get_shower_opdet_predicts( showers, adc_images)

    # Convert track to space points
    track_pmt_preds = get_track_opdet_predicts( tracks, adc_images)
    
    # Create visualization
    if args.save_output:
        fudge_factor = 1000.0
        print("\nCreating visualization...")
        import lardly
        from lardly.ubdl.pmtpos import getOpChannelFromOpDet, getOpDetFromOpChannel
        
        # PMT response histogram (predicted)
        # We make a histogram for each object and then a stack histogram as well 
        h_pmt_pred_stack = rt.THStack("h_pmt_pred_stack", "Total Neutrino Candidate OpDet Prediction")
        hists_v = []
        colors = []
        iprong = 0
        tlen = rt.TLegend(0.7,0.7,0.9,0.9)
        for itrack, photons_per_pmt in enumerate( track_pmt_preds ):
            h_pmt_pred = rt.TH1F(f"h_pmt_pred_track{itrack}", "PMT Light Response;PMT ID;Normalized Response", 32, 0, 32)

            rand_rgb = np.random.randint(0,255,3)
            rand_tcolor = rt.TColor.GetColor( rand_rgb[0], rand_rgb[1], rand_rgb[2] )
            pesum = 0.0
            for pmt in range(32):
                h_pmt_pred.SetBinContent(pmt+1, photons_per_pmt[pmt]*fudge_factor)
                pesum += photons_per_pmt[pmt]*fudge_factor
            xcolor = int(51+(iprong//5)+10*(iprong%5))
            print(f"xcolor[{iprong}]: {xcolor}")
            h_pmt_pred.SetFillColor( xcolor )
            h_pmt_pred.SetFillStyle(3001)
            h_pmt_pred.SetLineColor( xcolor )
            h_pmt_pred_stack.Add( h_pmt_pred )
            tlen.AddEntry( h_pmt_pred, f"Track[{itrack}]: {pesum:0.2f}" )
            hists_v.append( h_pmt_pred )
            colors.append( xcolor )
            iprong += 1

        for ishower, photons_per_pmt in enumerate( shower_pmt_preds ):
            h_pmt_pred = rt.TH1F(f"h_pmt_pred_shower{ishower}", "PMT Light Response;Optical Channel;Normalized Response", 32, 0, 32)

            rand_rgb = np.random.randint(0,255,3)
            rand_tcolor = rt.TColor.GetColor( rand_rgb[0], rand_rgb[1], rand_rgb[2] )
            pesum = 0.0
            for pmt in range(32):
                # predictions are indexed by opchannel
                opchid = pmt
                #opdetid = getOpDetFromOpChannel(opchid)
                h_pmt_pred.SetBinContent(opchid+1, photons_per_pmt[pmt]*fudge_factor)
                pesum += photons_per_pmt[pmt]*fudge_factor
            xcolor = int(51+(iprong//5)+10*(iprong%5))
            print(f"xcolor[{iprong}]: {xcolor}")
            h_pmt_pred.SetFillColor( xcolor )
            h_pmt_pred.SetFillStyle(3001)
            h_pmt_pred_stack.Add( h_pmt_pred )
            tlen.AddEntry( h_pmt_pred, f"Shower[{ishower}]: {pesum:0.2f}" )
            hists_v.append( h_pmt_pred )
            colors.append( xcolor )
            iprong += 1
        
        # PMT response histogram (observed) - only if we have opflash data
        h_pmt_obs = None
        if opflashes is not None and len(opflashes) > 0:
            flash = opflashes[0]  # Use first flash
            h_pmt_obs = rt.TH1F("h_pmt_obs", "PMT Light Response;PMT ID;Normalized Response", 32, 0, 32)
            h_pmt_obs.SetLineColor( rt.kBlack )
            h_pmt_obs.SetLineWidth( 2 )          

            # Fill observed PMT response and normalize to 1
            total_obs_pe = 0
            pe_per_pmt = []
            for pmt in range(32):
                # opflash is indexed by channel?
                #opdetid = getOpDetFromOpChannel( pmt )
                opchanid = pmt
                pe = flash.PE(opchanid)
                # opflash is indexed by opdet and pmt pred by channel?
                #opchid = getOpChannelFromOpDet( pmt )
                #pe = flash.PE(opchid)                
                pe_per_pmt.append(pe)
                total_obs_pe += pe
            
            # Normalize observed data so total = 1 for shape comparison
            if total_obs_pe > 0:
                for pmt in range(32):
                    normalized_pe = pe_per_pmt[pmt] / total_obs_pe
                    #h_pmt_obs.SetBinContent(pmt+1, normalized_pe)
                    h_pmt_obs.SetBinContent(pmt+1, pe_per_pmt[pmt])
            
            print(f"\nObserved flash: Total PE = {total_obs_pe:.1f}")
            print("Observed PE per PMT (>0.1):")
            for pmt in range(32):
                if pe_per_pmt[pmt] > 0.1:
                    print(f"  PMT {pmt:2d}: {pe_per_pmt[pmt]:6.1f} PE")
    
        
        # # shower 3D visualization
        # h_shower_xz = rt.TH2F("h_shower_xz", "shower Projection X-Z;X [cm];Z [cm]", 
        #                      100, 0, 200, 100, 100, 800)
        # h_shower_yz = rt.TH2F("h_shower_yz", "shower Projection Y-Z;Y [cm];Z [cm]", 
        #                      100, -50, 150, 100, 100, 800)
        
        # # Fill shower histograms
        # for i in range(shower.size()):
        #     hit = shower.at(i)
        #     h_shower_xz.Fill(hit[0],hit[2])
        #     h_shower_yz.Fill(hit[1],hit[2])
        
        # # Space points with charge
        # h_charge_xz = rt.TH2F("h_charge_xz", "Space Points with Charge;X [cm];Z [cm]", 
        #                       100, 0, 200, 100, 100, 800)
        # for sp in spacepoints:
        #     if sp.charge > 0:
        #         h_charge_xz.Fill(sp.position.X(), sp.position.Z(), sp.charge)
        
        # Create canvas with space for comparison plot
        canvas = rt.TCanvas("c1", "shower Light Estimation", 1600, 1000)
        #canvas.Divide(3, 3)
        
        canvas.cd(1)
        pred_max = h_pmt_pred_stack.GetMaximum()
        obs_max  = h_pmt_obs.GetMaximum()

        if pred_max>obs_max:
            h_pmt_pred_stack.Draw("hist")
        else:
            h_pmt_obs.Draw("hist")

        h_pmt_pred_stack.Draw("histsame")
        h_pmt_obs.Draw("E1same")
        tlen.Draw()
        
        # canvas.cd(2)
        # h_shower_yz.SetMarkerStyle(20)
        # h_shower_yz.SetMarkerColor(rt.kBlue)
        # h_shower_yz.Draw("P")
        # rt.gPad.SetTitle("shower Path (Y-Z)")
        
        # canvas.cd(3)
        # h_charge_xz.SetMarkerStyle(20)
        # h_charge_xz.SetMarkerColor(rt.kRed)
        # h_charge_xz.Draw("COLZ")
        # rt.gPad.SetTitle("Charge Deposits")
        
        # # PMT comparison plot
        # canvas.cd(4)
        # h_pmt_pred.SetLineColor(rt.kRed)
        # h_pmt_pred.SetLineWidth(2)
        # h_pmt_pred.SetTitle("PMT Response Comparison;PMT ID;Normalized Response")
        # max_val = h_pmt_pred.GetMaximum()
        
        # if h_pmt_obs is not None:
        #     h_pmt_obs.SetLineColor(rt.kBlue) 
        #     h_pmt_obs.SetLineWidth(2)
        #     max_val = max(max_val, h_pmt_obs.GetMaximum())
        #     h_pmt_pred.GetYaxis().SetRangeUser(0, max_val * 1.1)
        #     h_pmt_pred.Draw("HIST")
        #     h_pmt_obs.Draw("HIST SAME")
            
        #     # Add legend
        #     legend = rt.TLegend(0.6, 0.7, 0.89, 0.89)
        #     legend.AddEntry(h_pmt_pred, "Predicted", "l")
        #     legend.AddEntry(h_pmt_obs, "Observed", "l")
        #     legend.Draw()
        # else:
        #     h_pmt_pred.SetFillColor(rt.kCyan-9)
        #     h_pmt_pred.Draw("HIST")
        
        # canvas.cd(5)
        # # Individual predicted PMT response (unnormalized)
        # h_pmt_pred_raw = rt.TH1F("h_pmt_pred_raw", "Predicted PMT Response;PMT ID;Photoelectrons", 
        #                          32, 0, 32)
        # for pmt in range(32):
        #     h_pmt_pred_raw.SetBinContent(pmt+1, photons_per_pmt[pmt])
        # h_pmt_pred_raw.SetFillColor(rt.kGreen-9)
        # h_pmt_pred_raw.Draw("HIST")
        
        # canvas.cd(6)
        # if h_pmt_obs is not None:
        #     # Individual observed PMT response (unnormalized)  
        #     h_pmt_obs_raw = rt.TH1F("h_pmt_obs_raw", "Observed PMT Response;PMT ID;Photoelectrons", 
        #                             32, 0, 32)
        #     for pmt in range(32):
        #         h_pmt_obs_raw.SetBinContent(pmt+1, pe_per_pmt[pmt])
        #     h_pmt_obs_raw.SetFillColor(rt.kYellow-9)
        #     h_pmt_obs_raw.Draw("HIST")
        
        # # Summary text
        # canvas.cd(7)
        # text = rt.TText()
        # text.SetTextSize(0.08)
        # text.DrawText(0.1, 0.8, f"Space Points: {len(spacepoints)}")
        # text.DrawText(0.1, 0.7, f"Total Charge: {converter.getTotalCharge():.0f} ADC")
        # text.DrawText(0.1, 0.6, f"Emitted Photons: {total_photons_added:.0f}")
        # text.DrawText(0.1, 0.5, f"Detected Photons: {total_detected:.0f}")
        # text.DrawText(0.1, 0.4, f"Efficiency: {efficiency*100:.2f}%")
        # if opflashes is not None and len(opflashes) > 0:
        #     text.DrawText(0.1, 0.3, f"Observed PE: {total_obs_pe:.0f}")
        #     text.DrawText(0.1, 0.2, f"Flash Time: {flash.Time():.1f} us")
        
        canvas.SaveAs("nucand_light_example.png")
        print("Saved visualization to nucand_light_example.png")
        print("[enter] to exit.")
        input()
    
    print("\nExample completed successfully!")

if __name__ == "__main__":
    main()

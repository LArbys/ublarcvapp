/**
 * Example C++ usage of PhotonVisibilityEstimator
 * 
 * Compile with:
 * g++ -o test_photon_estimator test_photon_estimator.cxx \
 *     `root-config --cflags --libs` \
 *     -I$UBLARCVAPP_INCDIR -L$UBLARCVAPP_LIBDIR -lLArCVApp_UBPhotonLib
 */

#include <iostream>
#include <vector>
#include <map>

#include "ublarcvapp/UBPhotonLib/PhotonVisibilityEstimator.h"

using namespace ublarcvapp::ubphotonlib;

int main() {
    
    // Create estimator
    PhotonVisibilityEstimator estimator;
    
    std::cout << "PhotonVisibilityEstimator Example" << std::endl;
    std::cout << "Number of optical channels: " << estimator.getNumOpticalChannels() << std::endl;
    
    // Example 1: Track-like distribution
    std::cout << "\n--- Example 1: Track-like distribution ---" << std::endl;
    
    // Add points along a track
    float x_start = 128.0, y_start = 0.0, z_start = 500.0;  // cm
    float x_end = 128.0, y_end = 50.0, z_end = 600.0;       // cm
    int n_points = 10;
    float photons_per_point = 1000.0;
    
    for (int i = 0; i < n_points; ++i) {
        float t = float(i) / (n_points - 1);
        float x = x_start + t * (x_end - x_start);
        float y = y_start + t * (y_end - y_start);
        float z = z_start + t * (z_end - z_start);
        
        estimator.addPhotonSource(x, y, z, photons_per_point);
    }
    
    std::cout << "Added " << estimator.getNumSources() << " point sources" << std::endl;
    std::cout << "Total emitted photons: " << estimator.getTotalEmittedPhotons() << std::endl;
    
    // Calculate detected photons
    auto photons_per_opdet = estimator.calculateDetectedPhotons(true);
    
    // Print results
    std::cout << "\nPhotons detected by each optical detector:" << std::endl;
    float total_detected = 0.0;
    for (const auto& entry : photons_per_opdet) {
        int opch = entry.first;
        float n_photons = entry.second;
        if (n_photons > 0) {
            std::cout << "  OpDet " << opch << ": " << n_photons << " photons" << std::endl;
            total_detected += n_photons;
        }
    }
    
    std::cout << "\nTotal detected photons: " << total_detected << std::endl;
    std::cout << "Collection efficiency: " << estimator.getCollectionEfficiency(true) << std::endl;
    
    // Example 2: Using PhotonSource struct
    std::cout << "\n--- Example 2: Using PhotonSource struct ---" << std::endl;
    
    estimator.clear();
    
    // Create a vector of sources
    std::vector<PhotonVisibilityEstimator::PhotonSource> sources;
    
    // Add a few discrete sources
    sources.emplace_back(100.0, 0.0, 300.0, 5000.0);
    sources.emplace_back(150.0, 20.0, 400.0, 3000.0);
    sources.emplace_back(200.0, -10.0, 500.0, 4000.0);
    
    // Add all at once
    estimator.addPhotonSources(sources);
    
    std::cout << "Added " << estimator.getNumSources() << " discrete sources" << std::endl;
    std::cout << "Total emitted photons: " << estimator.getTotalEmittedPhotons() << std::endl;
    
    // Compare with and without trilinear interpolation
    float total_trilinear = estimator.getTotalDetectedPhotons(true);
    float total_no_trilinear = estimator.getTotalDetectedPhotons(false);
    
    std::cout << "\nWith trilinear interpolation: " << total_trilinear << " photons detected" << std::endl;
    std::cout << "Without trilinear interpolation: " << total_no_trilinear << " photons detected" << std::endl;
    std::cout << "Difference: " << (total_trilinear - total_no_trilinear) << " photons" << std::endl;
    
    return 0;
}
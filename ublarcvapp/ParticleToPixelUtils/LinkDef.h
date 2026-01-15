/** \defgroup ParticleToPixelUtils
 *
 * \brief Utilities for converting particle representations to pixel/space point data
 *
 * cint script to generate libraries and python bindings.
 * Declare namespace & classes you defined
 * pragma statement: order matters! Google it ;)
 *
 */
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace ublarcvapp;
#pragma link C++ namespace ublarcvapp::pixelutils;

#pragma link C++ class ublarcvapp::pixelutils::TrackToSpacePoints+;
#pragma link C++ class ublarcvapp::pixelutils::TrackToSpacePoints::SpacePointCharge+;
#pragma link C++ class ublarcvapp::pixelutils::TrackToSpacePoints::PixelData+;
#pragma link C++ class ublarcvapp::pixelutils::ShowerToSpacePoints+;
#pragma link C++ class ublarcvapp::pixelutils::ShowerToSpacePoints::SpacePointCharge+;
#pragma link C++ class ublarcvapp::pixelutils::ShowerToSpacePoints::PixelData+;

#endif
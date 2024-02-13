/** \defgroup MCTools
 *
 * \brief Algorithms that take in Monte Carlo (i.e. simulation) truth information
 *
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

#pragma link C++ namespace ublarcvapp+;
#pragma link C++ namespace ublarcvapp::mctools+;
#pragma link C++ class ublarcvapp::mctools::MCPixelPGraph+;
#pragma link C++ class ublarcvapp::mctools::MCPixelPMap+;
#pragma link C++ class ublarcvapp::mctools::CrossingPointsAnaMethods+;
#pragma link C++ class ublarcvapp::mctools::NeutrinoVertex+;
#pragma link C++ class ublarcvapp::mctools::LArbysMC+;
#pragma link C++ class ublarcvapp::mctools::NeutrinoPixelFilter+;
#pragma link C++ class ublarcvapp::mctools::TruthTrackSCE+;
#pragma link C++ class ublarcvapp::mctools::TruthShowerTrunkSCE+;
#pragma link C++ class ublarcvapp::mctools::RecoFlash_t+;
#pragma link C++ class ublarcvapp::mctools::FlashMatcher+;
#pragma link C++ class ublarcvapp::mctools::FlashMatcherV2+;
#endif





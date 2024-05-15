/** \defgroup UBPhotonLib
 *
 * \brief Interface and utilities for utilizing the UBPhotonLibrary
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
#pragma link C++ namespace ublarcvapp::ubphotonlib+;

#pragma link C++ class ublarcvapp::ubphotonlib::UBPhotonLib+;

#endif





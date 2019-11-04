//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace ublarcvapp;
#pragma link C++ namespace ublarcvapp::ncpi0;

// ncpi0
#pragma link C++ class ublarcvapp::ncpi0::Utils+;
#pragma link C++ class ublarcvapp::ncpi0::SaveProbabilities+;
#pragma link C++ class ublarcvapp::ncpi0::SaveCutVariables+;

//ADD_NEW_CLASS ... do not change this line
#endif

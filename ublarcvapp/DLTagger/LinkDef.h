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
#pragma link C++ namespace ublarcvapp::dltagger;

// mrcnnmatch
#pragma link C++ class ublarcvapp::dltagger::MaskMatchData+;
#pragma link C++ class ublarcvapp::dltagger::MaskCombo+;
#pragma link C++ class ublarcvapp::dltagger::CropMaskCombo+;
#pragma link C++ class ublarcvapp::dltagger::FeaturesMaskCombo+;
#pragma link C++ class ublarcvapp::dltagger::Gen3DEndpoints+;
#pragma link C++ class ublarcvapp::dltagger::GenGraphPoints+;
#pragma link C++ class ublarcvapp::dltagger::AStarMaskCombo+;
#pragma link C++ class ublarcvapp::dltagger::MRCNNMatch+;

// dltagger
#pragma link C++ class ublarcvapp::dltagger::DLTagger+;
#pragma link C++ class ublarcvapp::dltagger::DLTaggerProcess+;
#pragma link C++ class ublarcvapp::dltagger::DLTaggerProcessFactory+;

//ADD_NEW_CLASS ... do not change this line
#endif

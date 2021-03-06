/**
 * \file Run3DTracker.h
 *
 * \ingroup Package_Name
 *
 * \brief Class def header for a class Run3DTracker
 *
 * @author hourlier
 */

/** \addtogroup Package_Name

 @{*/
#ifndef __RUN3DTRACKER_H__
#define __RUN3DTRACKER_H__

#include "TH2D.h"
#include "TH1D.h"
#include "TTree.h"
#include "TVector3.h"

#include "DataFormat/storage_manager.h"
#include "DataFormat/track.h"

#include "larcv/core/Processor/ProcessBase.h"
#include "larcv/core/Processor/ProcessFactory.h"
#include "AStarTracker.h"

namespace ublarcvapp {
namespace reco3d {

    /**
     \class ProcessBase
     User defined class Run3DTracker ... these comments are used to generate
     doxygen documentation!
     */
    class Run3DTracker : public larcv::ProcessBase {

    public:

        /// Default constructor
        Run3DTracker(const std::string name="Run3DTracker");

        /// Default destructor
        ~Run3DTracker(){}

        void configure(const larcv::PSet&);
        void initialize();
        bool process(larcv::IOManager& mgr);

        void SetSplineLocation(const std::string& fpath);
        void SetLLOutName(const std::string& foutll) { _foutll = foutll; }
        void advance_larlite();
        void FillMC(const std::vector<larcv::ROI>& mc_roi_v);
        void ClearEvent();
        void ClearVertex();
        void SetOutDir(std::string s){out_dir = s;}
        void MakeTTree();
        void MakeTTree_SCEadded();

        bool IsGoodVertex(int run, int subrun, int event/*, int ROIid*/, int vtxID);
        bool IsGoodEntry(int run, int subrun, int event);
        void ReadVertexFile(std::string filename);
        std::vector<TVector3> GetJarretVertex(int run, int subrun, int event);

        void finalize();

        void MCevaluation();

        private :


        int iTrack;
        larlite::storage_manager _storage;
        ublarcvapp::reco3d::AStarTracker tracker;
        //std::vector< std::vector<int> > _vertexInfo;
        double NeutrinoEnergyTh;
        //double NeutrinoEnergyReco;
        double Ep_t;
        double Em_t;
        double MinLengthAllowed;

        int _run;
        int _subrun;
        int _event;
        int _nentry;
        int _entry;
        int _Nreco;
        int randomSeed;

        bool IsMCC9;

        int _vtx_id;
        int NtracksReco;
        int NtracksReco_sceadded;
        std::vector<std::string> checkEvents;
        std::string _filename;

        TTree *_recoTree;
        TTree *_recoTree_SCEadded;

        std::vector<int>    _trk_id_v;
        std::vector<double> _E_muon_v;
        std::vector<double> _E_proton_v;
        std::vector<double> _Length_v;
        std::vector<double> _Avg_Ion_v;
        std::vector<double> _Avg_IonY_v;
        std::vector<double> _vertexPhi;
        std::vector<double> _vertexTheta;

        std::vector<double> _vertexPhi_2cm;
        std::vector<double> _vertexTheta_2cm;
        std::vector<double> _vertexPhi_5cm;
        std::vector<double> _vertexTheta_5cm;
        std::vector<double> _vertexPhi_7cm;
        std::vector<double> _vertexTheta_7cm;
        std::vector<double> _vertexPhi_10cm;
        std::vector<double> _vertexTheta_10cm;
        std::vector<double> _vertexPhi_12cm;
        std::vector<double> _vertexTheta_12cm;
        std::vector<double> _vertexPhi_15cm;// same as simply _vertexPhi, I just wanted to make the distance averaged on explicit
        std::vector<double> _vertexTheta_15cm;// same as simply _vertexTheta, I just wanted to make the distance averaged on explicit
        std::vector<double> _vertexPhi_17cm;
        std::vector<double> _vertexTheta_17cm;
        std::vector<double> _vertexPhi_20cm;
        std::vector<double> _vertexTheta_20cm;
        std::vector<double> _vertexPhi_30cm;
        std::vector<double> _vertexTheta_30cm;
        std::vector<double> _closestWall;
        std::vector<double> _Ion_5cm_v;
        std::vector<double> _Ion_10cm_v;
        std::vector<double> _Ion_tot_v;


        std::vector<double> recoEndPoints_x;
        std::vector<double> recoEndPoints_y;
        std::vector<double> recoEndPoints_z;

        std::vector<double> _IonY_5cm_v;
        std::vector<double> _IonY_10cm_v;
        std::vector<double> _IonY_tot_v;

        std::vector<double> _Trunc_dQdX1_v;
        std::vector<double> _Trunc_dQdX3_v;
        std::vector<double> _IondivLength_v;

        std::vector< std::vector<double> > _trackQ3_v;
        std::vector< std::vector<double> > _trackQ5_v;
        std::vector< std::vector<double> > _trackQ10_v;
        std::vector< std::vector<double> > _trackQ20_v;
        std::vector< std::vector<double> > _trackQ30_v;
        std::vector< std::vector<double> > _trackQ50_v;
        std::vector< std::vector<double> > _TotalADCvalues_v;
        std::vector< std::vector<double> > _Angle_v;
        std::vector< std::vector<double> > _DeadWireList;
        std::vector<int> _DeadWireList_U;
        std::vector<int> _DeadWireList_V;
        std::vector<int> _DeadWireList_Y;
        std::vector<int>   _Reco_goodness_v;
        std::vector<int>  _track_Goodness_v;

        std::vector<TVector3> MCVertices;
        std::vector<TVector3> recoEndPoints;
        TVector3 MCvertex;
        TVector3 RecoVertex;
        TVector3 _RecoVertex;


        bool GoodVertex;
        bool _missingTrack;
        bool _nothingReconstructed;
        bool _tooShortDeadWire;
        bool _tooShortFaintTrack;
        bool _tooManyTracksAtVertex;
        bool _possibleCosmic;
        bool _possiblyCrossing;
        bool _branchingTracks;
        bool _jumpingTracks;

        float _vtx_x;
        float _vtx_y;
        float _vtx_z;

        double _Ep_t;
        double _Em_t;
        double _Ee_t;

        //
        // start point
        //
        double _MuonStartPoint_X;
        double _ProtonStartPoint_X;
        double _ElectronStartPoint_X;

        double _MuonStartPoint_Y;
        double _ProtonStartPoint_Y;
        double _ElectronStartPoint_Y;

        double _MuonStartPoint_Z;
        double _ProtonStartPoint_Z;
        double _ElectronStartPoint_Z;

        //
        // end point
        //
        double _MuonEndPoint_X;
        double _ProtonEndPoint_X;
        double _ElectronEndPoint_X;

        double _MuonEndPoint_Y;
        double _ProtonEndPoint_Y;
        double _ElectronEndPoint_Y;

        double _MuonEndPoint_Z;
        double _ProtonEndPoint_Z;
        double _ElectronEndPoint_Z;


        //////////////////////////////////////////////
        // values for tracks with added SCE
        //////////////////////////////////////////////

        std::vector<double> _E_muon_v_sceadded;
        std::vector<double> _E_proton_v_sceadded;
        std::vector<double> _Length_v_sceadded;
        std::vector<double> _Avg_Ion_v_sceadded;
        std::vector<double> _Avg_IonY_v_sceadded;
        std::vector<double> _vertexPhi_sceadded;
        std::vector<double> _vertexTheta_sceadded;
        std::vector<double> _vertexPhi_2cm_sceadded;
        std::vector<double> _vertexTheta_2cm_sceadded;
        std::vector<double> _vertexPhi_5cm_sceadded;
        std::vector<double> _vertexTheta_5cm_sceadded;
        std::vector<double> _vertexPhi_7cm_sceadded;
        std::vector<double> _vertexTheta_7cm_sceadded;
        std::vector<double> _vertexPhi_10cm_sceadded;
        std::vector<double> _vertexTheta_10cm_sceadded;
        std::vector<double> _vertexPhi_12cm_sceadded;
        std::vector<double> _vertexTheta_12cm_sceadded;
        std::vector<double> _vertexPhi_15cm_sceadded;
        std::vector<double> _vertexTheta_15cm_sceadded;
        std::vector<double> _vertexPhi_17cm_sceadded;
        std::vector<double> _vertexTheta_17cm_sceadded;
        std::vector<double> _vertexPhi_20cm_sceadded;
        std::vector<double> _vertexTheta_20cm_sceadded;
        std::vector<double> _vertexPhi_30cm_sceadded;
        std::vector<double> _vertexTheta_30cm_sceadded;
        std::vector<double> _closestWall_sceadded;
        std::vector<double> _Ion_5cm_v_sceadded;
        std::vector<double> _Ion_10cm_v_sceadded;
        std::vector<double> _Ion_tot_v_sceadded;
        std::vector<double> recoEndPoints_x_sceadded;
        std::vector<double> recoEndPoints_y_sceadded;
        std::vector<double> recoEndPoints_z_sceadded;
        std::vector<double> _IonY_5cm_v_sceadded;
        std::vector<double> _IonY_10cm_v_sceadded;
        std::vector<double> _IonY_tot_v_sceadded;
        std::vector<double> _Trunc_dQdX1_v_sceadded;
        std::vector<double> _Trunc_dQdX3_v_sceadded;
        std::vector<double> _IondivLength_v_sceadded;
        std::vector<double> trackAvg15cm_x;
        std::vector<double> trackAvg15cm_y;
        std::vector<double> trackAvg15cm_z;
        std::vector<double> trackAvg15cm_x_sceadded;
        std::vector<double> trackAvg15cm_y_sceadded;
        std::vector<double> trackAvg15cm_z_sceadded;

        std::vector< std::vector<double> > _trackQ3_v_sceadded;
        std::vector< std::vector<double> > _trackQ5_v_sceadded;
        std::vector< std::vector<double> > _trackQ10_v_sceadded;
        std::vector< std::vector<double> > _trackQ20_v_sceadded;
        std::vector< std::vector<double> > _trackQ30_v_sceadded;
        std::vector< std::vector<double> > _trackQ50_v_sceadded;
        std::vector< std::vector<double> > _TotalADCvalues_v_sceadded;
        std::vector< std::vector<double> > _Angle_v_sceadded;

        TVector3 RecoVertex_sceadded;
        TVector3 _RecoVertex_sceadded;


        std::vector<TVector3> recoEndPoints_sceadded;
        std::vector<TVector3> trackAvg15cm;
        std::vector<TVector3> trackAvg15cm_sceadded;

        float _vtx_x_sceadded;
        float _vtx_y_sceadded;
        float _vtx_z_sceadded;

        ////////////////////////////////////////////
        ////////////////////////////////////////////
        ////////////////////////////////////////////

        std::string _input_pgraph_producer;
        std::string _img2d_producer;
        std::string _par_pix_producer;
        std::string _true_roi_producer;
        std::string _spline_file;
        std::string _foutll;
        std::string out_dir;
        bool _mask_shower;

	// NAMES OF OUTPUT TREES
	std::string _anatree_name;
	std::string _anatree_sce_name;
	std::string _out_track_prodname;
	std::string _out_vertex_prodname;
	std::string _out_assoc_prodname;
	std::string _out_track_sce_prodname;
	std::string _out_assoc_sce_prodname;
    };

    /**
     \class larcv::Run3DTrackerFactory
     \brief A concrete factory class for larcv::Run3DTracker
     */
    class Run3DTrackerProcessFactory : public larcv::ProcessFactoryBase {
    public:
        /// ctor
        Run3DTrackerProcessFactory() { larcv::ProcessFactory::get().add_factory("Run3DTracker",this); }
        /// dtor
        ~Run3DTrackerProcessFactory() {}
        /// creation method
        larcv::ProcessBase* create(const std::string instance_name) { return new Run3DTracker(instance_name); }

    };

}
}

#endif
/** @} */ // end of doxygen group

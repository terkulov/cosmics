#include <memory>
#include <iostream>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include "TPRegexp.h"
#include <TLorentzVector.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

//////////////trigger info////////////////////////////////////

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTConfigData.h"

#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h" 
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h" 
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"

#include "CondFormats/HcalObjects/interface/HcalRespCorrs.h"
#include "CondFormats/DataRecord/interface/HcalRespCorrsRcd.h"

#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CondFormats/HcalObjects/interface/HcalPedestalWidth.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"

#include "Calibration/IsolatedParticles/interface/CaloPropagateTrack.h"
#include "Calibration/IsolatedParticles/interface/eECALMatrix.h" 
#include "Calibration/IsolatedParticles/interface/eHCALMatrix.h" 

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Calibration/IsolatedParticles/interface/CaloConstants.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"



#define EDM_ML_DEBUG
using namespace std;


class HBHEMuonCosmic2 : public edm::EDAnalyzer {

public:
  explicit HBHEMuonCosmic2(const edm::ParameterSet&);
  ~HBHEMuonCosmic2();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup& );
  virtual void endJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  void   clearVectors();
  int    matchId(const HcalDetId&, const HcalDetId&);
  double activeLength(const DetId&);
  bool   isGoodVertex(const reco::Vertex& vtx);
  double respCorr(const DetId& id);
  double gainFactor(const edm::ESHandle<HcalDbService>&, const HcalDetId& id);
  double pedWidth(const edm::ESHandle<HcalDbService>&, const HcalDetId& id);
  double pedEffective(const edm::ESHandle<HcalDbService>&, const HcalDetId& id);
  double pedestal(const edm::ESHandle<HcalDbService>&, const HcalDetId& id);
  double darkCurrent(const edm::ESHandle<HcalDbService>&, const HcalDetId& id);
  double FCbyPE(const edm::ESHandle<HcalDbService>&, const HcalDetId& id);
  int    depth16HE(int ieta, int iphi);
  spr::propagatedTrackID propagateCALO1(const reco::Track* pTrack, const CaloGeometry* geo,
                                        const MagneticField* bField, bool debug);
  std::pair<bool,HcalDetId> propagateHCALBack1(const reco::Track*, const CaloGeometry*, const MagneticField*, 
                                               bool debug);
  void energyHCALCellmy(HcalDetId, edm::Handle<HFRecHitCollection>&,
              std::vector<std::pair<double,int> >&, double, double, double,
               double, double, double, bool);

  spr::propagatedTrack propagateCalo2(const GlobalPoint &,const GlobalVector &,int charge,
		const MagneticField *,float zdist,float radius,float corner,bool debug ); 		
  GlobalPoint myfind(const reco::Track* pTrack, double depth);
  int myfind_iphi(const reco::Track* pTrack, double depth);
  int myfind_ieta(GlobalPoint p);
  double mycorrection(const reco::Track* pTrack, int ieta, int iphi);
  double mycorrection1(const reco::Track* pTrack, int ieta, int iphi);
//  int myfind(const reco::Track* pTrack, double depth);
  std::pair<double, double> mycheckz(int ieta, double depth);
  std::pair<double, double> mycheckx(int iphi, double depth);
  std::pair<double, double> mychecky(int iphi, double depth);

  // ----------member data ---------------------------
  HLTConfigProvider          hltConfig_;
  edm::Service<TFileService> fs;
  edm::InputTag              HLTriggerResults_, HLTriggerEvent_;
  edm::InputTag              labelHBHERecHit_, labelHFRecHit_, labelEBRecHit_, labelEERecHit_, labelTracks_;
  std::string                labelVtx_, labelMuon_;
  int                        verbosity_, maxDepth_, kount_;
  bool                       useRaw_;
  const int                  MaxDepth=7;
  const HcalTopology        *theHBHETopology;
  HcalRespCorrs             *respCorrs_;

  edm::EDGetTokenT<edm::TriggerResults>                   tok_trigRes_;
  edm::EDGetTokenT<trigger::TriggerEvent>                 tok_trigEvt_;
  edm::EDGetTokenT<reco::VertexCollection>                tok_Vtx_;
  edm::EDGetTokenT<EcalRecHitCollection>                  tok_EB_;
  edm::EDGetTokenT<EcalRecHitCollection>                  tok_EE_;
  edm::EDGetTokenT<HBHERecHitCollection>                  tok_HBHE_;
  edm::EDGetTokenT<HFRecHitCollection>                    tok_HF_;
  edm::EDGetTokenT<reco::MuonCollection>                  tok_Muon_;
  edm::EDGetTokenT<reco::TrackCollection>                 tok_Track_;

//////////////////////////////////////////////////////
  std::vector<HcalDDDRecConstants::HcalActiveLength> actHB, actHE;
  std::vector<std::string>  all_triggers;
//  const std::vector<std::string> trigNames_;
  ////////////////////////////////////////////////////////////
  TFile *out_file;
  TTree *ftree;

  double uhcalo1, uhcalo2, uhcalo3, uhcalo4, uhcalo5, uhcalo6, uhcalo7, weight;
  int accept1, ishot, match, HCal_ieta, HCal_iphi, HCal_ieta_hot, HCal_iphi_hot, myentries, mytight1, mytight2;
  int mymuons1, mymuon1, mymuons2, mymuon2, isize1, isize2;
  double isoR04, isoR03, trkmatch, vx, vy, vz;
  double z_vertex, x_vertex, y_vertex, outereta, mycorr;
  double Zmass, eCal1, eCal3, eCal5, eCal15, eCal25, eHCal1;
  double h1, h2, h3, h4, h5, h6, h7, phi1, phi2, phi3, phi4, phi5, phi6, phi7;



  double charge1, charge2, charge3, charge4, charge5, charge6, charge7;
  double chargec1, chargec2, chargec3, chargec4, chargec5, chargec6, chargec7;
  double chargeg1, chargeg2, chargeg3, chargeg4, chargeg5, chargeg6, chargeg7;
  double chgr1, chgr2, chgr3, chgr4, chgr5, chgr6, chgr7;
//  double chargel1, chargel2, chargel3, chargel4, chargel5, chargel6, chargel7;
  double gain1, gain2, gain3, gain4, gain5, gain6, gain7;
  double rz_depth1, rz_depth2, rz_depth3, rz_depth4; 
  double depth1_start, depth1_end;
  double depth2_start[4], depth2_end[4];
  double depth3_start[5], depth3_end[5];
  double depth4_start[7], depth4_end[7];
  int ndepth2 = 4, ndepth3 = 5, ndepth4 = 7;
  double eta, phi, muon_charge, theta, pt, p, hc1, hc2, hc3, hc4, hc5, hc6, hc7, px, py, pz ;
  double iphi1, iphi2, iphi3, iphi4, iphi5, iphi6, iphi7;
  double ieta1, ieta2, ieta3, ieta4, ieta5, ieta6, ieta7, ieta, iphi, ieta_old;
  double cross1, cross2, cross3, cross4, cross5, cross6, cross7;
  double hcalc1, hcalc2, hcalc3, hcalc4, hcalc5, hcalc6, hcalc7;
  double hcalg1, hcalg2, hcalg3, hcalg4, hcalg5, hcalg6, hcalg7;
  double hcr1, hcr2, hcr3, hcr4, hcr5, hcr6, hcr7;
  double hcal1, hcal2, hcal3, hcal4, hcal5, hcal6, hcal7;
  double resp1, resp2, resp3, resp4, resp5, resp6, resp7;
  int trigger;
  int ifirst, isecond, ithird, iforth;
  int stations, myhits, myhits1;
  int is41, is42, is43, is44, is45, is46, is47;
  int is21, is22, is23, is24;
  int nmuons;
  double cor;
  double xmin1, ymin1, zmin1, xmax1, ymax1, zmax1, x1, y1, z1;
  double xmin2, ymin2, zmin2, xmax2, ymax2, zmax2, x2, y2, z2;
  double xmin3, ymin3, zmin3, xmax3, ymax3, zmax3, x3, y3, z3;
  double xmin4, ymin4, zmin4, xmax4, ymax4, zmax4, x4, y4, z4;
  double rmin, r0;

  TTree                    *tree_;
  std::vector<bool>         muon_is_good_, muon_global_, muon_tracker_;// ismediummuon_; //isloosemuon_, istightmuon_;
  std::vector<int>          hltresults;
  unsigned int              runNumber_, nvtx_, nvtx_notFake_, nvtx_good_,eventNumber_ , lumiNumber_, bxNumber_;
 };

HBHEMuonCosmic2::HBHEMuonCosmic2(const edm::ParameterSet& iConfig): theHBHETopology(nullptr), respCorrs_(nullptr) {
  //now do what ever initialization is needed
  kount_            = 0;
  HLTriggerResults_ = iConfig.getParameter<edm::InputTag>("HLTriggerResults");
  HLTriggerEvent_ = iConfig.getParameter<edm::InputTag>("HLTriggerEvent");
  labelVtx_         = iConfig.getParameter<std::string>("LabelVertex");
  //labelEBRecHit_    = iConfig.getParameter<std::string>("LabelEBRecHit");
  //labelEERecHit_    = iConfig.getParameter<std::string>("LabelEERecHit");
  //labelHBHERecHit_  = iConfig.getParameter<std::string>("LabelHBHERecHit");
  labelEBRecHit_    = iConfig.getParameter<edm::InputTag>("LabelEBRecHit");
  labelEERecHit_    = iConfig.getParameter<edm::InputTag>("LabelEERecHit");
  labelHBHERecHit_  = iConfig.getParameter<edm::InputTag>("LabelHBHERecHit");
  labelHFRecHit_    = iConfig.getParameter<edm::InputTag>("LabelHFRecHit");
  labelMuon_        = iConfig.getParameter<std::string>("LabelMuon");
  labelTracks_      = iConfig.getParameter<edm::InputTag>("LabelTracks");
  verbosity_        = iConfig.getUntrackedParameter<int>("Verbosity",0);
  maxDepth_         = iConfig.getUntrackedParameter<int>("MaxDepth",7);
  if (maxDepth_ > MaxDepth) maxDepth_ = MaxDepth;
  else if (maxDepth_ < 1)   maxDepth_ = 4;
  std::string modnam = iConfig.getUntrackedParameter<std::string>("ModuleName","");
  std::string procnm = iConfig.getUntrackedParameter<std::string>("ProcessName","");
  useRaw_            = iConfig.getParameter<bool>("UseRaw");

  tok_trigRes_  = consumes<edm::TriggerResults>(HLTriggerResults_);
  tok_trigEvt_  = consumes<trigger::TriggerEvent>(HLTriggerEvent_);
  tok_EB_       = consumes<EcalRecHitCollection>(labelEBRecHit_);
  tok_EE_       = consumes<EcalRecHitCollection>(labelEERecHit_);
  tok_HBHE_     = consumes<HBHERecHitCollection>(labelHBHERecHit_);
  tok_HF_       = consumes<HFRecHitCollection>(labelHFRecHit_);
  tok_Track_     = consumes<reco::TrackCollection>(labelTracks_);

  if (modnam == "") {
    tok_Vtx_      = consumes<reco::VertexCollection>(labelVtx_);
    tok_Muon_     = consumes<reco::MuonCollection>(labelMuon_);
    edm::LogVerbatim("HBHEMuon")  << "Labels used " << HLTriggerResults_ << " "
                                  << labelVtx_ << " " << labelEBRecHit_ << " "
                                  << labelEERecHit_ << " " << labelHBHERecHit_
                                  << " " << labelMuon_;
  } else {
    tok_Vtx_      = consumes<reco::VertexCollection>(edm::InputTag(modnam,labelVtx_,procnm));
    tok_Muon_     = consumes<reco::MuonCollection>(edm::InputTag(modnam,labelMuon_,procnm));
    edm::LogVerbatim("HBHEMuon")   << "Labels used "   << HLTriggerResults_
                                   << "\n            " << edm::InputTag(modnam,labelVtx_,procnm)
                                   << "\n            " << labelEBRecHit_
                                   << "\n            " << labelEERecHit_
                                   << "\n            " << labelHBHERecHit_
                                   << "\n            " << edm::InputTag(modnam,labelMuon_,procnm);
  }
//  cout<<"Muob Collection "<<labelMuon_<<endl;

}

HBHEMuonCosmic2::~HBHEMuonCosmic2() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void HBHEMuonCosmic2::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  ++kount_;
  z_vertex = 1000;
  cout<<"-------------- analyse "<<kount_<<endl;
  clearVectors();
  // depthHE is the first depth index for HE for |ieta| = 16
  // It used to be 3 for all runs preceding 2017 and 4 beyond that
  int depthHE = (maxDepth_ <= 6) ? 3 : 4;
//  cout<<"depthHE "<<depthHE<<endl;
  runNumber_   = iEvent.id().run();
  eventNumber_ = iEvent.id().event();
  lumiNumber_  = iEvent.id().luminosityBlock();
  bxNumber_    = iEvent.bunchCrossing();
//#ifdef EDM_ML_DEBUG
// cout<<"************* define "<<useRaw_<<endl;
//#endif



#ifdef EDM_ML_DEBUG
//  cout<<"EDM_ML_DEBUG"<<endl;
  edm::LogInfo("HBHEMuon") << "Run " << runNumber_ << " Event " << eventNumber_
			   << " Lumi " << lumiNumber_ << " BX " << bxNumber_
			   << std::endl;
#endif  
  edm::Handle<edm::TriggerResults> _Triggers;
  iEvent.getByToken(tok_trigRes_, _Triggers); 
#ifdef EDM_ML_DEBUG
  edm::LogInfo("HBHEMuon") << "Size of all triggers "  
			   << all_triggers.size() << std::endl;
#endif
  int Ntriggers = all_triggers.size();
#ifdef EDM_ML_DEBUG
  edm::LogInfo("HBHEMuon") << "Size of HLT MENU: " << _Triggers->size()
			   << std::endl;
#endif
  if (_Triggers.isValid()) {
    const edm::TriggerNames &triggerNames_ = iEvent.triggerNames(*_Triggers);
    std::vector<int> index;
    for (int i=0; i<Ntriggers; i++) {
      index.push_back(triggerNames_.triggerIndex(all_triggers[i]));
      int triggerSize = int( _Triggers->size());
#ifdef EDM_ML_DEBUG
      edm::LogInfo("HBHEMuon") << "outside loop " << index[i]
			       << "\ntriggerSize " << triggerSize
			       << std::endl;
#endif
      if (index[i] < triggerSize) {
	hltresults.push_back(_Triggers->accept(index[i]));
#ifdef EDM_ML_DEBUG
	edm::LogInfo("HBHEMuon") << "Trigger_info " << triggerSize
				 << " triggerSize " << index[i]
				 << " trigger_index " << hltresults.at(i)
				 << " hltresult" << std::endl;
#endif
      } else {
	edm::LogInfo("HBHEMuon") << "Requested HLT path \"" 
				 << "\" does not exist\n";
      }
    }
  }
  #ifdef EDM_ML_DEBUG
  edm::LogInfo("HBHEMuon") << "after trigger ------------ " << std::endl;
  #endif

  edm::ESHandle<CaloGeometry> pG;
  iSetup.get<CaloGeometryRecord>().get(pG);
  const CaloGeometry* geo = pG.product();
  
  edm::ESHandle<MagneticField> bFieldH;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldH);
  const MagneticField* bField = bFieldH.product();
//  MagneticField* bField;
  cout<<"magnet "<<bField->nominalValue()<<endl;
  
  
  edm::ESHandle<EcalChannelStatus> ecalChStatus;
  iSetup.get<EcalChannelStatusRcd>().get(ecalChStatus);
  const EcalChannelStatus* theEcalChStatus = ecalChStatus.product();
  
  edm::ESHandle<EcalSeverityLevelAlgo> sevlv;
  iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevlv);
  
  edm::ESHandle<CaloTopology> theCaloTopology;
  iSetup.get<CaloTopologyRecord>().get(theCaloTopology);
  const CaloTopology *caloTopology = theCaloTopology.product();
 
  /* 
  edm::ESHandle<HcalTopology> htopo;
  iSetup.get<HcalRecNumberingRecord>().get(htopo);
  const HcalTopology* theHBHETopology = htopo.product();
  */

  edm::ESHandle<HcalDbService> conditions;
  iSetup.get<HcalDbRecord>().get(conditions);

  // Relevant blocks from iEvent
  edm::Handle<reco::VertexCollection> vtx;
  iEvent.getByToken(tok_Vtx_, vtx);
  
  edm::Handle<EcalRecHitCollection> barrelRecHitsHandle;
  iEvent.getByToken(tok_EB_, barrelRecHitsHandle);
  edm::Handle<EcalRecHitCollection> endcapRecHitsHandle;
  iEvent.getByToken(tok_EE_, endcapRecHitsHandle);
  
  edm::Handle<HBHERecHitCollection> hbhe;
  iEvent.getByToken(tok_HBHE_, hbhe);

  edm::Handle<HFRecHitCollection> hf;
  iEvent.getByToken(tok_HF_, hf);
//  cout<<"HF valid "<<hf.isValid()<<endl;

  edm::Handle<reco::TrackCollection> _Track;
  iEvent.getByToken(tok_Track_, _Track);
//  cout<<"track valid "<<_Track.isValid()<<endl;
  
  
  edm::Handle<reco::MuonCollection> _Muon;
  iEvent.getByToken(tok_Muon_, _Muon);

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(tok_trigRes_, triggerResults);

  trigger::TriggerEvent triggerEvent;
  edm::Handle<trigger::TriggerEvent> triggerEventHandle;
  iEvent.getByToken(tok_trigEvt_, triggerEventHandle);
//  cout<<"triggerevent valid "<<triggerEventHandle.isValid()<<endl;

//  cout<<"muon valid "<<_Muon.isValid()<<endl;
  // get handles to calogeometry and calotopology
  #ifdef EDM_ML_DEBUG
  edm::LogInfo("HBHEMuon") << "before checking vtx.isValid() ------------ " << std::endl;
  #endif
//  cout<<" before vertex "<<endl;
  if (!(vtx.isValid())) {                 return;}
  reco::VertexCollection::const_iterator firstGoodVertex = vtx->end();
  for (reco::VertexCollection::const_iterator it = vtx->begin(); it != firstGoodVertex; it++)
  {
//   cout<<"vertex iterator "<<it->ndof()<<" "<<it->isFake()<<" "<<it->position().Z()<<" "<<it->position().Rho()<<endl;
    if (isGoodVertex(*it)) {
//      cout<<"good "<<endl;
      firstGoodVertex = it;
      break;
    }
  }
// my test part
  trigger = 0;
  if (triggerResults.isValid()) 
  {
   const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults);
   const std::vector<std::string>& names = triggerNames.triggerNames();
   for (unsigned int iHLT = 0; iHLT < triggerResults->size(); iHLT++) 
   {
    int hlt = triggerResults->accept(iHLT);
    if (hlt > 0)
    {
     cout<<"This trigger "<<names[iHLT]<<" Flag "<<hlt<<" "<<iHLT<<" "<<triggerResults->size()<<endl;
     if ( names[iHLT] == "DST_Physics_v1" )
     {
//      cout<<"trigger 25"<<endl;
      trigger = 1;
     }
//     if ( names[iHLT] == "HLT_L1SingleMuCosmics_v1" ) cout<<"good "<<endl;
    }
   }
  }

  if (triggerEventHandle.isValid()) 
  {
//   cout<<"size "<<triggerEventHandle->size()<<endl;
   triggerEvent = *(triggerEventHandle.product());
//   cout<<"name "<<triggerEvent.usedProcessName()<<" "<<triggerEvent.sizeFilters()<<endl;
   for (unsigned int ifilter = 0; ifilter < triggerEvent.sizeFilters(); ++ifilter) 
   {
    std::vector<int> Keys;
//    cout<<"filter "<<triggerEvent.filterTag(ifilter).label()<<endl;
   }
  }


  if (hf.isValid())
  {
   for (HFRecHitCollection::const_iterator j=hf->begin(); j != hf->end(); j++)
   {
    HcalDetId cell(j->id());
//    if ( (cell.ieta() == 30 || cell.ieta() == 29 ) && cell.iphi() == 1 )
//    if ( fabs(cell.ieta()) == 29  )
//     cout<<"hfcoll eta depth phi ene "<<cell.ieta()<<" "<<cell.depth()<<" "<<cell.iphi()<<" "<<j->energy()<<endl;
//     cout<<"HF valid "<<endl;
   }
  }
  if (hbhe.isValid())
  {
   for (HBHERecHitCollection::const_iterator j=hbhe->begin(); j != hbhe->end(); j++)
   {
    HcalDetId cell(j->id());
//    if ( fabs(cell.ieta()) < 17  )
//     cout<<"hbhe eta/depth/phi/ene/eraw "<<cell.ieta()<<" "<<cell.depth()<<" "<<cell.iphi()<<" "<<j->energy()<<"/"<<j->eraw()<<endl;
   }
  } 


  #ifdef EDM_ML_DEBUG
  edm::LogInfo("HBHEMuon") << "before vertex: ------------ " << std::endl;
  #endif

  // require a good vertex
  //my comment comment next line for MC with noPU
//  if (firstGoodVertex == vtx->end()) {  return;}
  #ifdef EDM_ML_DEBUG
  edm::LogInfo("HBHEMuon") << "after vertex: ------------ " << std::endl;
  #endif

   nvtx_=0; nvtx_notFake_=0; nvtx_good_=0;
   for (reco::VertexCollection::const_iterator it = vtx->begin(); it != vtx->end(); it++) {
    nvtx_++;
    if (!it->isFake()) nvtx_notFake_++;
    if (isGoodVertex(*it)) { nvtx_good_++; z_vertex = it->z(); x_vertex = it->x(); y_vertex = it->y();}
//  my comment uncomment next line for MC with noPU
    z_vertex = it->z(); x_vertex = it->x(); y_vertex = it->y();
//    cout<<"good vertex "<<nvtx_good_<<" "<<it->z()<<" "<<it->position().Z()<<" "<<it->x()<<" "<<it->y()<<endl;
   }
//  cout<<"good vertex "<<nvtx_good_<<endl;
//  bool accept(false);
  bool accept2(false);
  
  #ifdef EDM_ML_DEBUG
  //edm::LogInfo("HBHEMuon") << "accept: " << accept<<std::endl;
  #endif
//  cout << "check before muon loop"<<std::endl;
//  cout << "hbhe valid/barrel valid/endcap valid/muon valid "<<hbhe.isValid()<<"/"<<barrelRecHitsHandle.isValid()
//  <<"/"<<endcapRecHitsHandle.isValid()<<"/"<<_Muon.isValid()<<endl;
  mymuons1 = 0;
  mymuons2 = 0;
  if (_Muon.isValid() && barrelRecHitsHandle.isValid() && 
      endcapRecHitsHandle.isValid() && hbhe.isValid()) 
  { 
//    cout << "before loop into muons size "<<_Muon->size()<<std::endl;
    for (reco::MuonCollection::const_iterator RecMuon = _Muon->begin(); RecMuon!= _Muon->end(); ++RecMuon)  
    {
     nmuons = _Muon->size();
     int i_start1, i_start2, i_start3, i_start4;
     isoR04 = (RecMuon->pfIsolationR04().sumChargedHadronPt + std::max(0.,RecMuon->pfIsolationR04().sumNeutralHadronEt + RecMuon->pfIsolationR04().sumPhotonEt - (0.5 *RecMuon->pfIsolationR04().sumPUPt))) / RecMuon->pt();
     stations = RecMuon->numberOfMatchedStations();
//     for ( int myjj = -1; myjj < 2; myjj += 2 )
     for ( int myjj = -1; myjj < 2; myjj += 2 )
     {
      cout<<"---------- new track"<<endl;
      ifirst = isecond = ithird = iforth = 0;
      i_start1 = 0;
      i_start2 = i_start3 = i_start4 = 0;
      ieta = iphi = -1000;
      int isHot(0);
      AlgebraicSymMatrix55 mycov;
      mycov(0,0) = 1.;
      mycov(1,1) = 1.;
      mycov(2,2) = 1.;
      mycov(3,3) = 1.;
      mycov(4,4) = 1.;
      
//      if ( RecMuon->py() > 0 ) myjj = -1;
//      else myjj = 1;
      const math::XYZVector myv(myjj*RecMuon->px(),myjj*RecMuon->py(),myjj*RecMuon->pz());
      const math::XYZPoint myp(RecMuon->vx(),RecMuon->vy(),RecMuon->vz());
      r0 = sqrt(RecMuon->vx()*RecMuon->vx()+RecMuon->vy()*RecMuon->vy());
//      const math::XYZVector myv(1,20,1);
//      const math::XYZPoint myp(0,0,0);
      const reco::Track pTrack(1., 2., myp, myv, 1, mycov);
      p = pTrack.p();
      px = pTrack.px();
      py = pTrack.py();
      pz = pTrack.pz();
      eta = pTrack.eta();
      theta = pTrack.theta();
      phi = pTrack.phi();
      pt = pTrack.pt();
      vx = pTrack.vx();
      vy = pTrack.vy();
      vz = pTrack.vz();
//      const reco::Track* pTrack2 = (RecMuon->outerTrack()).get();
      cout<<"muon p/eta/theta/phi/px/py/pz/vx/vy/vz "<<RecMuon->p()<<"/"<<RecMuon->eta()
      <<"/"<<RecMuon->theta()*180/TMath::Pi()<<"/"<<RecMuon->phi()*180/TMath::Pi()
      <<"/"<<RecMuon->px()<<"/"<<RecMuon->py()<<"/"<<RecMuon->pz()<<"/"<<RecMuon->vx()
      <<"/"<<RecMuon->vy()<<"/"<<RecMuon->vz()<<endl;
      cout<<"track reco eta/theta/phi/px/py/pz "<<pTrack.eta()<<"/"<<pTrack.theta()*180/TMath::Pi()<<"/"
      <<pTrack.phi()*180/TMath::Pi()<<"/"<<pTrack.px()<<"/"<<pTrack.py()<<"/"<<pTrack.pz()<<endl;
//      cout<<"track2 reco eta/theta/phi/px/py/pz "<<pTrack2->eta()<<"/"<<pTrack2->theta()*180/TMath::Pi()<<"/"
//      <<pTrack2->phi()*180/TMath::Pi()<<"/"<<pTrack2->px()<<"/"<<pTrack2->py()<<"/"<<pTrack2->pz()<<endl;
      GlobalPoint mypoint;
      i_start4 = 0;
      step4:
      hcal1=hcal2=hcal3=hcal4=hcal5=hcal6=hcal7=-10000;
      hcalc1=hcalc2=hcalc3=hcalc4=hcalc5=hcalc6=hcalc7=-10000;
      gain1=gain2=gain3=gain4=gain5=gain6=gain7=-10000;
      iphi1=iphi2=iphi3=iphi4=iphi5=iphi6=iphi7=-10000;
      ieta1=ieta2=ieta3=ieta4=ieta5=ieta6=ieta7=-10000;
      charge1=charge2=charge3=charge4=charge5=charge6=charge7=-10000;
      chargec1=chargec2=chargec3=chargec4=chargec5=chargec6=chargec7=-10000;
      hcalg1=hcalg2=hcalg3=hcalg4=hcalg5=hcalg6=hcalg7=-10000;
      hcr1=hcr2=hcr3=hcr4=hcr5=hcr6=hcr7=-10000;
      chargeg1=chargeg2=chargeg3=chargeg4=chargeg5=chargeg6=chargeg7=-10000;
      chgr1=chgr2=chgr3=chgr4=chgr5=chgr6=chgr7=-10000;
      match = 0;
      cout<<"depth4 istart4 "<<i_start4<<"/"<<ndepth4-i_start4-1<<endl;
//      mypoint = myfind(&pTrack, depth1_start);
      mypoint = myfind(&pTrack, depth4_start[ndepth4-i_start4-1]);
      int iphi_my = myfind_iphi(&pTrack, depth4_start[ndepth4-i_start4-1]);
      int ieta_my = myfind_ieta(mypoint);
//      cout<<"myieta/myiphi "<<ieta_my<<"/"<<iphi_my<<endl;

//m      cout<<"mypoint depth1_end "<<mypoint.x()<<"/"<<mypoint.y()<<"/"<<mypoint.z()<<"/"<<depth1_end<<endl;

      if (RecMuon->outerTrack().isNonnull()) 
       cout<<"outer p/theta/phi "<<RecMuon->outerTrack()->p()<<"/"<<RecMuon->outerTrack()->theta()
       <<"/"<<RecMuon->outerTrack()->phi()<<endl;
      if (RecMuon->innerTrack().isNonnull()) 
       cout<<"inner p/theta/phi "<<RecMuon->innerTrack()->p()<<"/"<<RecMuon->innerTrack()->theta()
       <<"/"<<RecMuon->innerTrack()->phi()<<endl;
      if (RecMuon->globalTrack().isNonnull()) 
       cout<<"global p/theta/phi "<<RecMuon->globalTrack()->p()<<"/"<<RecMuon->globalTrack()->theta()
       <<"/"<<RecMuon->globalTrack()->phi()<<endl;

      if (RecMuon->outerTrack().isNonnull()) 
      {
       const CaloSubdetectorGeometry* mygHB = geo->getSubdetectorGeometry(DetId::Hcal,HcalBarrel);
//m       cout<<"found x/y/z/eta/phi "<<mypoint.x()<<"/"<<mypoint.y()<<"/"<<mypoint.z()
//m       <<"/"<<sqrt(mypoint.x()*mypoint.x()+mypoint.y()*mypoint.y())
//m       <<"/"<<mypoint.eta()<<"/"<<mypoint.phi()<<endl;           
       spr::propagatedTrackID trackID = propagateCALO1(&pTrack, geo, bField, false);
       const DetId closestCell2(trackID.detIdHCAL);
       cout<<"from propagateCALO1 ieta/iphi "<<HcalDetId(closestCell2).ieta()<<"/"<<HcalDetId(closestCell2).iphi()
           <<endl;
       const DetId closestCell(mygHB->getClosestCell(mypoint));
       GlobalPoint myp1 = myfind(&pTrack, depth4_end[ndepth4-i_start4-1]);
       const DetId closestCell1(mygHB->getClosestCell(myp1));
       if ( HcalDetId(closestCell).ieta() == HcalDetId(closestCell1).ieta() &&
            HcalDetId(closestCell).iphi() == HcalDetId(closestCell1).iphi() )
       {
        myhits1 = 1;
        cout<<"0000000000000000 depth 4 crossed ieta/ieta1/iphi/iphi1 "<<HcalDetId(closestCell).ieta()<<"/"<<
        HcalDetId(closestCell1).ieta()<<"/"<<HcalDetId(closestCell).iphi()<<"/"<<
        HcalDetId(closestCell1).iphi()<<endl;
       }       
       HcalDetId hcidt(closestCell.rawId());  
//m       if ( fabs(hcidt.ieta()) <= 30)
//m        cout<<"init ieta/iphi "<<hcidt.ieta()<<"/"<<hcidt.iphi()<<endl;

       HcalDetId check;
//	std::pair<bool,HcalDetId> info = spr::propagateHCALBack(pTrack,  geo, bField, (((verbosity_/100)%10>0)));
       std::pair<bool,HcalDetId> info = propagateHCALBack1(&pTrack,  geo, bField, false);
       if (info.first) { 
	check = info.second;
       }	
//m       if ( fabs(hcidt.ieta()) <= 30)
//        cout<<"check ieta1/ieta2/iphi1/iphi2 "<<hcidt.ieta()<<"/"<<check.ieta()
//        <<"/"<<hcidt.iphi()<<"/"<<check.iphi()<<endl;
       
//       if (myvdet.okHCAL) 
       if (closestCell) 
       {
	if ((hcidt.ieta() == check.ieta()) && (hcidt.iphi() == check.iphi()))
        {
         cout<<"check1 ieta/ietac/iphi/iphic "<<hcidt.ieta()<<"/"<<check.ieta()<<"/"
         <<hcidt.iphi()<<"/"<<check.iphi()<<endl;
	 match= 1;
        }
        else match = 0;
	HcalSubdetector subdet = HcalDetId(closestCell).subdet();
	ieta   = HcalDetId(closestCell).ieta();
	iphi   = HcalDetId(closestCell).iphi();
        HcalDetId           hotCell;
        spr::eHCALmatrix(geo, theHBHETopology, closestCell, hbhe, 1,1, hotCell, false, useRaw_, false);
        int isHot = matchId(closestCell,hotCell);
        int ietahot = hotCell.ieta();
        int iphihot = hotCell.iphi();
        cout<<"hot cell ieta/ietahot/iphi/iphihot "<<ieta<<"/"<<ietahot<<"/"<<iphi<<"/"<<iphihot<<endl;

        HCal_ieta = ieta;
        HCal_iphi = iphi;
        if ( fabs(ieta) == 1 )
         cout<<"????????????????????????? ieta "<<ieta<<endl;
        std::pair<double, double> par = mycheckz(ieta_my,depth4_start[ndepth4-i_start4-1]);
//m        if ( mypoint.z() < par.first || mypoint.z() > par.second )
        zmin4 = par.first;
        zmax4 = par.second;
        z4 = mypoint.z();
        cout<<"zmin zmax mypoint depth4 "<<par.first<<"/"<<par.second<<"/"<<mypoint.z()<<endl;
        par = mycheckx(iphi_my,depth4_start[ndepth4-i_start4-1]);
        xmin4 = par.first;
        xmax4 = par.second;
        x4 = mypoint.x();
//m        if ( mypoint.x() < par.first || mypoint.x() > par.second )
        cout<<"xmin xmax depth4 "<<par.first<<"/"<<par.second<<"/"<<mypoint.x()<<endl;
        par = mychecky(iphi_my,depth4_start[ndepth4-i_start4-1]);
        ymin4 = par.first;
        ymax4 = par.second;
        y4 = mypoint.y();
//m        if ( mypoint.y() < par.first || mypoint.y() > par.second )
        cout<<"ymin ymax depth4 "<<par.first<<"/"<<par.second<<"/"<<mypoint.y()<<endl;

        cout<<"subdet step 4 ieta/iphi/depth "<<ieta<<"/"<<iphi<<"/"<<HcalDetId(closestCell).depth()<<endl;
        
        cor = mycorrection(&pTrack, ieta, iphi);
        double cor1 = mycorrection1(&pTrack, ieta, iphi);
        cout<<"cor/cor1/phi depth 1 "<<cor<<"/"<<cor1<<"/"<<phi<<endl;
	bool            hborhe = (std::abs(ieta) == 16);
        bool muon_1 = RecMuon->p()>1.0;
//          bool muon_1 = RecMuon->pt()>20.0;
        bool muon_2 = muon::isMediumMuon(*RecMuon);
        bool muon_3 = (RecMuon->pfIsolationR04().sumChargedHadronPt + std::max(0.,RecMuon->pfIsolationR04().sumNeutralHadronEt + RecMuon->pfIsolationR04().sumPhotonEt - (0.5 *RecMuon->pfIsolationR04().sumPUPt))) / RecMuon->pt()<0.15;
        if ( fabs(ieta)< 30 && fabs(ieta) >= 1 && muon_1 )
        {
//           bool okE = trackID.okECAL;
         HcalDetId my_hcal(HcalBarrel, ieta_my, iphi_my, 1);
	 std::vector<std::pair<double,int> > ehdepth;
	 spr::energyHCALCell((HcalDetId) closestCell, hbhe, ehdepth, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depthHE, false);
//	 spr::energyHCALCell(my_hcal, hbhe, ehdepth, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depthHE, false);

         isize1 = ehdepth.size();
         cross1 = 0;
         cross2 = cross3 = cross4 = 0;
         is41 = is42 = is43 = is44 = is45 = is46 = is47 = 0;
         ifirst = isecond = ithird = iforth = 0;
         if ( ehdepth.size() != 0 ) ifirst = 1;
         cout<<"depth 4 "<<isize1<<endl;
         myhits = ehdepth.size();
//         if ( ehdepth.size() == 0 )
//         {
//          if ( i_start4 <= (ndepth4-1) )
//          {
//           i_start4++;
//           goto step4;
//          }
//          else
//           goto step3;
//           goto myend;
//         }

         if ( myhits == 4 )
          cout<<"+++++++++++++++++++ myhits "<<myhits<<endl;
	 for (unsigned int i=0; i<ehdepth.size(); ++i) 
         {
          HcalSubdetector subdet0 = (hborhe) ? ((ehdepth[i].second >= depthHE) ? HcalEndcap : HcalBarrel) : subdet;

          HcalDetId hcid0(subdet0,ieta,iphi,ehdepth[i].second);
          //get charge 
          double ene  = ehdepth[i].first;
          double chg(ene), enec(ene);
          double corr = respCorr(DetId(hcid0));
          if (corr != 0) 
          {
           enec /= corr;
          }
          double gain = gainFactor(conditions,hcid0);
          double pWidth = pedWidth(conditions,hcid0);
          double peds = pedEffective(conditions,hcid0);
          double peds1 = pedestal(conditions,hcid0);
          double dark_current = darkCurrent(conditions,hcid0);
          double fc_by_pe = FCbyPE(conditions,hcid0);
// for me
//            if ( (fabs(ieta) == 25 && iphi == 25 ) )
          if ( fabs(ieta) <= 30)
          {
           cout<<"close depth4 i/ieta/iphi/phi/ene/gain/mymuons1 "<<ehdepth[i].second-1
           <<"/"<<ieta<<"/"<<iphi<<"/"<<pTrack.phi()*180/TMath::Pi()<<"/"
           <<setprecision(8)<<ehdepth[i].first
           <<"/"<<gain<<"/"<<mymuons1<<endl;
//             cout<<"close px/py "<<pTrack.px()<<"/"<<pTrack.py()<<endl;
          }
          if (gain  != 0) chg  /= gain;
          if ( ehdepth[i].second-1 == 0 ) 
          {
           gain1=gain; charge1=chg; chargec1=charge1*cor; chargeg1=chargec1; chgr1=chargec1/corr; hcal1=ehdepth[i].first; 
           hcalc1=hcal1*cor; hcalg1=hcalc1; resp1=corr; cross1 = 1; 
          }
          if ( ehdepth[i].second-1 == 1 ) 
          {
           gain2=gain; charge2=chg; chargec2=charge2*cor; chargeg2=chargec2; chgr2=chargec2/corr; hcal2=ehdepth[i].first; 
           hcalc2=hcal2*cor; hcalg2=hcalc2; resp2=corr;
          }
          if ( ehdepth[i].second-1 == 2 ) 
          {
           gain3=gain; charge3=chg; chargec3=charge3*cor; chargeg3=chargec3; chgr3=chargec3/corr; hcal3=ehdepth[i].first; 
           hcalc3=hcal3*cor; hcalg3=hcalc3; resp3=corr;
          }
          if ( ehdepth[i].second-1 == 3 ) 
          {
           gain4=gain; charge4=chg; chargec4=charge4*cor; chargeg4=chargec4; chgr4=chargec4/corr; hcal4=ehdepth[i].first; 
           hcalc4=hcal4*cor; hcalg4=hcalc4; resp4=corr;
          }
          if ( ehdepth[i].second-1 == 4 ) 
          {
           gain5=gain; charge5=chg; chargec5=charge5*cor; chargeg5=chargec5; chgr5=chargec5/corr; hcal5=ehdepth[i].first; 
           hcalc5=hcal5*cor; hcalg5=hcalc5; resp5=corr;
          }
          if ( ehdepth[i].second-1 == 5 ) 
          {
           gain6=gain; charge6=chg; chargec6=charge6*cor; chargeg6=chargec6; chgr6=chargec6/corr; hcal6=ehdepth[i].first; 
           hcalc6=hcal6*cor; hcalg6=hcalc6; resp6=corr;
          }
          if ( ehdepth[i].second-1 == 6 ) 
          {
           gain7=gain; charge7=chg; chargec7=charge7*cor; chargeg7=chargec7; chgr7=chargec7/corr; hcal7=ehdepth[i].first; 
           hcalc7=hcal7*cor; hcalg7=hcalc7; resp7=corr;
          }
         }
//m         cout<<"hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl; 
//////////////  propagate back
         GlobalPoint p2;
         is41 = is42 = is43 = is44 = is45 = is46 = is47 = 0;
         is21 = is22 = is23 = is24 = 0;
         i_start2 = 0;
//m         cout<<"depth1 hcal2 istart "<<hcal2<<"/"<<i_start2<<endl;
         if ( hcal2 > 0 )         
         {
          cross2 = 0;
          for ( int i = 0; i < ndepth2; i++ )
          {
           p2 = myfind(&pTrack, depth2_start[ndepth2-i-1]);
//m           cout<<"p2 "<<p2.x()<<"/"<<p2.y()<<"/"<<p2.z()<<endl;
           const DetId closestCell1(mygHB->getClosestCell(p2));
           if ( ieta == HcalDetId(closestCell1).ieta() &&
               iphi == HcalDetId(closestCell1).iphi() )
           {
            cross2++;
            i_start2 = i+1;
            if ( i == 0 ) is21 = 1;
            if ( i == 1 ) is22 = 1;
            if ( i == 2 ) is23 = 1;
            if ( i == 3 ) is24 = 1;
           }
           cout<<"subdet2 back cosmic2 i/ieta/iphi/depth/cross2/istart2 "<<i<<"/"<<ieta<<"/"<<iphi<<"/"
           <<HcalDetId(closestCell1).ieta()<<"/"<<
           HcalDetId(closestCell1).iphi()<<"/"<<HcalDetId(closestCell1).depth()<<"/"<<cross2
           <<"/"<<i_start2<<endl;
          }
         }
//         cross2 /= ndepth2;
         i_start3 = 0;
         if ( hcal3 > 0 )
         {
          cross3 = 0;
//m          cout<<"hcal3 "<<hcal3<<endl;
          for ( int i = 0; i < ndepth3; i++ )
          {
           p2 = myfind(&pTrack, depth3_start[ndepth3-i-1]);
           const DetId closestCell1(mygHB->getClosestCell(p2));
           if ( HcalDetId(closestCell).ieta() == HcalDetId(closestCell1).ieta() &&
               HcalDetId(closestCell).iphi() == HcalDetId(closestCell1).iphi() )
           {
            cross3++;
            i_start3 = i+1;
           }
           cout<<"subdet3 back cosmic2 i/ieta/iphi/depth/cross3 "<<i<<"/"<<ieta<<"/"<<iphi<<"/"
           <<HcalDetId(closestCell1).ieta()<<"/"<<
           HcalDetId(closestCell1).iphi()<<"/"<<HcalDetId(closestCell1).depth()<<"/"<<cross3<<endl;
          }
         }
//         cross3 /= ndepth3;
//         i_start4 = 0;
         if ( hcal4 > 0 )
         {
          cross4 = 0;
          for ( int i = 0; i < ndepth4; i++ )
          {
           p2 = myfind(&pTrack, depth4_start[ndepth4-i-1]);
           const DetId closestCell1(mygHB->getClosestCell(p2));
           if ( HcalDetId(closestCell).ieta() == HcalDetId(closestCell1).ieta() &&
               HcalDetId(closestCell).iphi() == HcalDetId(closestCell1).iphi() )
           {
            cross4++;
            if ( i == 0 ) is41++;
            if ( i == 1 ) is42++;
            if ( i == 2 ) is43++;
            if ( i == 3 ) is44++;
            if ( i == 4 ) is45++;
            if ( i == 5 ) is46++;
            if ( i == 6 ) is47++;
            i_start4 = i+1;
           }
           cout<<"subdet4 back cosmic2 i/ieta/iphi/depth/cross4 "<<i<<"/"<<ieta<<"/"<<iphi<<"/"
           <<HcalDetId(closestCell1).ieta()<<"/"<<
           HcalDetId(closestCell1).iphi()<<"/"<<HcalDetId(closestCell1).depth()<<"/"<<cross4<<endl;
          }
         }
//         cross4 /= ndepth4;
         if ( hcal4 > 0 && cross4 > 6 ) 
         {
          cout<<"****fill4-4 hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl;
          cout<<"****fill4-4 cross1/cross2/cross3/cross4/istart "<<cross1<<"/"<<cross2<<"/"<<cross3<<"/"<<cross4
          <<"/"<<i_start4<<endl;
          hcalg4=hcalc4/cross4; chargeg4=chargec4/cross4; hcr4=hcr4/cross4; chgr4=chgr4/cross4;
          ftree->Fill();
//          goto myend;
          goto step3;
         }
         if ( hcal4 > 0 && cross4 < 7 && cross4 > 0)
         {
          cout<<"****fill4-4 hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl;
          cout<<"****fill4-4 cross1/cross2/cross3/cross4/i_start "<<cross1<<"/"<<cross2<<"/"<<cross3<<"/"<<cross4
          <<"/"<<i_start4<<endl;
          cout<<"is1/is2/is3/is4/is5/is6/is7 "<<is41<<"/"<<is42<<"/"<<is43<<"/"<<is44<<"/"<<is45<<"/"<<is46<<"/"<<is47<<endl;
          hcalg4=hcalc4/cross4; chargeg4=chargec4/cross4; hcr4=hcr4/cross4; chgr4=chgr4/cross4;
          ftree->Fill();
          if ( i_start4 <= (ndepth4-1) && i_start4 > 0 )
          {
//           ftree->Fill();
           goto step4;
          }
          else
           goto step3;
//           goto myend;
         }
         if ( hcal4 < 0 )
         {
          if ( i_start4 < (ndepth4-1) )
          {
           cout<<"depth 4 less "<<i_start4<<endl;
           i_start4++;
           goto step4;
          }
          else
          {
//           goto myend;
           i_start3 = 0;
           goto step3;
          }
         }
//         if ( hcal3 > 0 && cross3 > 0.9 )
//         {
//          cout<<"****fill4-3 hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl;
//          cout<<"****fill4-3 cross1/cross2/cross3/cross4/istart "<<cross1<<"/"<<cross2<<"/"<<cross3<<"/"<<cross4
//          <<"/"<<i_start3<<endl;
//          hcalg3=hcalg3/cross3; chargeg3=chargeg3/cross3; hcr3=hcr3/cross3; chgr3=chgr3/cross3;
//          ftree->Fill();
//          goto step4;
//          goto myend;
//         }
//         if ( hcal3 > 0 && cross3 < 0.9 && cross3 > 0 ) 
//         {
//          cout<<"****fill4-3 hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl;
//          cout<<"****fill4-3 cross1/cross2/cross3/cross4/istart "<<cross1<<"/"<<cross2<<"/"<<cross3<<"/"<<cross4
//          <<"/"<<i_start3<<endl;
//          hcalg3=hcalg3/cross3; chargeg3=chargeg3/cross3; hcr3=hcr3/cross3; chgr3=chgr3/cross3;
//          ftree->Fill();
//          if ( i_start3 <= (ndepth3-1) && i_start3 > 0 )
//          {
//           ftree->Fill();
//           goto step3;
//           if ( i_start4 <= (ndepth4-1) && i_start4 > 0 )
//           {
//            i_start4++;
//            goto step4;
//           }
//           else
//            goto myend;
////           goto myend;
//          }
//          else
//           goto myend;
//           goto step2;
//         }
//         if ( hcal2 > 0 && cross2 > 0.9)
//         {
//          cout<<"****fill1-2 hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl;
//          cout<<"****fill1-2 cross1/cross2/cross3/cross4/istart "<<cross1<<"/"<<cross2<<"/"<<cross3<<"/"<<cross4
//          <<"/"<<i_start2<<endl;
//          hcalg2=hcalg2/cross2; chargeg2=chargeg2/cross2; hcr2=hcr2/cross2; chgr2=chgr2/cross2;
//          ftree->Fill();
//          goto myend;
//          goto step3;
//         }
//         if ( hcal2 > 0 && cross2 < 0.9 && cross2 > 0 ) 
//         {
//          cout<<"****fill1-2 hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl;
//          cout<<"****fill1-2 cross1/cross2/cross3/cross4/istart "<<cross1<<"/"<<cross2<<"/"<<cross3<<"/"<<cross4
//          <<"/"<<i_start2<<endl;
//          hcalg2=hcalg2/cross2; chargeg2=chargeg2/cross2; hcr2=hcr2/cross2; chgr2=chgr2/cross2;
//          ftree->Fill();
//          goto myend;
//          goto step2;
//         }
//         if ( hcal2 < 0 ) 
//         {
//          if ( hcal1 > 0 && cross1 > 0) 
//          {
//           cout<<"****fill1-1 hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl;
//           cout<<"****fill1-1 cross1/cross2/cross3/cross4/istart "<<cross1<<"/"<<cross2<<"/"<<cross3<<"/"<<cross4
//           <<"/"<<i_start2<<endl;
//           hcalg1=hcalg1/cross1; chargeg1=chargeg1/cross1; hcr1=hcr1/cross1; chgr1=chgr1/cross1;
//           ftree->Fill(); 
//           goto myend;
//           goto step2; 
//          }
//          else goto myend;
//         }
//         else goto myend1;
/////////////////////// end propagate back
        } // end fabs(ieta)< 30 && fabs(ieta)>16 && muon_1 
       } //end trackID.okHCAL
      }
///////////////// depth 2
      step2:
      if ( i_start2 > (ndepth2-1) ) goto step1;
      hcal1=hcal2=hcal3=hcal4=hcal5=hcal6=hcal7=-10000;
      hcalc1=hcalc2=hcalc3=hcalc4=hcalc5=hcalc6=hcalc7=-10000;
      gain1=gain2=gain3=gain4=gain5=gain6=gain7=-10000;
      iphi1=iphi2=iphi3=iphi4=iphi5=iphi6=iphi7=-10000;
      ieta1=ieta2=ieta3=ieta4=ieta5=ieta6=ieta7=-10000;
      charge1=charge2=charge3=charge4=charge5=charge6=charge7=-10000;
      chargec1=chargec2=chargec3=chargec4=chargec5=chargec6=chargec7=-10000;
      hcalg1=hcalg2=hcalg3=hcalg4=hcalg5=hcalg6=hcalg7=-10000;
      hcr1=hcr2=hcr3=hcr4=hcr5=hcr6=hcr7=-10000;
      chargeg1=chargeg2=chargeg3=chargeg4=chargeg5=chargeg6=chargeg7=-10000;
      chgr1=chgr2=chgr3=chgr4=chgr5=chgr6=chgr7=-10000;

//      mypoint = myfind(&pTrack, depth2_end[ndepth2-1]);
      cout<<"enter step 2 i_start "<<i_start2<<"/"<<ndepth2-i_start2-1<<endl;
      mypoint = myfind(&pTrack, depth2_start[ndepth2-i_start2-1]);
      if (RecMuon->outerTrack().isNonnull()) 
      {
////////////////// my part 
       const CaloSubdetectorGeometry* mygHB = geo->getSubdetectorGeometry(DetId::Hcal,HcalBarrel);
       const DetId closestCell(mygHB->getClosestCell(mypoint));
       HcalDetId hcidt(closestCell.rawId());  
       HcalDetId check;
//	std::pair<bool,HcalDetId> info = spr::propagateHCALBack(pTrack,  geo, bField, (((verbosity_/100)%10>0)));
       std::pair<bool,HcalDetId> info = propagateHCALBack1(&pTrack,  geo, bField, false);
       if (info.first) { 
	check = info.second;
       }	
       if (closestCell) 
       {
	if ((hcidt.ieta() == check.ieta()) && (hcidt.iphi() == check.iphi()))
        {
	 match= 1;
        }
        else match = 0;
	HcalSubdetector subdet = HcalDetId(closestCell).subdet();
	ieta   = HcalDetId(closestCell).ieta();
	iphi   = HcalDetId(closestCell).iphi();
        HCal_ieta = ieta;
        HCal_iphi = iphi;
        std::pair<double, double> par = mycheckz(ieta,depth2_start[ndepth2-i_start2-1]);
        zmin2 = par.first;
        zmax2 = par.second;
        z2 = mypoint.z();
        cout<<"2 zmin zmax mypoint "<<par.first<<"/"<<par.second<<"/"<<mypoint.z()<<endl;
        par = mycheckx(iphi,depth2_start[ndepth2-i_start2-1]);
        xmin2 = par.first;
        xmax2 = par.second;
        x2 = mypoint.z();
        cout<<"2 xmin xmax "<<par.first<<"/"<<par.second<<"/"<<mypoint.x()<<endl;
        par = mychecky(iphi,depth2_start[ndepth2-i_start2-1]);
        ymin2 = par.first;
        ymax2 = par.second;
        y2 = mypoint.z();
        cout<<"2 ymin ymax "<<par.first<<"/"<<par.second<<"/"<<mypoint.y()<<endl;
        cout<<"subdet depth 2 ieta/iphi/depth "<<ieta<<"/"<<iphi<<"/"<<HcalDetId(closestCell).depth()<<endl;

        cor = mycorrection(&pTrack, ieta, iphi);
	bool            hborhe = (std::abs(ieta) == 16);
        bool muon_1 = RecMuon->p()>1.0;
        bool muon_2 = muon::isMediumMuon(*RecMuon);
        bool muon_3 = (RecMuon->pfIsolationR04().sumChargedHadronPt + std::max(0.,RecMuon->pfIsolationR04().sumNeutralHadronEt + RecMuon->pfIsolationR04().sumPhotonEt - (0.5 *RecMuon->pfIsolationR04().sumPUPt))) / RecMuon->pt()<0.15;
        if ( fabs(ieta)< 40 && fabs(ieta) >= 1 && muon_1 )
        {
	 std::vector<std::pair<double,int> > ehdepth;
	 spr::energyHCALCell((HcalDetId) closestCell, hbhe, ehdepth, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depthHE, false);

         isize1 = ehdepth.size();
         cout<<"depth 2 size "<<isize1<<endl;
         cross1 = cross2 = cross3 = cross4 = 0;
         ifirst = isecond = ithird = iforth = 0;
//         if ( ehdepth.size() == 0 ) goto step3;
         if ( ehdepth.size() != 0 ) isecond = 1;
         myhits = ehdepth.size();
	 for (unsigned int i=0; i<ehdepth.size(); ++i) 
         {
          HcalSubdetector subdet0 = (hborhe) ? ((ehdepth[i].second >= depthHE) ? HcalEndcap : HcalBarrel) : subdet;

          HcalDetId hcid0(subdet0,ieta,iphi,ehdepth[i].second);
          //get charge 
          double ene  = ehdepth[i].first;
          double chg(ene), enec(ene);
          double corr = respCorr(DetId(hcid0));
          if (corr != 0) 
          {
           enec /= corr;
          }
          double gain = gainFactor(conditions,hcid0);
          double pWidth = pedWidth(conditions,hcid0);
          double peds = pedEffective(conditions,hcid0);
          double peds1 = pedestal(conditions,hcid0);
          double dark_current = darkCurrent(conditions,hcid0);
          double fc_by_pe = FCbyPE(conditions,hcid0);
// for me
//            if ( (fabs(ieta) == 25 && iphi == 25 ) )
          if ( fabs(ieta) <= 30)
          {
           cout<<"close depth2 i/ieta/iphi/phi/ene/gain/mymuons1 "<<ehdepth[i].second-1
           <<"/"<<ieta<<"/"<<iphi<<"/"<<pTrack.phi()*180/TMath::Pi()<<"/"
           <<setprecision(8)<<ehdepth[i].first
           <<"/"<<gain<<"/"<<mymuons1<<endl;
//             cout<<"close px/py "<<pTrack.px()<<"/"<<pTrack.py()<<endl;
          }
          if (gain  != 0) chg  /= gain;
          if ( ehdepth[i].second-1 == 0 ) 
          {
           gain1=gain; charge1=chg; chargec1=charge1*cor; chargeg1=chargec1; chgr1=chargec1/corr; hcal1=ehdepth[i].first; 
           hcalc1=hcal1*cor; hcalg1=hcalc1; hcr1=hcalc1/corr; iphi1=iphi; ieta1=ieta; cross1 = 1;
          }
          if ( ehdepth[i].second-1 == 1 ) 
          {
           gain2=gain; charge2=chg; chargec2=charge2*cor; chargeg2=chargec2; chgr2=chargec2/corr; hcal2=ehdepth[i].first; 
           hcalc2=hcal2*cor; hcalg2=hcalc2; hcr2=hcalc2/corr; resp2=corr;
          }
          if ( ehdepth[i].second-1 == 2 ) 
          {
           gain3=gain; charge3=chg; chargec3=charge3*cor; chargeg3=chargec3; chgr3=chargec3/corr; hcal3=ehdepth[i].first; 
           hcalc3=hcal3*cor; hcalg3=hcalc3; hcr3=hcalc3/corr; resp3=corr;
          }
          if ( ehdepth[i].second-1 == 3 ) 
          {
           gain4=gain; charge4=chg; chargec4=charge4*cor; chargeg4=chargec4; chgr4=chargec4/corr; hcal4=ehdepth[i].first; 
           hcalc4=hcal4*cor; hcalg4=hcalc4; hcr4=hcalc4/corr; resp4=corr;
          }
          if ( ehdepth[i].second-1 == 4 ) 
          {
           gain5=gain; charge5=chg; chargec5=charge5*cor; chargeg5=chargec5; chgr5=chargec5/corr; hcal5=ehdepth[i].first; 
           hcalc5=hcal5*cor; hcalg5=hcalc5; hcr5=hcalc5/corr; resp5=corr;
          }
          if ( ehdepth[i].second-1 == 5 ) 
          {
           gain6=gain; charge6=chg; chargec6=charge6*cor; chargeg6=chargec6; chgr6=chargec6/corr; hcal6=ehdepth[i].first; 
           hcalc6=hcal6*cor; hcalg6=hcalc6; hcr6=hcalc6/corr; resp6=corr;
          }
          if ( ehdepth[i].second-1 == 6 ) 
          {
           gain7=gain; charge7=chg; chargec7=charge7*cor; chargeg7=chargec7; chgr7=chargec7/corr; hcal7=ehdepth[i].first; 
           hcalc7=hcal7*cor; hcalg7=hcalc7; hcr7=hcalc7/corr; resp7=corr;
          }
         }
//m         cout<<"hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl; 
//////////////  propagate back
         GlobalPoint p2;
         is21 = is22 = is23 = is24 = 0;
         is41 = is42 = is43 = is44 = is45 = is46 = is47 = 0;
//         i_start2 = 0;
         if ( hcal2 > 0 )
         {
          cross2 = 0;
          for ( int i = 0; i < ndepth2; i++ )
          {
           p2 = myfind(&pTrack, (depth2_start[ndepth2-i-1]));
           const DetId closestCell1(mygHB->getClosestCell(p2));
           if ( HcalDetId(closestCell).ieta() == HcalDetId(closestCell1).ieta() &&
               HcalDetId(closestCell).iphi() == HcalDetId(closestCell1).iphi() )
           {
            cross2++;
            i_start2 = i+1;
            if ( i == 0 ) is21 = 1;
            if ( i == 1 ) is22 = 1;
            if ( i == 2 ) is23 = 1;
            if ( i == 3 ) is24 = 1;
           }
           cout<<"subdet2 step2 back i/ieta/iphi/depth/cross2 "<<i<<"/"<<ieta<<"/"<<iphi<<"/"
           <<HcalDetId(closestCell1).ieta()<<"/"<<
           HcalDetId(closestCell1).iphi()<<"/"<<HcalDetId(closestCell1).depth()<<"/"<<cross2<<endl;
          }
         }
//         cross2 /= ndepth2;
         i_start3 = 0;
         if ( hcal3 > 0 )
         {
          cross3 = 0;
          for ( int i = 0; i < ndepth3; i++ )
          {
           p2 = myfind(&pTrack, depth3_start[ndepth3-i-1]);
           const DetId closestCell1(mygHB->getClosestCell(p2));
           if ( HcalDetId(closestCell).ieta() == HcalDetId(closestCell1).ieta() &&
               HcalDetId(closestCell).iphi() == HcalDetId(closestCell1).iphi() )
           {
            cross3++;
            i_start3 = i+1;
           }
           cout<<"subdet3 step2 back i/ieta/iphi/depth/cross3 "<<i<<"/"<<ieta<<"/"<<iphi<<"/"
           <<HcalDetId(closestCell1).ieta()<<"/"<<
           HcalDetId(closestCell1).iphi()<<"/"<<HcalDetId(closestCell1).depth()<<"/"<<cross3<<endl;
          }
         }
//         cross3 /= ndepth3;
         i_start4 = 0;
         if ( hcal4 > 0 )
         {
          cross4 = 0;
          for ( int i = 0; i < ndepth4; i++ )
          {
           p2 = myfind(&pTrack, depth4_end[ndepth4-i-1]);
           const DetId closestCell1(mygHB->getClosestCell(p2));
           if ( HcalDetId(closestCell).ieta() == HcalDetId(closestCell1).ieta() &&
               HcalDetId(closestCell).iphi() == HcalDetId(closestCell1).iphi() )
           {
            cross4++;
            if ( i == 0 ) is41++;
            if ( i == 1 ) is42++;
            if ( i == 2 ) is43++;
            if ( i == 3 ) is44++;
            if ( i == 4 ) is45++;
            if ( i == 5 ) is46++;
            if ( i == 6 ) is47++;
            i_start4 = i+1;
           }
           cout<<"subdet4 step2 back i/ieta/iphi/depth/cross2 "<<i<<"/"<<ieta<<"/"<<iphi<<"/"
           <<HcalDetId(closestCell1).ieta()<<"/"<<
           HcalDetId(closestCell1).iphi()<<"/"<<HcalDetId(closestCell1).depth()<<"/"<<cross4<<endl;
          }
         }
//         cross4 /= ndepth4;
//         if ( hcal4 > 0 && cross4 > 0.9 ) 
//         {
//          cout<<"****fill2-4 hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl;
//          cout<<"****fill2-4 cross1/cross2/cross3/cross4/istart "<<cross1<<"/"<<cross2<<"/"<<cross3<<"/"<<cross4
//          <<"/"<<i_start4<<endl;
//          hcalg4=hcalg4/cross4; chargeg4=chargeg4/cross4; hcr4=hcr4/cross4; chgr4=chgr4/cross4;
//          ftree->Fill();
//          goto myend;
//         }
//         if ( hcal4 > 0 && cross4 < 0.9 && cross4 > 0 ) 
//         {
//          cout<<"****fill2-4 hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl;
//          cout<<"****fill2-4 cross1/cross2/cross3/cross4/istart "<<cross1<<"/"<<cross2<<"/"<<cross3<<"/"<<cross4
//          <<"/"<<i_start4<<endl;
//          hcalg4=hcalg4/cross4; chargeg4=chargeg4/cross4; hcr4=hcr4/cross4; chgr4=chgr4/cross4;
//          ftree->Fill();
//          goto step4;
//          goto myend;
//         }
//         if ( hcal3 > 0 && cross3 > 0.9 )
//         {
//          hcalg3=hcalg3/cross3; chargeg3=chargeg3/cross3; hcr3=hcr3/cross3; chgr3=chgr3/cross3;
//          cout<<"****fill2-3 hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl;
//          cout<<"****fill2-3 cross1/cross2/cross3/cross4/istart/hcalg3 "<<cross1<<"/"<<cross2<<"/"<<cross3<<"/"<<cross4
//          <<"/"<<i_start3<<"/"<<hcalg3<<endl;
//          ftree->Fill();
//          goto step4;
//         }
//         if ( hcal3 > 0 && cross3 < 0.9 && cross3 > 0 ) 
//         {
//          hcalg3=hcalg3/cross3; chargeg3=chargeg3/cross3; hcr3=hcr3/cross3; chgr3=chgr3/cross3;
//          cout<<"****fill2-3 hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl;
//          cout<<"****fill2-3 cross1/cross2/cross3/cross4/istart/hcalg3 "<<cross1<<"/"<<cross2<<"/"<<cross3<<"/"<<cross4
//          <<"/"<<i_start3<<"/"<<hcalg3<<endl;
//          ftree->Fill();
//          goto step3;
//         }
         if ( hcal2 > 0 && cross2 > 3)
         {
          cout<<"****fill2-2 0.9 hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl;
          cout<<"****fill2-2 cross1/cross2/cross3/cross4/istart "<<cross1<<"/"<<cross2<<"/"<<cross3<<"/"<<cross4
          <<"/"<<i_start2<<endl;
          cout<<"is1/is2/is3/is4 "<<is21<<"/"<<is22<<"/"<<is23<<"/"<<is24<<endl;
          hcalg2=hcalc2/cross2; chargeg2=chargec2/cross2; hcr2=hcr2/cross2; chgr2=chgr2/cross2;
          ftree->Fill();
//          goto step3;
          if ( hcal1 < 0 )
           goto step1;
          else 
           goto myend;
         }
         if ( hcal2 > 0 && cross2 < 4 && cross2 > 0 ) 
         {
          cout<<"fill2-2 <0.9 hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl;
          cout<<"fill2-2 cross1/cross2/cross3/cross4/istart "<<cross1<<"/"<<cross2<<"/"<<cross3<<"/"<<cross4
          <<"/"<<i_start2<<endl;
          cout<<"is1/is2/is3/is4 "<<is21<<"/"<<is22<<"/"<<is23<<"/"<<is24<<endl;
          hcalg2=hcalc2/cross2; chargeg2=chargec2/cross2; hcr2=hcr2/cross2; chgr2=chgr2/cross2;
          ftree->Fill();
          if ( i_start2 <= (ndepth2-1) && i_start2 > 0 ) 
          {
           cout<<"goto step2 "<<i_start2<<endl;
//           i_start2++;
           goto step2;
          }
          else
          {
           cout<<"goto end step2 "<<ndepth2-1<<endl;
           if ( hcal1 < 0 )
            goto step1;
           else
            goto myend;
          }
         }
         if ( hcal2 < 0 )
         {
          if ( i_start2 < (ndepth2-1) )
          {
           cout<<"depth 2 less "<<i_start2<<endl;
           i_start2++;
           goto step2;
          }
          else
          {
           if ( hcal1 < 0 )
            goto step1;
           else
           {
            cout<<"****fill2-1 hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl;
            cout<<"****fill2-1 cross1/cross2/cross3/cross4/istart "<<cross1<<"/"<<cross2<<"/"<<cross3<<"/"<<cross4
            <<"/"<<i_start2<<endl;
            cross1=1; hcalg1=hcalg1/cross1; chargeg1=chargeg1/cross1; hcr1=hcr1/cross1; chgr1=chgr1/cross1;
            ftree->Fill(); 
            goto myend;
           }
          }
//           goto step3;
         }
//         if ( hcal2 < 0 ) 
//         {
//          if ( hcal1 > 0 && cross1 > 0 ) 
//          {
//           cout<<"****fill2-1 hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl;
//           cout<<"****fill2-1 cross1/cross2/cross3/cross4/istart "<<cross1<<"/"<<cross2<<"/"<<cross3<<"/"<<cross4
//           <<"/"<<i_start2<<endl;
//           hcalg1=hcalg1/cross1; chargeg1=chargeg1/cross1; hcr1=hcr1/cross1; chgr1=chgr1/cross1;
//           ftree->Fill(); 
//           goto step3; 
//          }
//          goto step3;
//         }
//         else goto myend1;
/////////////////////// end propagate back
        } // end fabs(ieta)< 30 && fabs(ieta)>16 && muon_1
       } //end trackID.okHCAL
      }
///////////////// end depth 2	
///////////////// depth 3
//      step3_1:
//      i_start3 = 0;
      step3:
      if ( i_start3 > (ndepth3-1) ) goto step2;
      hcal1=hcal2=hcal3=hcal4=hcal5=hcal6=hcal7=-10000;
      hcalc1=hcalc2=hcalc3=hcalc4=hcalc5=hcalc6=hcalc7=-10000;
      gain1=gain2=gain3=gain4=gain5=gain6=gain7=-10000;
      iphi1=iphi2=iphi3=iphi4=iphi5=iphi6=iphi7=-10000;
      ieta1=ieta2=ieta3=ieta4=ieta5=ieta6=ieta7=-10000;
      charge1=charge2=charge3=charge4=charge5=charge6=charge7=-10000;
      chargec1=chargec2=chargec3=chargec4=chargec5=chargec6=chargec7=-10000;
      hcalg1=hcalg2=hcalg3=hcalg4=hcalg5=hcalg6=hcalg7=-10000;
      hcr1=hcr2=hcr3=hcr4=hcr5=hcr6=hcr7=-10000;
      chargeg1=chargeg2=chargeg3=chargeg4=chargeg5=chargeg6=chargeg7=-10000;
      chgr1=chgr2=chgr3=chgr4=chgr5=chgr6=chgr7=-10000;

//      mypoint = myfind(&pTrack, depth3_end[ndepth3-1]);
      cout<<"enter step 3 i_start "<<i_start3<<"/"<<ndepth3-i_start3-1<<endl;
      mypoint = myfind(&pTrack, depth3_start[ndepth3-i_start3-1]);
      cout<<"mypoint step 3 "<<mypoint.x()<<"/"<<mypoint.y()<<"/"<<mypoint.z()<<endl;
      if (RecMuon->outerTrack().isNonnull()) 
      {
       const CaloSubdetectorGeometry* mygHB = geo->getSubdetectorGeometry(DetId::Hcal,HcalBarrel);
       const DetId closestCell(mygHB->getClosestCell(mypoint));
       HcalDetId hcidt(closestCell.rawId());  
//       if ( fabs(hcidt.ieta()) <= 30)
//        cout<<"init ieta/iphi "<<hcidt.ieta()<<"/"<<hcidt.iphi()<<endl;

       HcalDetId check;
//	std::pair<bool,HcalDetId> info = spr::propagateHCALBack(pTrack,  geo, bField, (((verbosity_/100)%10>0)));
       std::pair<bool,HcalDetId> info = propagateHCALBack1(&pTrack,  geo, bField, false);
       if (info.first) { 
	check = info.second;
       }	
       if (closestCell) 
       {
	if ((hcidt.ieta() == check.ieta()) && (hcidt.iphi() == check.iphi()))
        {
	 match= 1;
        }
        else match = 0;
	HcalSubdetector subdet = HcalDetId(closestCell).subdet();
	ieta   = HcalDetId(closestCell).ieta();
	iphi   = HcalDetId(closestCell).iphi();
        HCal_ieta = ieta;
        HCal_iphi = iphi;
        std::pair<double, double> par = mycheckz(ieta,depth3_start[ndepth3-i_start3-1]);
        zmin3 = par.first;
        zmax3 = par.second;
        z3 = mypoint.z();
        cout<<"depth 3 zmin zmax mypoint "<<par.first<<"/"<<par.second<<"/"<<mypoint.z()<<endl;
        par = mycheckx(iphi,depth3_start[ndepth3-i_start3-1]);
        xmin3 = par.first;
        xmax3 = par.second;
        x3 = mypoint.z();
        cout<<"depth 3 xmin xmax "<<par.first<<"/"<<par.second<<"/"<<mypoint.x()<<endl;
        par = mychecky(iphi,depth3_start[ndepth3-i_start3-1]);
        ymin3 = par.first;
        ymax3 = par.second;
        y3 = mypoint.z();
        cout<<"depth 3 ymin ymax "<<par.first<<"/"<<par.second<<"/"<<mypoint.y()<<endl;
        cout<<"subdet 3 ieta/iphi/depth "<<ieta<<"/"<<iphi<<"/"<<HcalDetId(closestCell).depth()<<endl;

        cor = mycorrection(&pTrack, ieta, iphi);
	bool            hborhe = (std::abs(ieta) == 16);
        bool muon_1 = RecMuon->p()>1.0;
//          bool muon_1 = RecMuon->pt()>20.0;
        bool muon_2 = muon::isMediumMuon(*RecMuon);
        bool muon_3 = (RecMuon->pfIsolationR04().sumChargedHadronPt + std::max(0.,RecMuon->pfIsolationR04().sumNeutralHadronEt + RecMuon->pfIsolationR04().sumPhotonEt - (0.5 *RecMuon->pfIsolationR04().sumPUPt))) / RecMuon->pt()<0.15;
        if ( fabs(ieta)< 30 && fabs(ieta) >= 1 && muon_1 )
        {
//           bool okE = trackID.okECAL;
	 std::vector<std::pair<double,int> > ehdepth;
	 spr::energyHCALCell((HcalDetId) closestCell, hbhe, ehdepth, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depthHE, false);

         isize1 = ehdepth.size();
         cout<<"depth 3 size "<<isize1<<endl;
         cross1 = cross2 = cross3 = cross4 = 0;
//         if ( ehdepth.size() == 0 ) goto step4;
         ifirst = isecond = ithird = iforth = 0;
         if ( ehdepth.size() != 0 ) ithird = 1;
         myhits = ehdepth.size();
	 for (unsigned int i=0; i<ehdepth.size(); ++i) 
         {
          HcalSubdetector subdet0 = (hborhe) ? ((ehdepth[i].second >= depthHE) ? HcalEndcap : HcalBarrel) : subdet;

          HcalDetId hcid0(subdet0,ieta,iphi,ehdepth[i].second);
          //get charge 
          double ene  = ehdepth[i].first;
          double chg(ene), enec(ene);
          double corr = respCorr(DetId(hcid0));
          if (corr != 0) 
          {
           enec /= corr;
          }
          double gain = gainFactor(conditions,hcid0);
          double pWidth = pedWidth(conditions,hcid0);
          double peds = pedEffective(conditions,hcid0);
          double peds1 = pedestal(conditions,hcid0);
          double dark_current = darkCurrent(conditions,hcid0);
          double fc_by_pe = FCbyPE(conditions,hcid0);
// for me
//            if ( (fabs(ieta) == 25 && iphi == 25 ) )
          if ( fabs(ieta) <= 30)
          {
           cout<<"close depth3 i/ieta/iphi/phi/ene/gain/mymuons1 "<<ehdepth[i].second-1
           <<"/"<<ieta<<"/"<<iphi<<"/"<<pTrack.phi()*180/TMath::Pi()<<"/"
           <<setprecision(8)<<ehdepth[i].first
           <<"/"<<gain<<"/"<<mymuons1<<endl;
//             cout<<"close px/py "<<pTrack.px()<<"/"<<pTrack.py()<<endl;
          }
          if (gain  != 0) chg  /= gain;
          if ( ehdepth[i].second-1 == 0 ) 
          {
           gain1=gain; charge1=chg; chargec1=charge1*cor; chargeg1=chargec1; chgr1=chargec1/corr; hcal1=ehdepth[i].first; 
           hcalc1=hcal1*cor; hcalg1=hcalc1; hcr1=hcalc1/corr; iphi1=iphi; ieta1=ieta;
          }
          if ( ehdepth[i].second-1 == 1 ) 
          {
           gain2=gain; charge2=chg; chargec2=charge2*cor; chargeg2=chargec2; chgr2=chargec2/corr; hcal2=ehdepth[i].first; 
           hcalc2=hcal2*cor; hcalg2=hcalc2; hcr2=hcalc2/corr; iphi2=iphi; ieta2=ieta;
          }
          if ( ehdepth[i].second-1 == 2 ) 
          {
           gain3=gain; charge3=chg; chargec3=charge3*cor; chargeg3=chargec3; chgr3=chargec3/corr; hcal3=ehdepth[i].first; 
           hcalc3=hcal3*cor; hcalg3=hcalc3; hcr3=hcalc3/corr; resp3=corr;
          }
          if ( ehdepth[i].second-1 == 3 ) 
          {
           gain4=gain; charge4=chg; chargec4=charge4*cor; chargeg4=chargec4; chgr4=chargec4/corr; hcal4=ehdepth[i].first; 
           hcalc4=hcal4*cor; hcalg4=hcalc4; hcr4=hcalc4/corr; resp4=corr;
          }
          if ( ehdepth[i].second-1 == 4 ) 
          {
           gain5=gain; charge5=chg; chargec5=charge5*cor; chargeg5=chargec5; chgr5=chargec5/corr; hcal5=ehdepth[i].first; 
           hcalc5=hcal5*cor; hcalg5=hcalc5; hcr5=hcalc5/corr; resp5=corr;
          }
          if ( ehdepth[i].second-1 == 5 ) 
          {
           gain6=gain; charge6=chg; chargec6=charge6*cor; chargeg6=chargec6; chgr6=chargec6/corr; hcal6=ehdepth[i].first; 
           hcalc6=hcal6*cor; hcalg6=hcalc6; hcr6=hcalc6/corr; resp6=corr;
          }
          if ( ehdepth[i].second-1 == 6 ) 
          {
           gain7=gain; charge7=chg; chargec7=charge7*cor; chargeg7=chargec7; chgr7=chargec7/corr; hcal7=ehdepth[i].first; 
           hcalc7=hcal7*cor; hcalg7=hcalc7; hcr7=hcalc7/corr; resp7=corr;
          }
         }
//m         cout<<"hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl; 
//////////////  propagate back
         GlobalPoint p1, p2, p3;
         is21 = is22 = is23 = is24 = 0;
//         i_start3 = 0;
//         i_start4 = 0;
         if ( hcal4 > 0 )
         {
          cross4 = 0;
          for ( int i = 0; i < ndepth4; i++ )
          {
           p2 = myfind(&pTrack, depth4_end[ndepth4-i-1]);
           const DetId closestCell1(mygHB->getClosestCell(p2));
           if ( HcalDetId(closestCell).ieta() == HcalDetId(closestCell1).ieta() &&
               HcalDetId(closestCell).iphi() == HcalDetId(closestCell1).iphi() )
           {
            cross4++;
            if ( i == 0 ) is41++;
            if ( i == 1 ) is42++;
            if ( i == 2 ) is43++;
            if ( i == 3 ) is44++;
            if ( i == 4 ) is45++;
            if ( i == 5 ) is46++;
            if ( i == 6 ) is47++;
            i_start4 = i+1;
           }
           cout<<"subdet4 step2 back i/ieta/iphi/depth/cross2 "<<i<<"/"<<ieta<<"/"<<iphi<<"/"
           <<HcalDetId(closestCell1).ieta()<<"/"<<
           HcalDetId(closestCell1).iphi()<<"/"<<HcalDetId(closestCell1).depth()<<"/"<<cross4<<endl;
          }
         }
         if ( hcal3 > 0 )
         {
          cross3 = 0;
          for ( int i = 0; i < ndepth3; i++ )
          {
           p2 = myfind(&pTrack, depth3_start[ndepth3-i-1]);
           const DetId closestCell1(mygHB->getClosestCell(p2));
           if ( HcalDetId(closestCell).ieta() == HcalDetId(closestCell1).ieta() &&
               HcalDetId(closestCell).iphi() == HcalDetId(closestCell1).iphi() )
           {
            cross3++;
            i_start3 = i+1;
           }
           cout<<"subdet3 back i/ieta/iphi/depth/cross3 "<<i<<"/"<<ieta<<"/"<<iphi<<"/"
           <<HcalDetId(closestCell1).ieta()<<"/"<<
           HcalDetId(closestCell1).iphi()<<"/"<<HcalDetId(closestCell1).depth()<<"/"<<cross3<<endl;
          }
         }
//         cross3 /= ndepth3;
         i_start2 = 0;
         if ( hcal2 > 0 )
         {
          cross2 = 0;
          for ( int i = 0; i < ndepth2; i++ )
          {
           p2 = myfind(&pTrack, depth2_end[ndepth2-i-1]);
           const DetId closestCell1(mygHB->getClosestCell(p2));
           if ( HcalDetId(closestCell).ieta() == HcalDetId(closestCell1).ieta() &&
               HcalDetId(closestCell).iphi() == HcalDetId(closestCell1).iphi() )
           {
            cross2++;
            i_start2 = i+1;
            if ( i == 0 ) is21 = 1;
            if ( i == 1 ) is22 = 1;
            if ( i == 2 ) is23 = 1;
            if ( i == 3 ) is24 = 1;
           }
           cout<<"subdet2 back i/ieta/iphi/depth/cross2 "<<i<<"/"<<ieta<<"/"<<iphi<<"/"
           <<HcalDetId(closestCell1).ieta()<<"/"<<
           HcalDetId(closestCell1).iphi()<<"/"<<HcalDetId(closestCell1).depth()<<"/"<<cross2<<endl;
          }
         }
//         cross2 /= ndepth2;
         if ( hcal3 > 0 && cross3 > 4 ) 
         {
          cout<<"fill3-3 hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl;
          cout<<"fill3-3 cross1/cross2/cross3/cross4/istart "<<cross1<<"/"<<cross2<<"/"<<cross3<<"/"<<cross4
          <<"/"<<i_start3<<endl;
          hcalg3=hcalc3/cross3; chargeg3=chargec3/cross3; hcr3=hcr3/cross3; chgr3=chgr3/cross3;
          ftree->Fill();
//          goto myend;
          goto step2;
         }
         if ( hcal3 > 0 && cross3 < 5 && cross3 > 0 ) 
         {
          cout<<"fill3-3 hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl;
          cout<<"fill3-3 cross1/cross2/cross3/cross4/istart "<<cross1<<"/"<<cross2<<"/"<<cross3<<"/"<<cross4
          <<"/"<<i_start3<<endl;
          hcalg3=hcalc3/cross4; chargeg3=chargec3/cross3; hcr3=hcr3/cross3; chgr3=chgr3/cross3;
          ftree->Fill();
          if ( i_start3 <= (ndepth3-1) && i_start3 > 0 )
          {
//           ftree->Fill();
           goto step3;
          }
          else
//           goto myend;
           goto step2;
         }
         if ( hcal3 < 0 )
         {
          if ( i_start3 < (ndepth3-1) )
          {
           cout<<"depth 3 less "<<i_start3<<endl;
           i_start3++;
           goto step3;
          }
          else
//           goto myend;
           goto step2;
         }
//         if ( hcal3 > 0 && cross3 > 0 )
//         {
//          cout<<"fill3-3 hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl;
//          cout<<"fill3-3 cross1/cross2/cross3/cross4/istart "<<cross1<<"/"<<cross2<<"/"<<cross3<<"/"<<cross4
//          <<"/"<<i_start3<<endl;
//          hcalg3=hcalg3/cross3; chargeg3=chargec3/cross3; hcr3=hcr3/cross3; chgr3=chgr3/cross3;
//          ftree->Fill();
//          goto step4;
//         }

//         else goto myend1;
/////////////////////// end propagate back
        } // end fabs(ieta)< 30 && fabs(ieta)>16 && muon_1
       } //end trackID.okHCAL
      }
///////////////// end depth 3	
///////////////// depth 1
      step1:
      hcal1=hcal2=hcal3=hcal4=hcal5=hcal6=hcal7=-10000;
      hcalc1=hcalc2=hcalc3=hcalc4=hcalc5=hcalc6=hcalc7=-10000;
      gain1=gain2=gain3=gain4=gain5=gain6=gain7=-10000;
      iphi1=iphi2=iphi3=iphi4=iphi5=iphi6=iphi7=-10000;
      ieta1=ieta2=ieta3=ieta4=ieta5=ieta6=ieta7=-10000;
      charge1=charge2=charge3=charge4=charge5=charge6=charge7=-10000;
      chargec1=chargec2=chargec3=chargec4=chargec5=chargec6=chargec7=-10000;
      hcalg1=hcalg2=hcalg3=hcalg4=hcalg5=hcalg6=hcalg7=-10000;
      hcr1=hcr2=hcr3=hcr4=hcr5=hcr6=hcr7=-10000;
      chargeg1=chargeg2=chargeg3=chargeg4=chargeg5=chargeg6=chargeg7=-10000;
      chgr1=chgr2=chgr3=chgr4=chgr5=chgr6=chgr7=-10000;

      cout<<"enter step 1 i_start "<<i_start4<<endl;
      mypoint = myfind(&pTrack, depth1_start);
//      cout<<"mypoint 2 "<<mypoint.x()<<"/"<<mypoint.y()<<"/"<<mypoint.z()<<endl;
      if (RecMuon->outerTrack().isNonnull()) 
      {
       const CaloSubdetectorGeometry* mygHB = geo->getSubdetectorGeometry(DetId::Hcal,HcalBarrel);
       const DetId closestCell(mygHB->getClosestCell(mypoint));
       HcalDetId hcidt(closestCell.rawId());  

       HcalDetId check;
//	std::pair<bool,HcalDetId> info = spr::propagateHCALBack(pTrack,  geo, bField, (((verbosity_/100)%10>0)));
       std::pair<bool,HcalDetId> info = propagateHCALBack1(&pTrack,  geo, bField, false);
       if (info.first) { 
	check = info.second;
       }	
       
       if (closestCell) 
       {
	if ((hcidt.ieta() == check.ieta()) && (hcidt.iphi() == check.iphi()))
        {
//           if ( (hcidt.ieta()) <= 10 && hcidt.ieta() > 0 && hcidt.iphi() <= 2 )
         cout<<"check1 ieta/ietac/iphi/iphic "<<hcidt.ieta()<<"/"<<check.ieta()<<"/"
         <<hcidt.iphi()<<"/"<<check.iphi()<<endl;
	 match= 1;
        }
        else match = 0;
        std::pair<double, double> par = mycheckz(ieta,depth1_start);
        zmin1 = par.first;
        zmax1 = par.second;
        z1 = mypoint.z();
        cout<<"1 zmin zmax mypoint "<<par.first<<"/"<<par.second<<"/"<<mypoint.z()<<endl;
        par = mycheckx(iphi,depth1_start);
        xmin1 = par.first;
        xmax1 = par.second;
        x1 = mypoint.z();
        cout<<"1 xmin xmax "<<par.first<<"/"<<par.second<<"/"<<mypoint.x()<<endl;
        par = mychecky(iphi,depth1_start);
        ymin1 = par.first;
        ymax1 = par.second;
        y1 = mypoint.z();
        cout<<"1 ymin ymax "<<par.first<<"/"<<par.second<<"/"<<mypoint.y()<<endl;
	HcalSubdetector subdet = HcalDetId(closestCell).subdet();
	ieta   = HcalDetId(closestCell).ieta();
	iphi   = HcalDetId(closestCell).iphi();
        HCal_ieta = ieta;
        HCal_iphi = iphi;
        cout<<"subdet 1 ieta/iphi/depth "<<ieta<<"/"<<iphi<<"/"<<HcalDetId(closestCell).depth()<<endl;

        cor = mycorrection(&pTrack, ieta, iphi);
	bool            hborhe = (std::abs(ieta) == 16);
        bool muon_1 = RecMuon->p()>1.0;
//          bool muon_1 = RecMuon->pt()>20.0;
        bool muon_2 = muon::isMediumMuon(*RecMuon);
        bool muon_3 = (RecMuon->pfIsolationR04().sumChargedHadronPt + std::max(0.,RecMuon->pfIsolationR04().sumNeutralHadronEt + RecMuon->pfIsolationR04().sumPhotonEt - (0.5 *RecMuon->pfIsolationR04().sumPUPt))) / RecMuon->pt()<0.15;
        if ( fabs(ieta)< 40 && fabs(ieta) >= 1 && muon_1 )
        {
//           bool okE = trackID.okECAL;
	 std::vector<std::pair<double,int> > ehdepth;
	 spr::energyHCALCell((HcalDetId) closestCell, hbhe, ehdepth, maxDepth_, -100.0, -100.0, -100.0, -100.0, -500.0, 500.0, useRaw_, depthHE, false);

         isize1 = ehdepth.size();
         cout<<"depth 1 size "<<isize1<<endl;
         cross2 = cross3 = cross4 = 0;
//         if ( ehdepth.size() == 0 ) goto step2;
         ifirst = isecond = ithird = iforth = 0;
         if ( ehdepth.size() != 0 ) iforth = 1;
         myhits = ehdepth.size();
	 for (unsigned int i=0; i<ehdepth.size(); ++i) 
         {
          HcalSubdetector subdet0 = (hborhe) ? ((ehdepth[i].second >= depthHE) ? HcalEndcap : HcalBarrel) : subdet;

          HcalDetId hcid0(subdet0,ieta,iphi,ehdepth[i].second);
          //get charge 
          double ene  = ehdepth[i].first;
          double chg(ene), enec(ene);
          double corr = respCorr(DetId(hcid0));
          if (corr != 0) 
          {
           enec /= corr;
          }
          double gain = gainFactor(conditions,hcid0);
          double pWidth = pedWidth(conditions,hcid0);
          double peds = pedEffective(conditions,hcid0);
          double peds1 = pedestal(conditions,hcid0);
          double dark_current = darkCurrent(conditions,hcid0);
          double fc_by_pe = FCbyPE(conditions,hcid0);
// for me
//            if ( (fabs(ieta) == 25 && iphi == 25 ) )
          if ( fabs(ieta) <= 30)
          {
           cout<<"close depth1 i/ieta/iphi/phi/ene/gain/mymuons1 "<<ehdepth[i].second-1
           <<"/"<<ieta<<"/"<<iphi<<"/"<<pTrack.phi()*180/TMath::Pi()<<"/"
           <<setprecision(8)<<ehdepth[i].first
           <<"/"<<gain<<"/"<<mymuons1<<endl;
//             cout<<"close px/py "<<pTrack.px()<<"/"<<pTrack.py()<<endl;
          }
          if (gain  != 0) chg  /= gain;
          if ( ehdepth[i].second-1 == 0 ) 
          {
           gain1=gain; charge1=chg; chargec1=charge1*cor; chargeg1=chargec1; chgr1=chargec1/corr; hcal1=ehdepth[i].first; 
           hcalc1=hcal1*cor; hcalg1=hcalc1; hcr1=hcalc1/corr; iphi1=iphi; ieta1=ieta; cross1 = 1;
          }
          if ( ehdepth[i].second-1 == 1 ) 
          {
           gain2=gain; charge2=chg; chargec2=charge2*cor; chargeg2=chargec2; chgr2=chargec2/corr; hcal2=ehdepth[i].first; 
           hcalc2=hcal2*cor; hcalg2=hcalc2; hcr2=hcalc2/corr; iphi2=iphi; ieta2=ieta;
          }
          if ( ehdepth[i].second-1 == 2 ) 
          {
           gain3=gain; charge3=chg; chargec3=charge3*cor; chargeg3=chargec3; chgr3=chargec3/corr; hcal3=ehdepth[i].first; 
           hcalc3=hcal3*cor; hcalg3=hcalc3; hcr3=hcalc3/corr; iphi3=iphi; ieta3=ieta;
          }
          if ( ehdepth[i].second-1 == 3 ) 
          {
           gain4=gain; charge4=chg; chargec4=charge4*cor; chargeg4=chargec4; chgr4=chargec4/corr; hcal4=ehdepth[i].first; 
           hcalc4=hcal4*cor; hcalg4=hcalc4; hcr4=hcalc4/corr; resp4=corr;
          }
          if ( ehdepth[i].second-1 == 4 ) 
          {
           gain5=gain; charge5=chg; chargec5=charge5*cor; chargeg5=chargec5; chgr5=chargec5/corr; hcal5=ehdepth[i].first; 
           hcalc5=hcal5*cor; hcalg5=hcalc5; hcr5=hcalc5/corr; resp5=corr;
          }
          if ( ehdepth[i].second-1 == 5 ) 
          {
           gain6=gain; charge6=chg; chargec6=charge6*cor; chargeg6=chargec6; chgr6=chargec6/corr; hcal6=ehdepth[i].first; 
           hcalc6=hcal6*cor; hcalg6=hcalc6; hcr6=hcalc6/corr; resp6=corr;
          }
          if ( ehdepth[i].second-1 == 6 ) 
          {
           gain7=gain; charge7=chg; chargec7=charge7*cor; chargeg7=chargec7; chgr7=chargec7/corr; hcal7=ehdepth[i].first; 
           hcalc7=hcal7*cor; hcalg7=hcalc7; hcr7=hcalc7/corr; resp7=corr;
          }

         }
//m         cout<<"hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl; 
//////////////  propagate back
         GlobalPoint p1, p2;
//         if ( hcal4 > 0 )
//         {
//          cross4 = 0;
//          for ( int i = 0; i < ndepth4; i++ )
//          for ( int i = i_start4; i < ndepth4; i++ )
//          {
//           p2 = myfind(&pTrack, depth4_end[i]);
//           const DetId closestCell1(mygHB->getClosestCell(p2));
//           if ( HcalDetId(closestCell).ieta() == HcalDetId(closestCell1).ieta() &&
//               HcalDetId(closestCell).iphi() == HcalDetId(closestCell1).iphi() )
//            cross4++;
//           else
//           {
//            if ( i_start4 == 0 ) i_start4 = i;
//           }
//m           cout<<"subdet4 back i/ieta/iphi/depth/corss4 "<<i<<"/"<<ieta<<"/"<<iphi<<"/"
//           <<HcalDetId(closestCell1).ieta()<<"/"<<
//           HcalDetId(closestCell1).iphi()<<"/"<<HcalDetId(closestCell1).depth()<<"/"<<cross4<<endl;
//          }
//         }
//         cross4 /= ndepth4;
         if ( hcal1 > 0)
         {
          cout<<"****fill1-1 hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl;
          cout<<"****fill1-1 cross1/cross2/cross3/cross4/istart1 "<<cross1<<"/"<<cross2<<"/"<<cross3<<"/"<<cross4
          <<"/"<<i_start1<<endl;
//          hcalg4=hcalg4/cross4; chargeg4=chargec4/cross4; hcr4=hcr4/cross4; chgr4=chgr4/cross4;
          ftree->Fill();
          goto myend;
         }
         else
          goto myend;
/////////////////////// end propagate back
        } // end fabs(ieta)< 30 && fabs(ieta)>16 && muon_1
       } //end trackID.okHCAL
      }
///////////////// end depth 3	
      myend:
      cout<<""<<endl;
//      cout<<"end this track hcal4/cross4 "<<hcal4<<"/"<<cross4<<endl;
//      cout<<"hcal1/hcal2/hcal3/hcal4 "<<hcal1<<"/"<<hcal2<<"/"<<hcal3<<"/"<<hcal4<<endl;
//      cout<<"cross1/cross2/cross3/cross4 "<<cross1<<"/"<<cross2<<"/"<<cross3<<"/"<<cross4<<endl;
//      myend1:
//m      cout<<"not end "<<endl;
     }
    } // muon collection
  } // muon valid
  

//  if (accept) ftree->Fill();
//  cout<<"------------------end accept "<<endl;
}

// ------------ method called once each job just before starting event loop  ------------
void HBHEMuonCosmic2::beginJob() {

//  out_file = new TFile(outname.c_str(), "RECREATE");
  out_file = new TFile("file.root", "RECREATE");

//  ftree = fs->make<TTree>("TREE", "TREE");
  ftree = new TTree("ftree", "JetHT dataset");
//  ftree1 = new TTree("ftree1", "depths");

  ftree->Branch("event",             &eventNumber_);
  ftree->Branch("run",               &runNumber_);
  ftree->Branch("eta",               &eta);
  ftree->Branch("theta",             &theta);
  ftree->Branch("phi",               &phi);
  ftree->Branch("p",                 &p);
  ftree->Branch("pt",                &pt);
  ftree->Branch("px",                &px);
  ftree->Branch("py",                &py);
  ftree->Branch("pz",                &pz);
  ftree->Branch("vx",                &vx);
  ftree->Branch("vy",                &vy);
  ftree->Branch("vz",                &vz);
  ftree->Branch("match",             &match);
  ftree->Branch("cor",               &cor);
  ftree->Branch("HCal_ieta",         &HCal_ieta);
  ftree->Branch("HCal_iphi",         &HCal_iphi);
  ftree->Branch("isoR04",            &isoR04);
  ftree->Branch("stations",          &stations);
  ftree->Branch("ifirst",            &ifirst);
  ftree->Branch("isecond",           &isecond);
  ftree->Branch("ithird",            &ithird);
  ftree->Branch("iforth",            &iforth);
  ftree->Branch("myhits",            &myhits);

  ftree->Branch("is41",              &is41);
  ftree->Branch("is42",              &is42);
  ftree->Branch("is43",              &is43);
  ftree->Branch("is44",              &is44);
  ftree->Branch("is45",              &is45);
  ftree->Branch("is46",              &is46);
  ftree->Branch("is47",              &is47);

  ftree->Branch("is21",              &is21);
  ftree->Branch("is22",              &is22);
  ftree->Branch("is23",              &is23);
  ftree->Branch("is24",              &is24);

  ftree->Branch("nmuons",            &nmuons);

  ftree->Branch("xmin1",             &xmin1);
  ftree->Branch("xmax1",             &xmax1);
  ftree->Branch("x1",                &x1);
  ftree->Branch("ymin1",             &ymin1);
  ftree->Branch("ymax1",             &ymax1);
  ftree->Branch("y1",                &y1);
  ftree->Branch("zmin1",             &zmin1);
  ftree->Branch("zmax1",             &zmax1);
  ftree->Branch("z1",                &z1);

  ftree->Branch("xmin2",             &xmin2);
  ftree->Branch("xmax2",             &xmax2);
  ftree->Branch("x2",                &x2);
  ftree->Branch("ymin2",             &ymin2);
  ftree->Branch("ymax2",             &ymax2);
  ftree->Branch("y2",                &y2);
  ftree->Branch("zmin2",             &zmin2);
  ftree->Branch("zmax2",             &zmax2);
  ftree->Branch("z2",                &z2);

  ftree->Branch("xmin3",             &xmin3);
  ftree->Branch("xmax3",             &xmax3);
  ftree->Branch("x3",                &x3);
  ftree->Branch("ymin3",             &ymin3);
  ftree->Branch("ymax3",             &ymax3);
  ftree->Branch("y3",                &y3);
  ftree->Branch("zmin3",             &zmin3);
  ftree->Branch("zmax3",             &zmax3);
  ftree->Branch("z3",                &z3);

  ftree->Branch("xmin4",             &xmin4);
  ftree->Branch("xmax4",             &xmax4);
  ftree->Branch("x4",                &x4);
  ftree->Branch("ymin4",             &ymin4);
  ftree->Branch("ymax4",             &ymax4);
  ftree->Branch("y4",                &y4);
  ftree->Branch("zmin4",             &zmin4);
  ftree->Branch("zmax4",             &zmax4);
  ftree->Branch("z4",                &z4);

  ftree->Branch("rmin",              &rmin);
  ftree->Branch("r0",                &r0);


  ftree->Branch("hcal1",             &hcal1);
  ftree->Branch("hcal2",             &hcal2);
  ftree->Branch("hcal3",             &hcal3);
  ftree->Branch("hcal4",             &hcal4);
  ftree->Branch("hcal5",             &hcal5);
  ftree->Branch("hcal6",             &hcal6);
  ftree->Branch("hcal7",             &hcal7);

  ftree->Branch("hcalc1",            &hcalc1);
  ftree->Branch("hcalc2",            &hcalc2);
  ftree->Branch("hcalc3",            &hcalc3);
  ftree->Branch("hcalc4",            &hcalc4);
  ftree->Branch("hcalc5",            &hcalc5);
  ftree->Branch("hcalc6",            &hcalc6);
  ftree->Branch("hcalc7",            &hcalc7);

  ftree->Branch("hcalg1",            &hcalg1);
  ftree->Branch("hcalg2",            &hcalg2);
  ftree->Branch("hcalg3",            &hcalg3);
  ftree->Branch("hcalg4",            &hcalg4);
  ftree->Branch("hcalg5",            &hcalg5);
  ftree->Branch("hcalg6",            &hcalg6);
  ftree->Branch("hcalg7",            &hcalg7);

  ftree->Branch("charge1",           &charge1);
  ftree->Branch("charge2",           &charge2);
  ftree->Branch("charge3",           &charge3);
  ftree->Branch("charge4",           &charge4);
  ftree->Branch("charge5",           &charge5);
  ftree->Branch("charge6",           &charge6);
  ftree->Branch("charge7",           &charge7);

  ftree->Branch("chargec1",          &chargec1);
  ftree->Branch("chargec2",          &chargec2);
  ftree->Branch("chargec3",          &chargec3);
  ftree->Branch("chargec4",          &chargec4);
  ftree->Branch("chargec5",          &chargec5);
  ftree->Branch("chargec6",          &chargec6);
  ftree->Branch("chargec7",          &chargec7);

  ftree->Branch("chargeg1",          &chargeg1);
  ftree->Branch("chargeg2",          &chargeg2);
  ftree->Branch("chargeg3",          &chargeg3);
  ftree->Branch("chargeg4",          &chargeg4);
  ftree->Branch("chargeg5",          &chargeg5);
  ftree->Branch("chargeg6",          &chargeg6);
  ftree->Branch("chargeg7",          &chargeg7);

/*
  ftree->Branch("ieta1",             &ieta1);
  ftree->Branch("ieta2",             &ieta2);
  ftree->Branch("ieta3",             &ieta3);
  ftree->Branch("ieta4",             &ieta4);
  ftree->Branch("ieta5",             &ieta5);
  ftree->Branch("ieta6",             &ieta6);
  ftree->Branch("ieta7",             &ieta7);

  ftree->Branch("iphi1",             &iphi1);
  ftree->Branch("iphi2",             &iphi2);
  ftree->Branch("iphi3",             &iphi3);
  ftree->Branch("iphi4",             &iphi4);
  ftree->Branch("iphi5",             &iphi5);
  ftree->Branch("iphi6",             &iphi6);
  ftree->Branch("iphi7",             &iphi7);
*/
  ftree->Branch("cross1",            &cross1);
  ftree->Branch("cross2",            &cross2);
  ftree->Branch("cross3",            &cross3);
  ftree->Branch("cross4",            &cross4);
  ftree->Branch("cross5",            &cross5);
  ftree->Branch("cross6",            &cross6);
  ftree->Branch("cross7",            &cross7);

  ftree->Branch("gain1",             &gain1);
  ftree->Branch("gain2",             &gain2);
  ftree->Branch("gain3",             &gain3);
  ftree->Branch("gain4",             &gain4);
  ftree->Branch("gain5",             &gain5);
  ftree->Branch("gain6",             &gain6);
  ftree->Branch("gain7",             &gain7);
/*
  ftree->Branch("resp1",             &resp1);
  ftree->Branch("resp2",             &resp2);
  ftree->Branch("resp3",             &resp3);
  ftree->Branch("resp4",             &resp4);
  ftree->Branch("resp5",             &resp5);
  ftree->Branch("resp6",             &resp6);
  ftree->Branch("resp7",             &resp7);
*/

}

// ------------ method called once each job just after ending the event loop  ------------

// ------------ method called when starting to processes a run  ------------
void HBHEMuonCosmic2::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {

  edm::ESHandle<HcalDDDRecConstants> pHRNDC;
  iSetup.get<HcalRecNumberingRecord>().get(pHRNDC);
  const HcalDDDRecConstants & hdc = (*pHRNDC);
  actHB.clear();
  actHE.clear();
  actHB = hdc.getThickActive(0);
  actHE = hdc.getThickActive(1);
  cout<<" active l "<<actHE.size()<<endl;  
  cout<<" active l "<<actHB.size()<<endl;  
  cout<<"getRZ "<<hdc.getRZ(1,1,1,1)<<" "<<hdc.getRZ(1,5,1,1)<<" "<<hdc.getRZ(1,1,1,2)<<" "<<hdc.getRZ(1,1,1,3)
  <<" "<<hdc.getRZ(1,1,1,4)<<endl;
  rmin = hdc.getRZ(1,1,1,1);
  cout<<"getPhiBin "<<hdc.getPhiBin(0)<<" "<<hdc.getPhiBin(1)<<"/"<<hdc.getPhiBin(2)<<endl;
  cout<<"getNPhi "<<hdc.getNPhi(1)<<endl;
  std::pair<double, double> pp = hdc.getEtaPhi(1,1,1);
  cout<<"eta phi 1 1 "<<pp.first<<"/"<<pp.second<<endl;
  cout<<"phi 1/2 "<<hdc.getPhiOff(1)<<"/"<<hdc.getPhiOff(2)<<endl;
  cout<<"depth 1 layerf/layerb "<<hdc.getLayerFront(HcalBarrel, 1, 1, 1)<<"/"<<hdc.getLayerBack(HcalBarrel, 1, 1, 1)<<endl;
  cout<<"depth 2 layerf/layerb "<<hdc.getLayerFront(HcalBarrel, 1, 1, 2)<<"/"<<hdc.getLayerBack(HcalBarrel, 1, 1, 2)<<endl;
  cout<<"depth 3 layerf/layerb "<<hdc.getLayerFront(HcalBarrel, 1, 1, 3)<<"/"<<hdc.getLayerBack(HcalBarrel, 1, 1, 3)<<endl;
  cout<<"depth 4 layerf/layerb "<<hdc.getLayerFront(HcalBarrel, 1, 1, 4)<<"/"<<hdc.getLayerBack(HcalBarrel, 1, 1, 4)<<endl;
  rz_depth1 = hdc.getRZ(1,1,1,1);
  rz_depth2 = hdc.getRZ(1,1,1,2);
  rz_depth3 = hdc.getRZ(1,1,1,3);
  rz_depth4 = hdc.getRZ(1,1,1,4);
  HcalDetId myId1(HcalBarrel,1,1,1);
  std::pair<double, double> myt1 = hdc.getRZ(myId1);
  cout<<"HcalDetId depth1 first/second "<<myt1.first<<"/"<<myt1.second<<endl;
  HcalDetId myId2(HcalBarrel,1,1,2);
  std::pair<double, double> myt2 = hdc.getRZ(myId2);
  cout<<"HcalDetId depth2 first/second "<<myt2.first<<"/"<<myt2.second<<endl;
  HcalDetId myId3(HcalBarrel,1,1,3);
  std::pair<double, double> myt3 = hdc.getRZ(myId3);
  cout<<"HcalDetId depth3 first/second "<<myt3.first<<"/"<<myt3.second<<endl;
  HcalDetId myId4(HcalBarrel,1,1,4);
  std::pair<double, double> myt4 = hdc.getRZ(myId4);
  cout<<"HcalDetId depth4 first/second "<<myt4.first<<"/"<<myt4.second<<endl;
  HcalDetId myId5(HcalBarrel,1,1,5);
  std::pair<double, double> myt5 = hdc.getRZ(myId5);
  cout<<"HcalDetId depth5 first/second  "<<myt5.first<<"/"<<myt5.second<<endl;
//  ndepth2 = 4;
//  ndepth3 = 5;
//  ndepth4 = 7;
  depth1_start = myt1.first; 
  depth1_end = myt1.second; 

  double met2 = (myt2.second-myt2.first-ndepth2*0.37)/(ndepth2-1);
  for ( int i = 0; i < ndepth2; i++ )
  {
   depth2_start[i] = myt2.first + (0.37+met2)*i; 
   depth2_end[i] = depth2_start[i] + 0.37;
   cout<<"depth2 i/start/end "<<i<<"/"<<depth2_start[i]<<"/"<<depth2_end[i]<<endl; 
  }  

  double met3 = (myt3.second-myt3.first-ndepth3*0.37)/(ndepth3-1);
  for ( int i = 0; i < ndepth3; i++ )
  {
   depth3_start[i] = myt3.first + (0.37+met3)*i; 
   depth3_end[i] = depth3_start[i] + 0.37;
   cout<<"depth3 i/start/end "<<i<<"/"<<depth3_start[i]<<"/"<<depth3_end[i]<<endl; 
  }  
//  depth3_start = myt3.first; 
//  depth3_end = myt3.second; 

  double met4 = (myt4.second-myt4.first-ndepth4*0.37)/(ndepth4-1);
  for ( int i = 0; i < ndepth4; i++ )
  {
   depth4_start[i] = myt4.first + (0.37+met4)*i; 
   depth4_end[i] = depth4_start[i] + 0.37;
   cout<<"depth4 i/start/end "<<i<<"/"<<depth4_start[i]<<"/"<<depth4_end[i]<<endl; 
  }  
//  depth4_start = myt4.first; 
//  depth4_end = myt4.second; 

  cout<<"met2/met3/met4 "<<met2<<"/"<<met3<<"/"<<met4<<endl;
  for (unsigned int i=0; i<actHE.size(); ++i) 
  {
   cout<<"actHE "<<i<<" "<<actHE[i].ieta<<" "<<actHE[i].depth<<endl;
  }
  for (unsigned int i=0; i<actHB.size(); ++i) 
  {
   cout<<"actHB "<<i<<" "<<actHB[i].ieta<<" "<<actHB[i].depth<<endl;
  }

 
  bool changed = true;
  all_triggers.clear();
  if (hltConfig_.init(iRun, iSetup,"HLT" , changed)) {
    // if init returns TRUE, initialisation has succeeded!
#ifdef EDM_ML_DEBUG
    edm::LogInfo("HBHEMuon") << "HLT config with process name " 
			     << "HLT" << " successfully extracted"
			     << std::endl;
#endif
//  std::string string_search[5]={"HLT_IsoMu_","HLT_L1SingleMu_","HLT_L2Mu","HLT_Mu","HLT_RelIso1p0Mu"};
    std::string string_search[6]={"HLT_L1SingleMuOpen_v2","HLT_L1SingleMuOpen_DT_v2","HLT_IsoMu24","HLT_IsoMu27","HLT_Mu45","HLT_Mu50"};
  
    unsigned int ntriggers = hltConfig_.size();
    for (unsigned int t=0;t<ntriggers;++t) {
      std::string hltname(hltConfig_.triggerName(t));
      for (unsigned int ik=0; ik<6; ++ik) {
	if (hltname.find(string_search[ik])!=std::string::npos ){
          cout<<"find trigger "<<hltname<<endl;
	  all_triggers.push_back(hltname);
	  break;
	}
      }
    }//loop over ntriggers
    edm::LogInfo("HBHEMuon") << "All triggers size in begin run " 
			     << all_triggers.size() << std::endl;
  } else {
    edm::LogError("HBHEMuon") << "Error! HLT config extraction with process name " 
			      << "HLT" << " failed";
  }

  edm::ESHandle<HcalTopology> htopo;
  iSetup.get<HcalRecNumberingRecord>().get(htopo);
  theHBHETopology = htopo.product();

  edm::ESHandle<HcalRespCorrs> resp;
  iSetup.get<HcalRespCorrsRcd>().get(resp);
  respCorrs_ = new HcalRespCorrs(*resp.product());
  respCorrs_->setTopo(theHBHETopology);
  
}


// ------------ method called when ending the processing of a run  ------------
void HBHEMuonCosmic2::endRun(edm::Run const&, edm::EventSetup const&) { }

// ------------ method called when starting to processes a luminosity block  ------------
void HBHEMuonCosmic2::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { }

// ------------ method called when ending the processing of a luminosity block  ------------
void HBHEMuonCosmic2::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { }

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HBHEMuonCosmic2::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("HLTriggerResults",edm::InputTag("TriggerResults","","HLT"));
  //desc.add<std::string>("LabelBeamSpot","offlineBeamSpot");
  desc.add<bool>("UseRaw",false);
  desc.add<std::string>("LabelVertex","offlinePrimaryVertices");
  desc.add<edm::InputTag>("LabelEBRecHit",edm::InputTag("ecalRecHit","EcalRecHitsEB"));
  desc.add<edm::InputTag>("LabelEERecHit",edm::InputTag("ecalRecHit","EcalRecHitsEE"));
//  desc.add<edm::InputTag>("LabelHBHERecHit",edm::InputTag("hbheprereco"));
  desc.add<edm::InputTag>("LabelHBHERecHit",edm::InputTag("hbhereco"));
  desc.add<edm::InputTag>("HLTriggerEvent", edm::InputTag("hltTriggerSummaryAOD", "", "HLT"));
  desc.add<edm::InputTag>("LabelHFRecHit",edm::InputTag("hfreco"));
  desc.add<std::string>("LabelMuon","muons");
//  desc.add<edm::InputTag>("LabelTracks",edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("LabelTracks",edm::InputTag("cosmicMuons"));
  desc.addUntracked<int>("Verbosity",0);
//  desc.addUntracked<int>("MaxDepth",4);
  desc.addUntracked<int>("MaxDepth",7);
  descriptions.add("HBHEMuonCosmic2",desc);
  //desc.setUnknown();
  //descriptions.addDefault(desc);
}

void HBHEMuonCosmic2::clearVectors() {
  ///clearing vectots
}

int HBHEMuonCosmic2::matchId(const HcalDetId& id1, const HcalDetId& id2) {

  HcalDetId kd1(id1.subdet(),id1.ieta(),id1.iphi(),1);
  HcalDetId kd2(id2.subdet(),id2.ieta(),id2.iphi(),1);
  int match = ((kd1 == kd2) ? 1 : 0);
  return match;
}

double HBHEMuonCosmic2::activeLength(const DetId& id_) {
  HcalDetId id(id_);
  int ieta = id.ietaAbs();
  int zside= id.zside();
  int iphi = id.iphi();
  std::vector<int> dpths;
  /*
  if (mergedDepth_) {
    std::vector<HcalDetId> ids;
    hdc_->unmergeDepthDetId(id,ids);
    for (auto idh : ids) 
      dpths.emplace_back(idh.depth());
  } else {
  */
    dpths.emplace_back(id.depth());
  //}
  double lx(0);
  if (id.subdet() == HcalBarrel) {
    cout<<"active len/ieta/iphi "<<actHB.size()<<" "<<ieta<<" "<<iphi<<endl;
    for (unsigned int i=0; i<actHB.size(); ++i) {
      if ((ieta == actHB[i].ieta) && (zside == actHB[i].zside) &&
          (std::find(dpths.begin(),dpths.end(),actHB[i].depth) != dpths.end())&&
          (std::find(actHB[i].iphis.begin(),actHB[i].iphis.end(),iphi) !=
           actHB[i].iphis.end())) {
        cout<<"i thick/lx/ieta/iphi "<<actHB[i].thick<<" "<<lx<<" "<<actHB[i].ieta<<endl;
        lx += actHB[i].thick;
      }
    }
  } else {
    for (unsigned int i=0; i<actHE.size(); ++i) {
      if ((ieta == actHE[i].ieta) && (zside == actHE[i].zside) && 
          (std::find(dpths.begin(),dpths.end(),actHE[i].depth) != dpths.end())&&
          (std::find(actHE[i].iphis.begin(),actHE[i].iphis.end(),iphi) !=
           actHE[i].iphis.end())) 
      {
//        cout<<" ll "<<i<<endl;
        lx += actHE[i].thick;
      }
    }
  }
  return lx;
}

double HBHEMuonCosmic2::respCorr(const DetId& id) {
  double cfac(1.0);
  /*
  if (useMyCorr_) {
    auto itr = corrValue_.find(id);
    if (itr != corrValue_.end()) cfac = itr->second;
  } else*/
  if (respCorrs_ != nullptr) {
    cfac = (respCorrs_->getValues(id))->getValue();
  }
  return cfac;
}

double HBHEMuonCosmic2::gainFactor(const edm::ESHandle<HcalDbService>& conditions, const HcalDetId& id) {
  double gain(0.0);
  const HcalCalibrations& calibs=conditions->getHcalCalibrations(id);
  for (int capid=0; capid<4; ++capid)
    gain += (0.25*calibs.respcorrgain(capid));
  return gain;
}

double HBHEMuonCosmic2::pedWidth(const edm::ESHandle<HcalDbService>& conditions, const HcalDetId& id) {
  double ped(0.0);
  const HcalPedestalWidth* calibs=conditions->getPedestalWidth(id);
  for (int capid=0; capid<4; ++capid)
    ped += (0.25*calibs->getWidth(capid));
  return ped;
}

double HBHEMuonCosmic2::pedEffective(const edm::ESHandle<HcalDbService>& conditions, const HcalDetId& id) {
  double ped(0.0);
  const HcalPedestal* calibs=conditions->getEffectivePedestal(id);
  for (int capid=0; capid<4; ++capid)
    ped += (0.25*calibs->getValue(capid));
  return ped;
}

double HBHEMuonCosmic2::pedestal(const edm::ESHandle<HcalDbService>& conditions, const HcalDetId& id) {
  double ped(0.0);
  const HcalPedestal* calibs=conditions->getPedestal(id);
  for (int capid=0; capid<4; ++capid)
    ped += (0.25*calibs->getValue(capid));
  return ped;
}

double HBHEMuonCosmic2::darkCurrent(const edm::ESHandle<HcalDbService>& conditions, const HcalDetId& id) {
  float ped(0.0);
  const HcalSiPMParameter* calibs=conditions->getHcalSiPMParameter(id);
  ped = calibs->getDarkCurrent();
  return ped;
}

double HBHEMuonCosmic2::FCbyPE(const edm::ESHandle<HcalDbService>& conditions, const HcalDetId& id) {
  float ped(0.0);
  const HcalSiPMParameter* calibs=conditions->getHcalSiPMParameter(id);
  ped = calibs->getFCByPE();
  return ped;
}

int HBHEMuonCosmic2::depth16HE(int ieta, int iphi) {

  int depth = (maxDepth_ <= 6) ? 3 : 4;
//  if ((ieta < 0) || (iphi < 62) || (iphi > 66)) depth = 3;
  return depth;
}


/*
double HBHEMuonCosmic2::activeLength(const DetId& id_) {
  HcalDetId id(id_);
  int ieta = id.ietaAbs();
  int depth= id.depth();
  double lx(0);
  if (id.subdet() == HcalBarrel) {
    for (unsigned int i=0; i<actHB.size(); ++i) {
      if (ieta == actHB[i].ieta && depth == actHB[i].depth) {
	lx = actHB[i].thick;
	break;
      }
    }
  } else {
    for (unsigned int i=0; i<actHE.size(); ++i) {
      if (ieta == actHE[i].ieta && depth == actHE[i].depth) {
	lx = actHE[i].thick;
	break;
      }
    }
  }
  return lx;
}
*/
bool HBHEMuonCosmic2::isGoodVertex(const reco::Vertex& vtx) {
  if (vtx.isFake())                   return false;
  if (vtx.ndof() < 4)                 return false;
  if (vtx.position().Rho() > 2.)      return false;
  if (fabs(vtx.position().Z()) > 24.) return false;
  return true;
}

void HBHEMuonCosmic2::endJob() {
  out_file->cd();
  cout<<" out file "<<endl;
  gDirectory->ls();

  ftree->Write();
//  ftree1->Write();
  ftree->Print();
  out_file->Close();

}

spr::propagatedTrackID HBHEMuonCosmic2::propagateCALO1(const reco::Track* pTrack, const CaloGeometry* geo, const MagneticField* bField, bool debug) 
{

    const CaloSubdetectorGeometry *barrelGeom = geo->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);
    const CaloSubdetectorGeometry *endcapGeom = geo->getSubdetectorGeometry(DetId::Ecal,EcalEndcap);
    const CaloSubdetectorGeometry* gHB =        geo->getSubdetectorGeometry(DetId::Hcal,HcalBarrel);

    spr::propagatedTrackID vdet;
    vdet.ok        = true;
    vdet.detIdECAL = DetId(0);
    vdet.detIdHCAL = DetId(0);
    vdet.detIdEHCAL= DetId(0);
#ifdef EDM_ML_DEBUG
    if (debug) std::cout << "Propagate track:  p " << pTrack->p() << " eta " << pTrack->eta() << " phi " << pTrack->phi() << " Flag  " << vdet.ok << std::endl;
#endif
    std::pair<math::XYZPoint,bool> info = spr::propagateECAL (pTrack, bField, debug);
    vdet.okECAL = info.second;
    if (vdet.okECAL) {
      const GlobalPoint point(info.first.x(),info.first.y(),info.first.z());
      vdet.etaECAL = point.eta();
      vdet.phiECAL = point.phi();
      if (std::abs(point.eta())<spr::etaBEEcal) {
        vdet.detIdECAL = barrelGeom->getClosestCell(point);
      } else {
        if (endcapGeom)
          vdet.detIdECAL = endcapGeom->getClosestCell(point);
        else
          vdet.okECAL    = false;
      }
      vdet.detIdEHCAL = gHB->getClosestCell(point);
    }
    info = spr::propagateHCAL (pTrack, bField, debug);
    vdet.okHCAL = info.second;
    if (vdet.okHCAL) {
      const GlobalPoint point(info.first.x(),info.first.y(),info.first.z());
      vdet.etaHCAL = point.eta();
      vdet.phiHCAL = point.phi();
      vdet.detIdHCAL = gHB->getClosestCell(point);
    }
#ifdef EDM_ML_DEBUG
    if (debug) {
      std::cout << "propagateCALO:: for 1 track" << std::endl;
      std::cout << "Track [0] Flag: " << vdet.ok << " ECAL (" << vdet.okECAL << ") ";
      if (vdet.detIdECAL.subdetId() == EcalBarrel) {
        std::cout << (EBDetId)(vdet.detIdECAL);
      } else {
        std::cout << (EEDetId)(vdet.detIdECAL);
      }
      std::cout << " HCAL (" << vdet.okHCAL << ") " << (HcalDetId)(vdet.detIdHCAL) << " Or " << (HcalDetId)(vdet.detIdEHCAL) << std::endl;
    }
#endif
    return vdet;
}

std::pair<bool,HcalDetId> HBHEMuonCosmic2::propagateHCALBack1(const reco::Track* track, const CaloGeometry* geo, const MagneticField* bField, bool debug) 
{
    const CaloSubdetectorGeometry* gHB = geo->getSubdetectorGeometry(DetId::Hcal,HcalBarrel);
    const GlobalPoint  vertex (track->vx(), track->vy(), track->vz());
    const GlobalVector momentum (track->px(), track->py(), track->pz());
    int charge (track->charge());
    spr::propagatedTrack info = propagateCalo2(vertex, momentum, charge, bField, spr::zBackHE, spr::rBackHB, spr::etaBEHcal, debug);
    if (info.ok) {
      const GlobalPoint point = GlobalPoint(info.point.x(),info.point.y(),info.point.z());
      return std::pair<bool,HcalDetId>(true,HcalDetId(gHB->getClosestCell(point)));
    } else {
      return std::pair<bool,HcalDetId>(false,HcalDetId());
    }
}



spr::propagatedTrack HBHEMuonCosmic2::propagateCalo2(const GlobalPoint& tpVertex,
                                     const GlobalVector& tpMomentum,
                                     int tpCharge, const MagneticField* bField,
                                     float zdist, float radius, float corner, bool debug ) 
{

    spr::propagatedTrack track;
#ifdef EDM_ML_DEBUG
    if (debug) std::cout << "propagateCalo:: Vertex " << tpVertex << " Momentum " << tpMomentum << " Charge " << tpCharge << " Radius " << radius << " Z " << zdist << " Corner " << corner << std::endl;
#endif
    FreeTrajectoryState fts (tpVertex, tpMomentum, tpCharge, bField);

    Plane::PlanePointer lendcap = Plane::build(Plane::PositionType (0, 0, -zdist), Plane::RotationType());
    Plane::PlanePointer rendcap = Plane::build(Plane::PositionType (0, 0,  zdist), Plane::RotationType());

    Cylinder::CylinderPointer barrel = Cylinder::build(Cylinder::PositionType (0, 0, 0), Cylinder::RotationType (), radius);

    AnalyticalPropagator myAP (bField, alongMomentum, 2*M_PI);

    TrajectoryStateOnSurface tsose;
    if (tpMomentum.eta() < 0) {
      tsose = myAP.propagate(fts, *lendcap);
    } else {
      tsose = myAP.propagate(fts, *rendcap);
    }

    TrajectoryStateOnSurface tsosb = myAP.propagate(fts, *barrel);

    track.ok=true;
    if (tsose.isValid() && tsosb.isValid()) {
      float absEta = std::abs(tsosb.globalPosition().eta());
      if (absEta < corner) {
        track.point.SetXYZ(tsosb.globalPosition().x(), tsosb.globalPosition().y(), tsosb.globalPosition().z());
        track.direction = tsosb.globalDirection();
      } else {
        track.point.SetXYZ(tsose.globalPosition().x(), tsose.globalPosition().y(), tsose.globalPosition().z());
        track.direction = tsose.globalDirection();
      }
    } else if (tsose.isValid()) {
      track.point.SetXYZ(tsose.globalPosition().x(), tsose.globalPosition().y(), tsose.globalPosition().z());
      track.direction = tsose.globalDirection();
    } else if (tsosb.isValid()) {
      track.point.SetXYZ(tsosb.globalPosition().x(), tsosb.globalPosition().y(), tsosb.globalPosition().z());
      track.direction = tsosb.globalDirection();
    } else {
      track.point.SetXYZ(-999., -999., -999.);
      track.direction = GlobalVector(0,0,1);
      track.ok = false;
    }
#ifdef EDM_ML_DEBUG
    if (debug) {
      std::cout << "propagateCalo:: Barrel " << tsosb.isValid() << " Endcap " << tsose.isValid() << " OverAll " << track.ok << " Point " << track.point << " Direction " << track.direction << std::endl;
      if (track.ok) {
        math::XYZPoint vDiff(track.point.x()-tpVertex.x(), track.point.y()-tpVertex.y(), track.point.z()-tpVertex.z());
        double dphi = track.direction.phi()-tpMomentum.phi();
        double rdist = std::sqrt(vDiff.x()*vDiff.x()+vDiff.y()*vDiff.y());
        double pt    = tpMomentum.perp();
        double rat   = 0.5*dphi/std::sin(0.5*dphi);
        std::cout << "RDist " << rdist << " pt " << pt << " r/pt " << rdist*rat/pt << " zdist " << vDiff.z() << " pz " << tpMomentum.z() << " z/pz " << vDiff.z()/tpMomentum.z() << std::endl;
      }
    }
#endif
    return track;
}


void HBHEMuonCosmic2::energyHCALCellmy(HcalDetId detId, edm::Handle<HFRecHitCollection>& hits, 
              std::vector<std::pair<double,int> >& energyCell, 
               double hbThr, double heThr, double hfThr, 
               double hoThr, double tMin, double tMax, bool
 #ifdef EDM_ML_DEBUG
  debug
 #endif
     ) 
{
// cout<<"inside "<<endl;
 energyCell.clear();
 int    subdet  = detId.subdet();    
 double eThr    = spr::eHCALThreshold(subdet, hbThr, heThr, hfThr, hoThr);
// bool   hbhe    = (detId.ietaAbs() == 16);
 #ifdef EDM_ML_DEBUG
  if (debug)
   std::cout << "energyHCALCell: input ID " << detId <<  " Threshold (E) " << eThr << " (T) " << tMin << ":" << tMax << std::endl;
 #endif
  HcalSubdetector subdet0 = detId.subdet();
  HcalDetId hcid(subdet0,detId.ieta(),detId.iphi(),1);
  DetId det(hcid.rawId());

  double energy(0);
  int ii = -10;
  if (hits.isValid())
  {
   for (HFRecHitCollection::const_iterator j=hits->begin(); j != hits->end(); j++)
   {
    HcalDetId cell(j->id());
    if (j->time() > tMin && j->time() < tMax && cell.ieta() == detId.ieta() && cell.iphi() == detId.iphi()) 
    {
//     cout<<"energy "<<j->energy()<<" "<<j->time()<<" "<<cell.ieta()<<" "<<cell.iphi()<<endl;
     energy = j->energy();
     ii = cell.depth();
    #ifdef EDM_ML_DEBUG
     if (debug)
      std::cout << "energyHCALCell:: Hit " << hcid << " E " << j->energy() << " t " << j->time() << std::endl;
    #endif
    #ifdef EDM_ML_DEBUG
     if (debug)
      std::cout << "energyHCALCell:: Cell " << hcid << " E " << energy << " threshold " << eThr << std::endl;
    #endif
     if (energy > eThr && ii >= 0 ) 
     {
      energyCell.push_back(std::pair<double,int>(energy,ii));
     }
    }
   }
  }
// }
 #ifdef EDM_ML_DEBUG
  if (debug) 
  {
   std::cout << "energyHCALCell:: " << energyCell.size();
   for (unsigned int i=0; i<energyCell.size(); ++i) 
   {
    std::cout << " [" << i << "] (" << energyCell[i].first << ":"
                     << energyCell[i].second << ")";
   }
   std::cout << std::endl;
  }
 #endif
}
 
GlobalPoint HBHEMuonCosmic2::myfind(const reco::Track* pTrack, double depth)
//int HBHEMuonCosmic2::myfind(const reco::Track* pTrack, double depth)
{
// double m = pTrack->px();
// double n = pTrack->py();
// double p = pTrack->pz();
 double m = sin(pTrack->theta())*cos(pTrack->phi());
 double n = sin(pTrack->theta())*sin(pTrack->phi());
 double p = cos(pTrack->theta());
 double x0 = pTrack->vx();
 double y0 = pTrack->vy();
 double z0 = pTrack->vz();
 double rz_min = 1000000;
// int i_min = 0;
 double xk, yk, zk;
 double xk1, yk1, zk1;
 for ( int i = 0; i < 72; i++)
 {
  double angle = 2.5*TMath::Pi()/180. + i*5.*TMath::Pi()/180.;
  double a = depth/cos(angle);
  double b = depth/cos(TMath::Pi()/2.-angle);
  double A = b;
  double B = a;
  double D = -a*b;
//  cout<<"A/B/D "<<A<<"/"<<B<<"/"<<D<<endl;
  double t0 = -(A*x0+B*y0+D)/(A*m+B*n);
//  if ( t0 < 0 ) t0 = -t0;
  xk = m*t0+x0;
  yk = n*t0+y0;
  zk = p*t0+z0;
//  if ( sqrt(xk*xk+yk*yk+zk*zk) < rz_min && ((xk-x0)/m > 0 && (yk-y0)/n > 0 && (zk-z0)/p > 0) ) 
  if ( sqrt((xk-x0)*(xk-x0)+(yk-y0)*(yk-y0)+(zk-z0)*(zk-z0)) < rz_min && 
     ((xk-x0)/m > 0 && (yk-y0)/n > 0 && (zk-z0)/p > 0) ) 
  {
   rz_min = sqrt((xk-x0)*(xk-x0)+(yk-y0)*(yk-y0)+(zk-z0)*(zk-z0)); 
//   i_min=i; 
   xk1=xk; yk1=yk; zk1=zk;
  }
//  cout<<"i/x/y/z/r/m/n/p/A/B/D/x0/y0/z0 "<<i<<"/"<<xk<<"/"<<yk<<"/"<<zk<<" "
//  <<sqrt((xk-x0)*(xk-x0)+(yk-y0)*(yk-y0)+(zk-z0)*(zk-z0))<<"/"<<rz_min<<endl;
 }
// cout<<"i_min/rz_min/depth/z0 "<<i_min<<"/"<<rz_min<<"/"<<depth<<"/"<<z0<<endl;
 GlobalPoint myp(xk1, yk1, zk1);
 return myp;
// return i_min;
}

double HBHEMuonCosmic2::mycorrection(const reco::Track* pTrack, int ieta, int iphi)
{
 double phi1 = (iphi-1)*5.*TMath::Pi()/180 + 2.5*TMath::Pi()/180;
 double theta1;
 if ( ieta > 0 )
 {
  theta1 = 0.087/2. + (ieta-1)*0.087;
  theta1 = exp(-theta1);
  theta1 = 2.*atan(theta1);
 }
 else
 {
  theta1 = -0.087/2. - (fabs(ieta)-1)*0.087;
  theta1 = exp(-theta1);
  theta1 = 2.*atan(theta1);
 }
 theta1 = TMath::Pi()/2.;
 double v2x, v2y, v2z;
 v2x = sin(theta1)*cos(phi1);
 v2y = sin(theta1)*sin(phi1);
 v2z = cos(theta1);

 double phi = pTrack->phi();
 double theta = pTrack->theta();
 double v1x, v1y, v1z;
// if ( phi > 0 )
// {
  v1x = sin(theta)*cos(phi);
  v1y = sin(theta)*sin(phi);
  v1z = cos(theta);
// }
// cout<<"iphi phi/phi1/theta "<<iphi<<"/"<<phi*180/TMath::Pi()<<"/"<<phi1*180/TMath::Pi()<<"/"<<theta*180/TMath::Pi()<<endl;
// else
// {
//  v1x = sin(theta)*cos(phi+2.*TMath::Pi());
//  v1y = sin(theta)*sin(phi+2.*TMath::Pi());
//  v1z = cos(theta);
// }
 double cross = v1x*v2x+v1y*v2y+v1z*v2z;
// double norm1 =sqrt(v1x*v1x+v1y*v1y+v1z*v1z);
// double norm2 =sqrt(v2x*v2x+v2y*v2y+v2z*v2z);
// double mycorr = cross/norm1/norm2;
 double mycorr = cross;
// if ( mycorr < 0.5 )
//  cout<<"--------- cross phi1/phi/theta/v2x/v2y/v1x/v1y "<<cross<<"/"<<phi1<<"/"<<phi<<"/"<<theta
//  <<"/"<<v2x<<"/"<<v2y<<"/"<<v1x<<"/"<<v1y<<endl;

 return fabs(mycorr);
}

double HBHEMuonCosmic2::mycorrection1(const reco::Track* pTrack, int ieta, int iphi)
{
// ieta = -ieta;
 double phi1 = (iphi-1)*5.*TMath::Pi()/180 + 2.5*TMath::Pi()/180;
 double v2x, v2y, v2z;
 v2x = cos(phi1);
 v2y = sin(phi1);
 v2z = 0;

 double phi = pTrack->phi();
 double theta = pTrack->theta();
 double v1x, v1y, v1z;
// if ( phi > 0 )
// {
  v1x = sin(theta)*cos(phi);
  v1y = sin(theta)*sin(phi);
  v1z = cos(theta);
// }
// cout<<"iphi phi/phi1/theta "<<iphi<<"/"<<phi*180/TMath::Pi()<<"/"<<phi1*180/TMath::Pi()<<"/"<<theta*180/TMath::Pi()<<endl;
// else
// {
//  v1x = sin(theta)*cos(phi+2.*TMath::Pi());
//  v1y = sin(theta)*sin(phi+2.*TMath::Pi());
//  v1z = cos(theta);
// }
 double cross = v1x*v2x+v1y*v2y+v1z*v2z;
// double norm1 =sqrt(v1x*v1x+v1y*v1y+v1z*v1z);
// double norm2 =sqrt(v2x*v2x+v2y*v2y+v2z*v2z);
// double mycorr = cross/norm1/norm2;
 double mycorr = cross;
// cout<<"cross/n1/n2 "<<cross<<"/"<<norm1<<"/"<<norm2<<endl;

 return fabs(mycorr);
}



std::pair<double, double> HBHEMuonCosmic2::mycheckz(int ieta, double depth)
{
 double theta1, theta2;
 double etamin, etamax;
 if ( ieta > 0 ) 
 {
  etamax = 0.087+(ieta-1)*0.087;
  etamin = (ieta-1)*0.087;
 }
 if ( ieta < 0 ) 
 {
  etamax = -0.087-(fabs(ieta)-1)*0.087;
  etamin = -(fabs(ieta)-1)*0.087;
 }
 theta1 = exp(-etamin); 
 theta2 = exp(-etamax);
 theta1 = 2.*atan(theta1); 
 theta2 = 2.*atan(theta2); 
// cout<<"theta1/theta2/pi "<<theta1<<"/"<<theta2<<"/"<<TMath::Pi()/2<<endl;
 double zmin, zmax;
 if ( ieta > 0 )
 {
  zmin = depth/tan(theta1);
  zmax = depth/tan(theta2);
 }
 else
 {
  zmax = depth/tan(theta1);
  zmin = depth/tan(theta2);
 }
 std::pair<double, double> par(zmin, zmax);
 return par;
}

std::pair<double, double> HBHEMuonCosmic2::mycheckx(int iphi, double depth)
{
 double phimin, phimax;
 phimax = iphi*5.*TMath::Pi()/180.;
 phimin = (iphi-1)*5.*TMath::Pi()/180.;
 double xmin, xmax;
 if ( phimin < TMath::Pi() )
 {
  xmax = depth*cos(phimin);
  xmin = depth*cos(phimax);
 }
 else
 {
  xmax = depth*cos(phimax);
  xmin = depth*cos(phimin);
 }
// cout<<"depth/phimin/phimax "<<depth<<"/"<<phimin<<"/"<<phimax<<endl;
 std::pair<double, double> par(xmin, xmax);
 return par;
}

std::pair<double, double> HBHEMuonCosmic2::mychecky(int iphi, double depth)
{
 double phimin, phimax;
 phimax = iphi*5.*TMath::Pi()/180.;
 phimin = (iphi-1)*5.*TMath::Pi()/180.;
 double xmin, xmax;
 if ( phimin < TMath::Pi()/2 || phimin > TMath::Pi()*3/2 )
 {
  xmax = depth*sin(phimax);
  xmin = depth*sin(phimin);
 }
 else
 {
  xmax = depth*sin(phimin);
  xmin = depth*sin(phimax);
 }
 if ( phimin > phimax )
 {
  xmax = depth*sin(phimax);
  xmin = depth*sin(phimin);
 }
// cout<<"depth/phimin/phimax "<<depth<<"/"<<phimin<<"/"<<phimax<<endl;
 std::pair<double, double> par(xmin, xmax);
 return par;
}


int HBHEMuonCosmic2::myfind_iphi(const reco::Track* pTrack, double depth)
{
 double m = pTrack->px();
 double n = pTrack->py();
 double p = pTrack->pz();
 double x0 = pTrack->vx();
 double y0 = pTrack->vy();
 double z0 = pTrack->vz();
 double rz_min = 1000000;
 int i_min = 0;
 double xk, yk, zk;
// double xk1, yk1, zk1;
 for ( int i = 1; i <= 72; i++)
 {
  double angle = 2.5*TMath::Pi()/180. + (i-1)*5.*TMath::Pi()/180.;
  double a = depth/cos(angle);
  double b = depth/cos(TMath::Pi()/2.-angle);
  double A = b;
  double B = a;
  double D = -a*b;
//  cout<<"A/B/D "<<A<<"/"<<B<<"/"<<D<<endl;
  double t0 = -(A*x0+B*y0+D)/(A*m+B*n);
//  if ( t0 < 0 ) t0 = -t0;
  xk = m*t0+x0;
  yk = n*t0+y0;
  zk = p*t0+z0;
//  if ( sqrt(xk*xk+yk*yk+zk*zk) < rz_min && ((xk-x0)/m > 0 && (yk-y0)/n > 0 && (zk-z0)/p > 0) ) 
  if ( sqrt((xk-x0)*(xk-x0)+(yk-y0)*(yk-y0)+(zk-z0)*(zk-z0)) < rz_min && 
     ((xk-x0)/m > 0 && (yk-y0)/n > 0 && (zk-z0)/p > 0) ) 
  {
   rz_min = sqrt((xk-x0)*(xk-x0)+(yk-y0)*(yk-y0)+(zk-z0)*(zk-z0)); 
   i_min=i; 
//   xk1=xk; yk1=yk; zk1=zk;
  }
//  cout<<"i/x/y/z/r/m/n/p/A/B/D/x0/y0/z0 "<<i<<"/"<<xk<<"/"<<yk<<"/"<<zk<<" "
//  <<sqrt((xk-x0)*(xk-x0)+(yk-y0)*(yk-y0)+(zk-z0)*(zk-z0))<<"/"<<rz_min<<endl;
 }
// cout<<"i_min/rz_min/depth/z0 "<<i_min<<"/"<<rz_min<<"/"<<depth<<"/"<<z0<<endl;
// GlobalPoint myp(xk1, yk1, zk1);
// return myp;
 return i_min;
}

int HBHEMuonCosmic2::myfind_ieta(GlobalPoint p)
{
// double theta1, theta2;
 double etamin, etamax;
// if ( ieta > 0 ) 
// {
//  etamax = 0.087+(ieta-1)*0.087;
//  etamin = (ieta-1)*0.087;
// }
// if ( ieta < 0 ) 
// {
//  etamax = -0.087-(fabs(ieta)-1)*0.087;
//  etamin = -(fabs(ieta)-1)*0.087;
// }
// theta1 = exp(-etamin); 
// theta2 = exp(-etamax);
// theta1 = 2.*atan(theta1); 
// theta2 = 2.*atan(theta2); 
// cout<<"theta1/theta2/pi "<<theta1<<"/"<<theta2<<"/"<<TMath::Pi()/2<<endl;
 double eta = p.eta();
 int ii;
 if ( eta > 0 )
 {
  for ( int i = 1; i < 16; i++ )
  {
   etamax = 0.087+(i-1)*0.087;
   etamin = (i-1)*0.087;
   if ( eta <= etamax && eta >= etamin ) ii = i;
  }
 }
 else
 {
  for ( int i = -1; i > -16; i-- )
  {
   etamin = -0.087-(fabs(i)-1)*0.087;
   etamax = -(fabs(i)-1)*0.087;
   if ( eta <= etamax && eta >= etamin ) ii = i;
  }
 }
 cout<<"eta/etamax/etamin/ieta "<<eta<<"/"<<etamax<<"/"<<etamin<<"/"<<ii<<endl;
 return ii;
}


//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(HBHEMuonCosmic2);




// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EVENT_DISPLAY_MAKER_H
#define EVENT_DISPLAY_MAKER_H

#include <fun4all/SubsysReco.h>

#include <KFParticle.h>

#include <trackbase/ActsGeometry.h>
#include <calotrigger/TriggerAnalyzer.h>
#include <calotrigger/TriggerRunInfo.h>
#include <ffarawobjects/Gl1Packet.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>
#include <kfparticle_sphenix/KFParticle_Container.h>
#include <trackbase/TrkrCluster.h>           // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <nlohmann/json.hpp>

#include <string>
#include <fstream>

using json = nlohmann::json;

class KFParticle_Container;
class PHCompositeNode;
class SvtxTrackMap;
class SvtxTrack;
class SvtxVertexMap;
class SvtxVertex;
class TriggerAnalyzer;

class event_display_maker : public SubsysReco
{
 public:

  event_display_maker(const std::string &name = "event_display_maker");

  ~event_display_maker() override;

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  void setKFParticleContainerName(const std::string &name) { m_kfparticleContainerName = name; }

  void setNumberDaughters(const unsigned int nDaught = 2) { m_number_of_daughters = nDaught; }

  void setEventDisplayPath(std::string path) { m_evt_display_path = path; }

  void setMassRange(float min = 0.48, float max = 0.51) { m_min_mass = min; m_max_mass = max; }

  void setMaxEvtDisplays(int max = 10) { m_max_displays = max; }

  void setMotherName(std::string name = "K_S0") { m_mother_name = name; }

  void setIntermediateNames(std::vector<std::string> name = {"Lambda0"}) { m_intermerdiate_names = name; }

 private:
  int counter = 0;
  int m_runNumber = 0;
  std::string m_evt_display_path = "./";
  std::string m_run_date = "2024-04-14";
  std::string m_mother_name = "K_S0";
  std::vector<std::string> m_intermerdiate_names;
  float m_min_mass = 0.4;
  float m_max_mass = 0.6;
  int m_max_displays = 100;

  int load_nodes(PHCompositeNode *topNode);
  std::string getParticleName(int ID);
  SvtxTrack *getTrack(unsigned int track_id, SvtxTrackMap *trackmap);
  void getJSONdata(std::vector<SvtxTrack*> tracks, std::vector<KFParticle*> particles, json &data, std::string jsonEntryName);
  std::string getDate();
  bool isInRange(float min, float value, float max);

  std::string m_kfparticleContainerName;
  std::string m_kfparticleContainerName_trackMap;
  std::string m_kfparticleContainerName_container;

  ActsGeometry *geometry {nullptr};
  SvtxTrack *m_track {nullptr};
  SvtxTrackMap *m_trackMap {nullptr};
  SvtxTrackMap *m_kfp_trackMap {nullptr};
  SvtxVertexMap *m_vertexMap {nullptr};
  SvtxVertex *m_vertex {nullptr};
  KFParticle_Container *m_kfp_container {nullptr};
  Gl1Packet *gl1packet {nullptr};
  TriggerAnalyzer *triggeranalyzer {nullptr};
  TriggerRunInfo *triggerruninfo {nullptr};
  TrkrClusterContainer *dst_clustermap {nullptr};

  KFParticle *mother {nullptr};
  std::vector<KFParticle*> intermediates;
  std::vector<KFParticle*> kfp_daughters;
  std::vector<SvtxTrack*> trigger_tracks;
  std::vector<SvtxTrack*> all_tracks;

  unsigned int m_number_of_daughters = 2; //Need something jsut now to ensure there are 2 candidates per BC

  std::ofstream json_output;

  std::string eventName = "EVENT";
  std::string metadataName = "META";
  std::string triggerTrackName = "TRIGGER";
  std::string allTrackName = "ALL";
  std::string trackName = "TRACKS";
  std::string hitsName = "HITS";
  std::string electronTriggerHitsName = "eHITS";
  std::string muonTriggerHitsName = "muHITS";
  std::string pionTriggerHitsName = "piHITS";
  std::string kaonTriggerHitsName = "kHITS";
  std::string protonTriggerHitsName = "pHITS";
  std::vector<std::string> hitsNameVector = {electronTriggerHitsName, muonTriggerHitsName, pionTriggerHitsName, kaonTriggerHitsName, protonTriggerHitsName, allTrackName, triggerTrackName};

  int electronColour = 16252672; //Yellow
  int muonColour = 16743424; //Orange
  int pionColour = 1834752; //Green
  int kaonColour = 16711680; //Red
  int protonColour = 12031; //Blue
  int triggerColour = 16711680; //Red
  int allColour = 16446450; //White
  int motherColour = 4251856; //Turquoise
  int intermediateColour = 16711935; //Purple

  std::map<int, int> pidMap{
  { 11, electronColour},
  {-11, electronColour},
  { 13, muonColour},
  {-13, muonColour},
  { 211, pionColour},
  {-211, pionColour},
  { 321, kaonColour},
  {-321, kaonColour},
  { 2212, protonColour},
  {-2212, protonColour}};

  std::map<int, std::string> pidToHitMap{
  { 11, electronTriggerHitsName},
  {-11, electronTriggerHitsName},
  { 13, muonTriggerHitsName},
  {-13, muonTriggerHitsName},
  { 211, pionTriggerHitsName},
  {-211, pionTriggerHitsName},
  { 321, kaonTriggerHitsName},
  {-321, kaonTriggerHitsName},
  { 2212, protonTriggerHitsName},
  {-2212, protonTriggerHitsName}};

  std::map< std::string, int> pidToColourMap{
  {electronTriggerHitsName, electronColour},
  {muonTriggerHitsName, muonColour},
  {pionTriggerHitsName, pionColour},
  {kaonTriggerHitsName, kaonColour},
  {protonTriggerHitsName, protonColour},
  {allTrackName, allColour},
  {triggerTrackName, triggerColour}};
};

#endif // EVENT_DISPLAY_MAKER_H

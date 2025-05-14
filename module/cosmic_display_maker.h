// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef COSMIC_DISPLAY_MAKER_H
#define COSMIC_DISPLAY_MAKER_H

#include <fun4all/SubsysReco.h>

#include <trackbase/ActsGeometry.h>
#include <calotrigger/TriggerAnalyzer.h>
#include <calotrigger/TriggerRunInfo.h>
#include <ffarawobjects/Gl1Packet.h>
#include <trackbase/TrkrCluster.h>           // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <nlohmann/json.hpp>

#include <string>
#include <fstream>

using json = nlohmann::json;

class PHCompositeNode;
class SvtxTrackMap;
class SvtxTrack;
class SvtxVertexMap;
class SvtxVertex;
class TriggerAnalyzer;

class cosmic_display_maker : public SubsysReco
{
 public:

  cosmic_display_maker(const std::string &name = "cosmic_display_maker");

  ~cosmic_display_maker() override;

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  void setEventDisplayPath(std::string path) { m_evt_display_path = path; }

  void setMaxEvtDisplays(int max = 10) { m_max_displays = max; }

  void setDecayTag(std::string tag) { m_decay_tag = tag; }

  void useTheseBCOs(std::vector<uint64_t> bcoList) { m_bco_list = bcoList; }

  void plotsApproved(bool approved = true) { m_approved = approved; }

 private:
  int counter = 0;
  int m_runNumber = 0;
  std::string m_evt_display_path = "./";
  std::string m_run_date = "2024-04-14";
  std::string m_decay_tag = "Cosmic ray";
  std::vector<uint64_t> m_bco_list;
  int m_max_displays = 100;
  bool m_approved = false;

  int load_nodes(PHCompositeNode *topNode);
  std::string getParticleName(int ID);
  SvtxTrack *getTrack(unsigned int track_id, SvtxTrackMap *trackmap);
  void getJSONdata(std::vector<SvtxTrack*> tracks, json &data);
  std::string getDate();
  bool isInRange(float min, float value, float max);

  ActsGeometry *geometry {nullptr};
  SvtxTrack *m_track {nullptr};
  SvtxTrackMap *m_trackMap {nullptr};
  Gl1Packet *gl1packet {nullptr};
  TrkrClusterContainer *dst_clustermap {nullptr};

  std::vector<SvtxTrack*> all_tracks;

  std::ofstream json_output;

  std::string eventName = "EVENT";
  std::string metadataName = "META";
  std::string triggerTrackName = "TRIGGER";
  std::string allTrackName = "ALL";
  std::string trackName = "TRACKS";
  std::string hitsName = "HITS";
  std::vector<std::string> hitsNameVector = {allTrackName, triggerTrackName};

  int yellow = 16252672;
  int orange = 16743424;
  int green = 1834752;
  int turquoise = 4251856;
  int blue = 12031;
  int red = 16711680;
  int pink = 16070058;
  int white = 16446450;

  int electronColour = yellow; 
  int muonColour = pink;
  int pionColour = green;
  int kaonColour = turquoise;
  int protonColour = orange;
  int motherColour = red; 
  int intermediateColour = blue;
  int triggerColour = red;
  int allColour = white;

  std::map<std::string, int> pidToColourMap{
  {allTrackName, allColour},
  {triggerTrackName, triggerColour}};
};

#endif // COSMIC_DISPLAY_MAKER_H

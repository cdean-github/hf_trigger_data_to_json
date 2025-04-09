// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef HF_TRIGGER_DATA_TO_JSON_H
#define HF_TRIGGER_DATA_TO_JSON_H

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

class hf_trigger_data_to_json : public SubsysReco
{
 public:

  hf_trigger_data_to_json(const std::string &name = "hf_trigger_data_to_json");

  ~hf_trigger_data_to_json() override;

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  void Print(const std::string &what = "ALL") const override;

  void setKFParticleContainerName(const std::string &name) { m_kfparticleContainerName = name; }

  void setNumberDaughters(const unsigned int nDaught = 2) { m_number_of_daughters = nDaught; }

  void setOutputFile(std::string name) { m_output_file = name; }

 private:
  int counter = 0;

  int load_nodes(PHCompositeNode *topNode);
  std::string getParticleName(int ID);
  void identify(const KFParticle &particle);
  SvtxTrack *getTrack(unsigned int track_id, SvtxTrackMap *trackmap);
  void getJSONdata(std::vector<SvtxTrack*> tracks, json &data, std::string jsonEntryName);

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

  std::vector<KFParticle*> kfp_daughters;
  std::vector<SvtxTrack*> trigger_tracks;
  std::vector<SvtxTrack*> all_tracks;

  unsigned int m_number_of_daughters = 2; //Need something jsut now to ensure there are 2 candidates per BC

  std::string m_output_file = "test.json";
  std::ofstream json_output;

  std::string metadataName = "Metadata";
  std::string mvtxClusterName = "MVTXClusters";
  std::string inttClusterName = "INTTClusters";
  std::string triggerTrackName = "TriggerTracks";
  std::string allTrackName = "AllTracks";
  std::string trackName = "Tracks";
};

#endif // HF_TRIGGER_DATA_TO_JSON_H

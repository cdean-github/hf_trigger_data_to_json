#include "hf_trigger_data_to_json.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <trackbase/InttDefs.h>              // for getLadderPhiId, getLad...
#include <trackbase/MvtxDefs.h>              // for getChipId, getStaveId
#include <trackbase/TpcDefs.h>               // for getSectorId, getSide
#include <trackbase/TrkrDefs.h>              // for getLayer, getTrkrId

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <TDatabasePDG.h>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

//____________________________________________________________________________..
hf_trigger_data_to_json::hf_trigger_data_to_json(const std::string &name)
 : SubsysReco(name)
 , m_kfparticleContainerName("reconstructedParticles")
{
}

//____________________________________________________________________________..
hf_trigger_data_to_json::~hf_trigger_data_to_json()
{
}

//____________________________________________________________________________..
int hf_trigger_data_to_json::Init(PHCompositeNode *topNode)
{
  load_nodes(topNode);

  triggeranalyzer = new TriggerAnalyzer();

  json_output.open(m_output_file);
  json_output << "{\n\"Events\" : [\n";

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int hf_trigger_data_to_json::process_event(PHCompositeNode *topNode)
{
  load_nodes(topNode);

  if (m_kfp_container->size() == 0) //No HF candidates
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  std::string trackableParticles[] = {"e-", "mu-", "pi+", "K+", "proton"};

  for (auto &vtx_iter : *m_vertexMap)
  {
    kfp_daughters.clear();
    trigger_tracks.clear();
    all_tracks.clear();
  
    m_vertex = vtx_iter.second;
  
    for (KFParticle_Container::Iter kfp_iter = m_kfp_container->begin(); kfp_iter != m_kfp_container->end(); ++kfp_iter)
    {
      KFParticle *myParticle = kfp_iter->second;
      std::string thisParticle = getParticleName(abs(myParticle->GetPDG())); //Only care about trackable particles

      if (std::find(std::begin(trackableParticles), std::end(trackableParticles), thisParticle) != std::end(trackableParticles))
      {
        m_track = getTrack(myParticle->Id(), m_trackMap);

        if (m_track->get_vertex_id() != m_vertex->get_id()) continue; //Trck isnt from vertex (need continue to account for mother)

        //identify(*myParticle);
        //m_track->identify();

        kfp_daughters.push_back(myParticle);
        trigger_tracks.push_back(m_track);
      }
    }

    if (kfp_daughters.size() != m_number_of_daughters) continue; //Check that we didnt have more than 1 candidate in this BC

    //Get all tracks associated to this vertex
    for (SvtxVertex::TrackIter all_track_iter =  m_vertex->begin_tracks(); all_track_iter != m_vertex->end_tracks(); ++all_track_iter)
    {
      m_track = m_trackMap->get(*all_track_iter);
      all_tracks.push_back(m_track);
      //m_track->identify();
    }

    triggeranalyzer->decodeTriggers(topNode);

    //Now fill the json file
    json data;

    uint64_t m_bco = gl1packet->lValue(0, "BCO") + trigger_tracks[0]->get_crossing();
    data[metadataName]["EventID"] = m_bco;
    data[metadataName]["CollisionVertex"] = {m_vertex->get_x(), m_vertex->get_y(), m_vertex->get_z()};

    getJSONdata(trigger_tracks, data, triggerTrackName);
    getJSONdata(all_tracks, data, allTrackName);

    if (counter != 0) json_output << ",\n";

    json_output << data.dump(2);

    ++counter;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int hf_trigger_data_to_json::End(PHCompositeNode *topNode)
{

  json_output << "\n]\n}";
  json_output.close();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void hf_trigger_data_to_json::Print(const std::string &what) const
{
}

int hf_trigger_data_to_json::load_nodes(PHCompositeNode *topNode)
{
  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!m_vertexMap)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error : missing SvtxVertexMap" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!m_trackMap)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error : missing SvtxTrackMap" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_kfparticleContainerName_container = m_kfparticleContainerName + "_KFParticle_Container";

  m_kfp_container = findNode::getClass<KFParticle_Container>(topNode, m_kfparticleContainerName_container);
  if (!m_kfp_container)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error : missing " << m_kfparticleContainerName_container << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  gl1packet = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
  if (!gl1packet)
  {
    gl1packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");
    if (!gl1packet)
    {
      std::cout << __PRETTY_FUNCTION__ << " Fatal Error : missing Gl1Packet" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }

  triggerruninfo = findNode::getClass<TriggerRunInfo>(topNode, "TriggerRunInfo");
  if (!triggerruninfo)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error : missing TriggerRunInfo" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  dst_clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!dst_clustermap)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error : missing TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!geometry)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error : missing ActsGeometry" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

std::string hf_trigger_data_to_json::getParticleName(int ID)
{
  return TDatabasePDG::Instance()->GetParticle(ID)->GetName();
}

void hf_trigger_data_to_json::identify(const KFParticle &particle)
{
  std::cout << "Track ID: " << particle.Id() << std::endl;
  std::cout << "PDG ID: " << particle.GetPDG() << ", charge: " << (int) particle.GetQ() << ", mass: " << particle.GetMass() << " GeV" << std::endl;
  std::cout << "(px,py,pz) = (" << particle.GetPx() << " +/- " << std::sqrt(particle.GetCovariance(3, 3)) << ", ";
  std::cout << particle.GetPy() << " +/- " << std::sqrt(particle.GetCovariance(4, 4)) << ", ";
  std::cout << particle.GetPz() << " +/- " << std::sqrt(particle.GetCovariance(5, 5)) << ") GeV" << std::endl;
  std::cout << "(x,y,z) = (" << particle.GetX() << " +/- " << std::sqrt(particle.GetCovariance(0, 0)) << ", ";
  std::cout << particle.GetY() << " +/- " << std::sqrt(particle.GetCovariance(1, 1)) << ", ";
  std::cout << particle.GetZ() << " +/- " << std::sqrt(particle.GetCovariance(2, 2)) << ") cm\n"
            << std::endl;
}

SvtxTrack *hf_trigger_data_to_json::getTrack(unsigned int track_id, SvtxTrackMap *trackmap)
{
  SvtxTrack *matched_track = nullptr;

  for (auto &iter : *trackmap)
  {
    if (iter.first == track_id)
    {
      matched_track = iter.second;
    }
  }

  return matched_track;
}

void hf_trigger_data_to_json::getJSONdata(std::vector<SvtxTrack*> tracks, json &jsonData, std::string jsonEntryName)
{
  json tracksJson, mvtxClusters, inttClusters;

  for (auto track : tracks)
  {
    std::vector<uint64_t> mvtxClusterKeys, inttClusterKeys;

    unsigned int trackID = track->get_id();

    std::vector<float> trackPosition = {track->get_x(), track->get_y(), track->get_z()};
    std::vector<float> trackMomentum = {track->get_px(), track->get_py(), track->get_pz()};

    tracksJson["TrackSequenceInEvent"] = trackID;
    tracksJson["TrackMomentum"] = trackMomentum;
    tracksJson["TrackPosition"] = trackPosition;

    TrackSeed *silseed = track->get_silicon_seed();
    if (!silseed) continue; //Track has no silicon seeds

    for (auto cluster_iter = silseed->begin_cluster_keys();
    cluster_iter != silseed->end_cluster_keys(); ++cluster_iter)
    {
      uint64_t clusterKey = *cluster_iter;

      TrkrCluster *cluster = dst_clustermap->findCluster(clusterKey);
      auto global = geometry->getGlobalPosition(clusterKey, cluster);

      std::vector<float> clusterPosition = {(float) global.x(), (float) global.y(), (float) global.z()};

      uint8_t id = TrkrDefs::getTrkrId(clusterKey);
      uint8_t layer = TrkrDefs::getLayer(clusterKey);

      uint16_t stave = 0;
      uint16_t chip = 0;
      uint16_t row = 0;
      uint16_t col = 0;
      uint16_t zID = 0;
      uint16_t phiID = 0; 
     
      json mvtxCluster;
      json inttCluster;

      switch (id)
      {
        case TrkrDefs::mvtxId:
          stave = MvtxDefs::getStaveId(clusterKey);
          chip = MvtxDefs::getChipId(clusterKey);
          row = MvtxDefs::getRow(clusterKey);
          col = MvtxDefs::getCol(clusterKey);

          mvtxCluster["ID"]["ClusKey"] = clusterKey;
          mvtxCluster["ID"]["Layer"] = layer;
          mvtxCluster["ID"]["Stave"] = stave;
          mvtxCluster["ID"]["Chip"] = chip;
          mvtxCluster["ID"]["Row"] = row;
          mvtxCluster["ID"]["Col"] = col;
          mvtxCluster["Coordinate"] = clusterPosition;
          mvtxClusters += mvtxCluster;

          mvtxClusterKeys.push_back(clusterKey);

          break;

        case TrkrDefs::inttId:
          zID = InttDefs::getLadderZId(clusterKey);
          phiID = InttDefs::getLadderPhiId(clusterKey);
          row = InttDefs::getRow(clusterKey);
          col = InttDefs::getCol(clusterKey);

          inttCluster["ID"]["ClusKey"] = clusterKey;
          inttCluster["ID"]["Layer"] = layer;
          inttCluster["ID"]["LadderZId"] = zID;
          inttCluster["ID"]["LadderPhiId"] = phiID;
          inttCluster["ID"]["Row"] = row;
          inttCluster["ID"]["Col"] = col;
          inttCluster["Coordinate"] = clusterPosition;
          inttClusters += inttCluster;

          inttClusterKeys.push_back(clusterKey);

          break;

        default:
          break; 
      } //End of switch
    } //End of state loop

    tracksJson["MVTXclusterKeys"] = mvtxClusterKeys;
    tracksJson["INTTclusterKeys"] = inttClusterKeys;
    jsonData[jsonEntryName][trackName] += tracksJson;

  } //End of track loop

  jsonData[jsonEntryName][mvtxClusterName] = mvtxClusters;
  jsonData[jsonEntryName][inttClusterName] = inttClusters;
}

#include "cosmic_display_maker.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <trackbase/InttDefs.h>              // for getLadderPhiId, getLad...
#include <trackbase/MvtxDefs.h>              // for getChipId, getStaveId
#include <trackbase/TpcDefs.h>               // for getSectorId, getSide
#include <trackbase/TrkrDefs.h>              // for getLayer, getTrkrId

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/recoConsts.h>

#include <trackbase/TrkrDefs.h>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

//____________________________________________________________________________..
cosmic_display_maker::cosmic_display_maker(const std::string &name)
 : SubsysReco(name)
{
}

//____________________________________________________________________________..
cosmic_display_maker::~cosmic_display_maker()
{
}

//____________________________________________________________________________..
int cosmic_display_maker::Init(PHCompositeNode *topNode)
{
  load_nodes(topNode);

  recoConsts *rc = recoConsts::instance();
  m_runNumber = rc->get_IntFlag("RUNNUMBER");

  m_run_date = getDate();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int cosmic_display_maker::process_event(PHCompositeNode *topNode)
{
  if (counter >= m_max_displays) return Fun4AllReturnCodes::EVENT_OK; //Made the max number of displays

  load_nodes(topNode);

  all_tracks.clear();

  if (m_trackMap->size() == 0) return Fun4AllReturnCodes::EVENT_OK;;

  for (auto &iter : *m_trackMap)
  {
    m_track = iter.second;
    all_tracks.push_back(m_track);
  }

  //Now fill the json file
  json data;

  uint64_t m_bco = gl1packet->lValue(0, "BCO") + all_tracks[0]->get_crossing();

  //Check if we have a list of BCOs we want a plot for, and this is one of them
  if (m_bco_list.size() != 0 && std::find(std::begin(m_bco_list), std::end(m_bco_list), m_bco) == std::end(m_bco_list)) return Fun4AllReturnCodes::EVENT_OK;

  std::string date_run_stamp = m_run_date + ", Run " + std::to_string(m_runNumber);
  std::string bco_stamp = "BCO: " + std::to_string(m_bco);

  data[eventName]["evtid"] = 1;
  data[eventName]["runid"] = m_runNumber;
  data[eventName]["type"] = "Collision";
  data[eventName]["s_nn"] = 0;
  data[eventName]["B"] = 1.38;

  if (m_approved)
  {
    data[eventName]["runstats"] = {"sPHENIX Tracking", date_run_stamp, m_decay_tag};
  }
  else
  {
    data[eventName]["runstats"] = {"sPHENIX Internal", date_run_stamp, bco_stamp,  m_decay_tag};
  }

  for (auto& hitsMetaSetup : hitsNameVector)
  {      
    data[metadataName][hitsName][hitsMetaSetup]["type"] = "3D";
    data[metadataName][hitsName][hitsMetaSetup]["options"]["size"] = 5;
    data[metadataName][hitsName][hitsMetaSetup]["options"]["color"] = pidToColourMap[hitsMetaSetup];
  }

  data[metadataName][trackName][triggerTrackName]["width"] = 0.3;

  data[trackName]["B"] = 0.000014;

  getJSONdata(all_tracks, data);

  if (!hasAll) return Fun4AllReturnCodes::EVENT_OK; 

  std::string m_output_file = m_evt_display_path + "/EvtDisplay_cosmic_" + std::to_string(m_runNumber) + "_" + std::to_string(m_bco) + ".json";

  json_output.open(m_output_file);
  json_output << data.dump(2);
  json_output.close();

  ++counter;
  //}

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int cosmic_display_maker::End(PHCompositeNode *topNode)
{
  std::cout << "Made " << counter << " event displays" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int cosmic_display_maker::load_nodes(PHCompositeNode *topNode)
{
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!m_trackMap)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error : missing SvtxTrackMap" << std::endl;
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

void cosmic_display_maker::getJSONdata(std::vector<SvtxTrack*> tracks, json &jsonData)
{
  json tracksJson, clustersJson, allClustersJson;

  for (auto& det : {TrkrDefs::TrkrId::mvtxId, TrkrDefs::TrkrId::inttId,
                    TrkrDefs::TrkrId::tpcId, TrkrDefs::TrkrId::micromegasId})
  {
    for (const auto& hitsetkey : dst_clustermap->getHitSetKeys(det))
    {
      auto range = dst_clustermap->getClusters(hitsetkey);
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        TrkrDefs::cluskey clusterKey = iter->first;
        TrkrCluster* cluster = dst_clustermap->findCluster(clusterKey);
        auto global = geometry->getGlobalPosition(clusterKey, cluster);
        allClustersJson["x"] = (float) global.x();
        allClustersJson["y"] = (float) global.y();
        allClustersJson["z"] = (float) global.z();
        allClustersJson["e"] = 0;
    
        jsonData[hitsName][allTrackName] += allClustersJson;
      }
    }
  }

  for (auto track : tracks)
  {
    float length = 100;

    TrackSeed *silseed = track->get_silicon_seed();
    if (silseed)
    {
      for (auto cluster_iter = silseed->begin_cluster_keys(); cluster_iter != silseed->end_cluster_keys(); ++cluster_iter)
      {
        uint64_t clusterKey = *cluster_iter;
    
        TrkrCluster *cluster = dst_clustermap->findCluster(clusterKey);
        auto global = geometry->getGlobalPosition(clusterKey, cluster);

        clustersJson["x"] = (float) global.x();
        clustersJson["y"] = (float) global.y();
        clustersJson["z"] = (float) global.z();
        clustersJson["e"] = 0;

        jsonData[hitsName][triggerTrackName] += clustersJson;
  
      }
    }

    TrackSeed *tpcseed = track->get_tpc_seed();
    if (tpcseed)
    {
      for (auto cluster_iter = tpcseed->begin_cluster_keys(); cluster_iter != tpcseed->end_cluster_keys(); ++cluster_iter)
      {
        uint64_t clusterKey = *cluster_iter;
    
        TrkrCluster *cluster = dst_clustermap->findCluster(clusterKey);
        auto global = geometry->getGlobalPosition(clusterKey, cluster);

        clustersJson["x"] = (float) global.x();
        clustersJson["y"] = (float) global.y();
        clustersJson["z"] = (float) global.z();
        clustersJson["e"] = 0;

        jsonData[hitsName][triggerTrackName] += clustersJson;
      }
    }
    
    std::vector<float> trackPosition = {track->get_x(), track->get_y(), track->get_z()};
    std::vector<float> trackMomentum = {track->get_px(), track->get_py(), track->get_pz()};

    tracksJson["pxyz"] = trackMomentum;
    tracksJson["xyz"] = trackPosition;
    tracksJson["trk_color"] = triggerColour;
    tracksJson["nh"] = 60;
    tracksJson["l"] = length;
    tracksJson["q"] = track->get_charge();

    jsonData[trackName][triggerTrackName] += tracksJson;    

  } //End of track loop
}

//https://stackoverflow.com/questions/997946/how-to-get-current-time-and-date-in-c
std::string cosmic_display_maker::getDate()
{
    std::time_t t = std::time(0);   // get time now
    std::tm* now = std::localtime(&t);

    std::stringstream date;
    date << (now->tm_year + 1900) << '-'
         << (now->tm_mon + 1) << '-'
         <<  now->tm_mday;
    return date.str();
}

bool cosmic_display_maker::isInRange(float min, float value, float max)
{
  return min <= value && value <= max;
}

#include "betabeta-data.cc"

#include <bayeux/dpp/chain_module.h>
#include <bayeux/datatools/service_manager.h>

#include <falaise/snemo/datamodels/event_header.h>
#include <falaise/snemo/datamodels/geomid_utils.h>
#include <falaise/snemo/datamodels/calibrated_calorimeter_hit.h>
#include <falaise/snemo/datamodels/calibrated_data.h>
#include <falaise/snemo/datamodels/line_trajectory_pattern.h>
#include <falaise/snemo/datamodels/particle_track_data.h>
#include <falaise/snemo/datamodels/precalibrated_data.h>
#include <falaise/snemo/datamodels/precalibrated_tracker_hit.h>
#include <falaise/snemo/datamodels/tracker_clustering_data.h>
#include <falaise/snemo/datamodels/tracker_trajectory_data.h>
// #include <falaise/snemo/services/geometry.h>
// #include <falaise/snemo/services/service_handle.h>
#include <falaise/property_set.h>

#include <cstdio>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"

////////////////////////////////////////////////////////////////

class betabeta_selection_module : public dpp::chain_module
{
public:
  betabeta_selection_module();
  ~betabeta_selection_module();

  void initialize(const datatools::properties &,
		  datatools::service_manager &,
		  dpp::module_handle_dict_type &);
  
  dpp::chain_module::process_status process(datatools::things &);

  void finalize();

private:
  std::string _ptd_label_;
  std::string _output_filename_;

  betabeta_data _betabeta_data_;

  TFile *_output_file_;
  TTree *_output_tree_;

  DPP_MODULE_REGISTRATION_INTERFACE(betabeta_selection_module);
};

DPP_MODULE_REGISTRATION_IMPLEMENT(betabeta_selection_module, "betabeta_selection_module");

////////////////////////////////////////////////////////////////

betabeta_selection_module::betabeta_selection_module() : dpp::chain_module()
{
  
}

betabeta_selection_module::~betabeta_selection_module()
{
  this->finalize();
}

void betabeta_selection_module::initialize(const datatools::properties & config,
			   datatools::service_manager &, dpp::module_handle_dict_type &)
{
  const falaise::property_set fps {config};

  _ptd_label_ = fps.get<std::string>("PTD_label", "PTD");
  _output_filename_ = fps.get<std::string>("output_filename", "output.root");

  // prepare ROOT output file/tree
  _output_file_ = new TFile(_output_filename_.c_str(), "RECREATE");
  _output_tree_ = new TTree("betabeta_tree", "");

  _output_tree_->Branch("betabeta", &_betabeta_data_);

  gROOT->cd();

  this->_set_initialized(true);
}

dpp::chain_module::process_status betabeta_selection_module::process(datatools::things &event)
{
  if (!event.has("PTD"))
    {
      std::cout << "*** no PTD bank ***" << std::endl;
      return dpp::base_module::PROCESS_ERROR;
    }

  const snemo::datamodel::event_header & EH = event.get<snemo::datamodel::event_header>("EH");
  const int eh_run = EH.get_id().get_run_number();
  const int eh_event = EH.get_id().get_event_number();

  // const snemo::datamodel::calibrated_data & CD = event.get<snemo::datamodel::calibrated_data>("CD");
  // const snemo::datamodel::tracker_clustering_solution & tcd_solution = TCD.get_default();
  // const snemo::datamodel::tracker_trajectory_data & TTD = event.get<snemo::datamodel::tracker_trajectory_data>("TTD");
  // const snemo::datamodel::tracker_trajectory_solution & ttd_solution = TTD.get_default_solution();



  ////////////////////////////////////////////////////
  // compute average anode time per tracker cluster //
  ////////////////////////////////////////////////////

  const snemo::datamodel::precalibrated_data & pCD = event.get<snemo::datamodel::precalibrated_data>("pCD");
  const snemo::datamodel::tracker_clustering_data & TCD = event.get<snemo::datamodel::tracker_clustering_data>("TCD");
  const snemo::datamodel::tracker_clustering_solution & tcd_solution = TCD.get_default();

  std::vector<double> cluster_mean_anodic_time;
  cluster_mean_anodic_time.reserve(tcd_solution.get_clusters().size());

  for (const datatools::handle<snemo::datamodel::tracker_cluster> & cluster : tcd_solution.get_clusters()) {

    double anodic_time_sum = 0;

    for (const datatools::handle<snemo::datamodel::calibrated_tracker_hit>  & tracker_hit : cluster->hits()) {
      const int pcd_index = tracker_hit->get_auxiliaries().fetch_integer("pCD.parent");
      anodic_time_sum += pCD.tracker_hits()[pcd_index]->get_anodic_time();
    }

    cluster_mean_anodic_time.push_back(anodic_time_sum / cluster->hits().size());
  }

  ///////////////////////////////////////////////////////////////////////
  // loop over PTD bank and classify particle as alpha, beta and gamma //
  ///////////////////////////////////////////////////////////////////////

  const snemo::datamodel::particle_track_data & PTD = event.get<snemo::datamodel::particle_track_data>("PTD");

  std::vector<alpha_data> alphas;
  std::vector<beta_data>   betas;
  std::vector<gamma_data> gammas;

  std::vector<betabeta_data> betabetas;

  for (const datatools::handle<snemo::datamodel::particle_track> & particle : PTD.particles())
    {
      beta_data new_beta;
      memset(&new_beta, 0, sizeof(new_beta));

      bool has_source_vertex = false;
      bool has_calo_vertex = false;

      for (const datatools::handle<snemo::datamodel::vertex> & vertex : particle->get_vertices())
	{
	  if (vertex->is_on_main_calorimeter())
	    {
	      const geomtools::vector_3d & vertex_position = vertex->get_spot().get_position();

	      for (int xyz=0; xyz<3; xyz++)
		new_beta.calo_vtx[xyz] = vertex_position[xyz] / CLHEP::m;

	      has_calo_vertex = true;
	    }

	  else if (vertex->is_on_reference_source_plane())
	    {
	      const geomtools::vector_3d & vertex_position = vertex->get_spot().get_position();

	      for (int xyz=0; xyz<3; xyz++)
		new_beta.source_vtx[xyz] = vertex_position[xyz] / CLHEP::m;

	      has_source_vertex = true;
	    }

	} // for (vertex)

      if (particle->get_trajectory().get_cluster().is_delayed() || (!has_calo_vertex) || (!has_source_vertex)) {

	// consider this particule as an alpha;

	alpha_data new_alpha;
	// new_alpha.length;
	// new_alpha.vtx0[3];
	// new_alpha.vtx1[3];
	// new_alpha.nb_cells;

	new_alpha.mean_anode = cluster_mean_anodic_time[particle->get_trajectory().get_cluster().get_cluster_id()];

	alphas.push_back(new_alpha);

	continue;
      }
      
      if (particle->get_associated_calorimeter_hits().size() == 0)
	continue;

      if (particle->get_associated_calorimeter_hits().size() > 1)
	printf("[%d_%d] *** particle track with > 1 calo\n", eh_run, eh_event);

      const auto & calo_hit = particle->get_associated_calorimeter_hits().front();
      new_beta.energy = calo_hit->get_energy() / CLHEP::MeV;
      new_beta.time = calo_hit->get_time() / CLHEP::ns;

      for (int xyz=0; xyz<3; xyz++)
	new_beta.length += std::pow(new_beta.calo_vtx[xyz]-new_beta.source_vtx[xyz],2);
      new_beta.length = std::sqrt(new_beta.length);

      new_beta.om_num = snemo::datamodel::om_num(calo_hit->get_geom_id());

      new_beta.nb_cells = particle->get_trajectory().get_cluster().size();
      new_beta.cluster_id = particle->get_trajectory().get_cluster().get_cluster_id();
      new_beta.particle_id = particle->get_track_id();

      // for (const auto & tracker_hit : particle->get_trajectory().get_cluster().hits())
      // 	new_beta.mean_anode += tracker_hit->get_anode_time() / CLHEP::ns; // utiliser pCD anodic time ...
      // new_beta.mean_anode /= particle->get_trajectory().get_cluster().hits().size();
      new_beta.mean_anode = cluster_mean_anodic_time[new_beta.cluster_id];

      betas.push_back(new_beta);
      // compute PTD length

      // calo time+energy

    } // for (particle)

  // printf("[%d_%d] %zd beta(s)\n", eh_run, eh_event, betas.size());


  if (betas.size() < 2)
    return dpp::base_module::PROCESS_SUCCESS;

  for (size_t i=0; i<betas.size(); i++) {

    const beta_data & beta_i = betas[i];

    for (size_t j=i+1; j<betas.size(); j++) {

      // printf("[%d_%d] (%zd,%zd)\n", eh_run, eh_event, i, j);

      const beta_data & beta_j = betas[j];

      if (beta_i.cluster_id == beta_j.cluster_id) {
	// printf("[%d_%d] skipping same cluster beta pair\n", eh_run, eh_event);
	continue;
      }

      if (beta_i.om_num == beta_j.om_num) {
	// printf("[%d_%d] skipping beta pair with same OM\n", eh_run, eh_event);
	continue;
      }

      if (std::abs(beta_i.time - beta_j.time) > 1000) {
	// printf("[%d_%d] skipping time separated beta pair\n", eh_run, eh_event);
	continue;
      }

      if (beta_i.om_num < beta_j.om_num) {
	std::memcpy(&_betabeta_data_.beta[0], &beta_i, sizeof(beta_i));
	std::memcpy(&_betabeta_data_.beta[1], &beta_j, sizeof(beta_j));
      } else {
	std::memcpy(&_betabeta_data_.beta[0], &beta_j, sizeof(beta_j));
	std::memcpy(&_betabeta_data_.beta[1], &beta_i, sizeof(beta_i));
      }

      for (int xyz=0; xyz<3; xyz++) {
	_betabeta_data_.source_vtx[xyz] = beta_i.source_vtx[xyz];
	_betabeta_data_.source_vtx[xyz] += beta_j.source_vtx[xyz];
	_betabeta_data_.source_vtx[xyz] /= 2;
      }	

      _betabeta_data_.deltay = _betabeta_data_.beta[1].source_vtx[1] - _betabeta_data_.beta[0].source_vtx[1];
      _betabeta_data_.deltaz = _betabeta_data_.beta[1].source_vtx[2] - _betabeta_data_.beta[0].source_vtx[2];

      const float deltar2 = _betabeta_data_.deltay*_betabeta_data_.deltay + _betabeta_data_.deltaz*_betabeta_data_.deltaz;

      if (deltar2 > 0.3*0.3) {
	// printf("[%d_%d] skip betabeta candidates with large deltar\n", eh_run, eh_event);
	continue;
      }

      _betabeta_data_.run = eh_run;
      _betabeta_data_.run = eh_event;

      betabetas.push_back(_betabeta_data_);

    } // for (j)

  } // for (i)

  ////////////////////////////////

  int nb_betabeta = 0;

  if (betabetas.size() == 1)
    _output_tree_->Fill();

  else if (betabetas.size() > 1) {

    // if (betabetas.size() > 3)
    //   printf("[%d_%d] %zd betabeta candidates\n", eh_run, eh_event, betabetas.size());

    std::vector<bool> skip_betabeta;
    skip_betabeta.reserve(betabetas.size());
    for (size_t bb=0; bb<betabetas.size(); bb++)
      skip_betabeta.push_back(false);

    for (size_t bb1=0; bb1<betabetas.size(); bb1++) {

      if (skip_betabeta[bb1])
	continue;

      betabeta_data *best_betabeta = & betabetas[bb1];

      // look in other bb candidate for possible duplicated cases
      // (when multiple trajectory solution exists on same cluster)

      for (size_t bb2=1+bb1; bb2<betabetas.size(); bb2++) {

	if (skip_betabeta[bb2])
	  continue;

	betabeta_data *new_betabeta = & betabetas[bb2];

	if ((best_betabeta->beta[0].cluster_id == new_betabeta->beta[0].cluster_id) &&
	    (best_betabeta->beta[1].cluster_id == new_betabeta->beta[1].cluster_id)) {

	  if (TMath::Abs(new_betabeta->deltay) < TMath::Abs(best_betabeta->deltay)) {

	    // printf("[%d_%d] better betabeta selected (%d,%d) -> (%d,%d)\n", eh_run, eh_event,
	    // 	   best_betabeta->beta[0].particle_id, best_betabeta->beta[1].particle_id,
	    // 	   new_betabeta->beta[0].particle_id, new_betabeta->beta[1].particle_id);

	    best_betabeta = new_betabeta;

	  } else {

	    // printf("[%d_%d] keep betabeta selected (%d,%d) <- (%d,%d)\n", eh_run, eh_event,
	    // 	   best_betabeta->beta[0].particle_id, best_betabeta->beta[1].particle_id,
	    // 	   new_betabeta->beta[0].particle_id, new_betabeta->beta[1].particle_id);

	  }

	  skip_betabeta[bb2] = true;

	}

	else if ((best_betabeta->beta[0].cluster_id == new_betabeta->beta[1].cluster_id) &&
		 (best_betabeta->beta[1].cluster_id == new_betabeta->beta[0].cluster_id)) {

	  printf("[%d_%d] !!!!!!!! swapped double betabeta candidates !!!!!!\n", eh_run, eh_event);

	}

      } // for (bb2)

      // iterate over all other clusters for possible
      // extra tracks with space/time correlation

      const double betabeta_mean_anode = (best_betabeta->beta[0].mean_anode + best_betabeta->beta[1].mean_anode) / 2.0;

      for (size_t cluster_id=0; cluster_id<tcd_solution.get_clusters().size(); cluster_id++) {

	if (cluster_id == best_betabeta->beta[0].cluster_id) continue;
	if (cluster_id == best_betabeta->beta[1].cluster_id) continue;

	printf("[%d_%d] extra cluster with deltat = %.3f us and deltar = xxx\n", eh_run, eh_event,
	       (cluster_mean_anodic_time[cluster_id] - betabeta_mean_anode)/CLHEP::microsecond);

      }

      std::memcpy(&_betabeta_data_, best_betabeta, sizeof(betabeta_data));
      _output_tree_->Fill();

      nb_betabeta++;

    } // for (bb1)

    if (nb_betabeta > 1)
      printf("[%d_%d] %d betabeta candidates !\n", eh_run, eh_event, nb_betabeta);
  }

  return dpp::base_module::PROCESS_SUCCESS;
}

void betabeta_selection_module::finalize()
{
  _output_file_->cd();
  _output_tree_->Write("", TObject::kOverwrite);
  _output_file_->Close();

  this->_set_initialized(false);
}

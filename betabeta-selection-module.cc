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
#include <limits>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
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

  float _deltat_extra_calo_;
  float _deltat_extra_cluster_;
  float _deltar_delayed_cluster_;

  TFile *_output_file_;
  TTree *_output_tree_;

  TH2F *_output_deltart_bb_cluster_;
  TH1F *_output_deltat_bb_calo_;

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

  _deltat_extra_calo_ = fps.get<double>("deltat_extra_calo", 25.0);
  _deltat_extra_cluster_ = fps.get<double>("deltat_extra_cluster", 5.0);
  _deltar_delayed_cluster_ = fps.get<double>("deltar_delayed_cluster", 0.5);

  // prepare ROOT output file/tree
  _output_file_ = new TFile(_output_filename_.c_str(), "RECREATE");
  _output_tree_ = new TTree("betabeta_tree", "");
  _output_tree_->Branch("betabeta", &_betabeta_data_);

  _output_deltart_bb_cluster_ = new TH2F( "deltart_bb_cluster", ";#DeltaT(cluster - betabeta) (us);#DeltaR(cluster - betabeta) (m)", 500, -125, 125, 500, 0, 5);
  _output_deltat_bb_calo_ = new TH1F( "deltat_bb_calo", ";#DeltaT(calo - betabeta) (ns);", 1000, -1000, 1000);

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

  // const snemo::datamodel::tracker_clustering_solution & tcd_solution = TCD.get_default();
  // const snemo::datamodel::tracker_trajectory_data & TTD = event.get<snemo::datamodel::tracker_trajectory_data>("TTD");
  // const snemo::datamodel::tracker_trajectory_solution & ttd_solution = TTD.get_default_solution();

  const snemo::datamodel::precalibrated_data & pCD = event.get<snemo::datamodel::precalibrated_data>("pCD");
  const snemo::datamodel::calibrated_data & CD = event.get<snemo::datamodel::calibrated_data>("CD");
  const snemo::datamodel::tracker_clustering_data & TCD = event.get<snemo::datamodel::tracker_clustering_data>("TCD");
  const snemo::datamodel::tracker_clustering_solution & tcd_solution = TCD.get_default();

  /////////////////////////////////////////////
  // pre-compute some parameters of clusters //
  /////////////////////////////////////////////

  // side
  std::vector<uint16_t> cluster_side;
  cluster_side.reserve(tcd_solution.get_clusters().size());

  // number of cells with reconstructed Z
  std::vector<uint16_t> cluster_nb_cells_withz;
  cluster_nb_cells_withz.reserve(tcd_solution.get_clusters().size());

  // cluster mean anodic time (use as "cluster time")
  std::vector<double> cluster_mean_anodic_time;
  cluster_mean_anodic_time.reserve(tcd_solution.get_clusters().size());

  for (const datatools::handle<snemo::datamodel::tracker_cluster> & cluster : tcd_solution.get_clusters()) {

    uint16_t nb_cells_withz = 0;
    double anodic_time_sum = 0;

    cluster_side.push_back(cluster->hits().front()->get_side());

    for (const datatools::handle<snemo::datamodel::calibrated_tracker_hit> & tracker_hit : cluster->hits()) {

      if (datatools::is_valid(tracker_hit->get_z()))
	nb_cells_withz++;

      const int pcd_index = tracker_hit->get_auxiliaries().fetch_integer("pCD.parent");
      anodic_time_sum += pCD.tracker_hits()[pcd_index]->get_anodic_time();

    } // for (tracker_hit)

    cluster_nb_cells_withz.push_back(nb_cells_withz);
    cluster_mean_anodic_time.push_back(anodic_time_sum / cluster->hits().size());

  } // for (cluster)

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
      bool has_calo_vertex = false;
      bool has_source_vertex = false;

      float calo_vtx[3];
      float source_vtx[3];

      for (const datatools::handle<snemo::datamodel::vertex> & vertex : particle->get_vertices())
	{
	  if (vertex->is_on_main_calorimeter())
	    {
	      const geomtools::vector_3d & vertex_position = vertex->get_spot().get_position();

	      for (int xyz=0; xyz<3; xyz++)
		calo_vtx[xyz] = vertex_position[xyz] / CLHEP::m;

	      has_calo_vertex = true;
	    }

	  else if (vertex->is_on_reference_source_plane())
	    {
	      const geomtools::vector_3d & vertex_position = vertex->get_spot().get_position();

	      for (int xyz=0; xyz<3; xyz++)
		source_vtx[xyz] = vertex_position[xyz] / CLHEP::m;

	      has_source_vertex = true;
	    }

	} // for (vertex)

      if (particle->get_trajectory().get_cluster().is_delayed() || (!has_calo_vertex) || (!has_source_vertex) 
	  || (particle->get_associated_calorimeter_hits().size() == 0)) {

	// if : - particle is delayed
	//      - particle does not have calo vertex and source vertex
	//      - particle is not associated to a calorimeter hit
	// => consider this particle as an alpha !

	alpha_data new_alpha;
	// new_alpha.length;
	// new_alpha.vtx0[3];
	// new_alpha.vtx1[3];
	// new_alpha.nb_cells;

	new_alpha.mean_anode = cluster_mean_anodic_time[particle->get_trajectory().get_cluster().get_cluster_id()];

	alphas.push_back(new_alpha);

	continue;
      }
      
      // if (particle->get_associated_calorimeter_hits().size() == 0)
      // 	continue;

      if (particle->get_associated_calorimeter_hits().size() > 1)
	printf("[%d_%d] *** particle track with > 1 calo\n", eh_run, eh_event);

      // this particle is considered beta => prepare data storage

      betas.push_back(beta_data());
      beta_data & new_beta = betas.back();
      // memset(&new_beta, 0, sizeof(new_beta));

      const auto & calo_hit = particle->get_associated_calorimeter_hits().front();
      new_beta.energy = calo_hit->get_energy() / CLHEP::MeV;
      new_beta.time = calo_hit->get_time() / CLHEP::ns;

      std::memcpy(new_beta.calo_vtx, calo_vtx, sizeof(calo_vtx));
      std::memcpy(new_beta.source_vtx, source_vtx, sizeof(source_vtx));

      for (int xyz=0; xyz<3; xyz++)
	new_beta.length += std::pow(new_beta.calo_vtx[xyz]-new_beta.source_vtx[xyz],2);
      new_beta.length = std::sqrt(new_beta.length);

      new_beta.cluster_id = particle->get_trajectory().get_cluster().get_cluster_id();
      new_beta.particle_id = particle->get_track_id();

      new_beta.flag = cluster_side[new_beta.cluster_id];
      new_beta.om_num = snemo::datamodel::om_num(calo_hit->get_geom_id());

      new_beta.nb_cells = particle->get_trajectory().get_cluster().size();
      new_beta.nb_cells_withz = cluster_nb_cells_withz[new_beta.cluster_id];

      new_beta.nb_cells = particle->get_trajectory().get_cluster().size();

      // for (const auto & tracker_hit : particle->get_trajectory().get_cluster().hits())
      // 	new_beta.mean_anode += tracker_hit->get_anode_time() / CLHEP::ns; // utiliser pCD anodic time ...
      // new_beta.mean_anode /= particle->get_trajectory().get_cluster().hits().size();
      new_beta.mean_anode = cluster_mean_anodic_time[new_beta.cluster_id];

    } // for (particle)

  // printf("[%d_%d] %zd beta(s)\n", eh_run, eh_event, betas.size());

  if (betas.size() < 2)
    return dpp::base_module::PROCESS_SUCCESS;

  ////////////////////////////////
  // perform betabeta selection //
  ////////////////////////////////

  for (size_t i=0; i<betas.size(); i++) {

    const beta_data & beta_i = betas[i];

    for (size_t j=i+1; j<betas.size(); j++) {

      const beta_data & beta_j = betas[j];

      // printf("[%d_%d] checking (%zd,%zd)\n", eh_run, eh_event, i, j);

      if (beta_i.cluster_id == beta_j.cluster_id) {
	// printf("[%d_%d] skipping beta pair with same cluster\n", eh_run, eh_event);
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

      // order the beta by om num
      const beta_data & beta_1 = (beta_i.om_num < beta_j.om_num) ? beta_i : beta_j;
      const beta_data & beta_2 = (beta_i.om_num < beta_j.om_num) ? beta_j : beta_i;

      // compute source vertex distance
      const float deltay = beta_2.source_vtx[1] - beta_1.source_vtx[1];
      const float deltaz = beta_2.source_vtx[2] - beta_1.source_vtx[2];
      const float deltar2 = deltay*deltay + deltaz*deltaz;

      if (deltar2 > 0.3*0.3) {
	// printf("[%d_%d] skip betabeta candidates with large deltar\n", eh_run, eh_event);
	continue;
      }

      betabetas.push_back(betabeta_data());
      betabeta_data & new_betabeta = betabetas.back();

      // store the betabeta candidates
      new_betabeta.run = eh_run;
      new_betabeta.event = eh_event;
      new_betabeta.flag = 0;

      if (beta_1.flag & 0x1)
	new_betabeta.flag |= (1 << 0); // 0x1 = side1 bit

      if (beta_2.flag & 0x1)
	new_betabeta.flag |= (1 << 1); // 0x2 = side2 bit

      if ((beta_1.flag & 0x1) == (beta_2.flag & 0x1))
	new_betabeta.flag |= (1 << 2); // 0x4 = same side bit

      new_betabeta.cluster1 = beta_1.cluster_id;
      new_betabeta.cluster2 = beta_2.cluster_id;

      new_betabeta.particle1 = beta_1.particle_id;
      new_betabeta.particle2 = beta_2.particle_id;

      new_betabeta.om1 = beta_1.om_num;
      new_betabeta.om2 = beta_2.om_num;

      new_betabeta.nb_gg1 = beta_1.nb_cells;
      new_betabeta.nb_gg2 = beta_2.nb_cells;

      new_betabeta.nb_gg1_withz = beta_1.nb_cells_withz;
      new_betabeta.nb_gg2_withz = beta_2.nb_cells_withz;

      new_betabeta.e1 = beta_1.energy;
      new_betabeta.e2 = beta_2.energy;

      new_betabeta.t1 = beta_1.time;
      new_betabeta.t2 = beta_2.time;

      new_betabeta.l1 = beta_1.length;
      new_betabeta.l2 = beta_2.length;

      for (int xyz=0; xyz<3; xyz++) {
        new_betabeta.calo_vtx1[xyz] = beta_1.calo_vtx[xyz];
	new_betabeta.calo_vtx2[xyz] = beta_2.calo_vtx[xyz];

        new_betabeta.source_vtx[xyz] = 0.5 * (beta_1.source_vtx[xyz] + beta_2.source_vtx[xyz]);
      }

      new_betabeta.source_deltay = beta_2.source_vtx[1] - beta_1.source_vtx[1];
      new_betabeta.source_deltaz = beta_2.source_vtx[2] - beta_1.source_vtx[2];

      // new_betabeta.alphas.clear();
      // new_betabeta.gammas.clear();

    } // for (j)

  } // for (i)

  ////////////////////////////////

  int nb_betabeta = 0;

  if (betabetas.size() > 0) {

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

	if ((best_betabeta->cluster1 == new_betabeta->cluster1) &&
	    (best_betabeta->cluster2 == new_betabeta->cluster2)) {

	  if (TMath::Abs(new_betabeta->source_deltay) < TMath::Abs(best_betabeta->source_deltay)) {

	    // printf("[%d_%d] better betabeta selected (%d,%d) -> (%d,%d)\n", eh_run, eh_event,
	    // 	   best_betabeta->particle1, best_betabeta->particle2, new_betabeta->particle1, new_betabeta->particle2);

	    best_betabeta = new_betabeta;

	  } else {

	    // printf("[%d_%d] keep betabeta selected (%d,%d) <- (%d,%d)\n", eh_run, eh_event,
	    // 	   best_betabeta->particle1, best_betabeta->particle2, new_betabeta->particle1, new_betabeta->particle2);

	  }

	  skip_betabeta[bb2] = true;

	}

	else if ((best_betabeta->cluster1 == new_betabeta->cluster2) &&
		 (best_betabeta->cluster2 == new_betabeta->cluster1)) {

	  printf("[%d_%d] !!!!!!!! swapped double betabeta candidates !!!!!!\n", eh_run, eh_event);

	}

      } // for (bb2)

      // iterate over all other clusters for possible space/time
      // correlation of betabeta candidate with an extra track

      const double betabeta_mean_anode = 0.5 * (cluster_mean_anodic_time[best_betabeta->cluster1] + cluster_mean_anodic_time[best_betabeta->cluster2]);

      bool has_intime_cluster = false;
      bool has_delayed_cluster = false;

      for (size_t cluster_id=0; cluster_id<tcd_solution.get_clusters().size(); cluster_id++) {

	if (cluster_id == best_betabeta->cluster1) continue;
	if (cluster_id == best_betabeta->cluster2) continue;

	// compute deltat between betabeta and cluster
	const float deltat_bb_cluster = (cluster_mean_anodic_time[cluster_id] - betabeta_mean_anode) / CLHEP::microsecond;

	// compute deltar between betabeta and cluster
	float deltar2_bb_cluster = std::numeric_limits<float>::max();

	for (const datatools::handle<snemo::datamodel::calibrated_tracker_hit> & tracker_hit : tcd_solution.get_clusters()[cluster_id]->hits()) {

	  float tmp_deltar2_bb_cluster = 0;

	  // deltar cell - source vertex
	  tmp_deltar2_bb_cluster = 0;
	  tmp_deltar2_bb_cluster += std::pow(tracker_hit->get_x()/CLHEP::m - best_betabeta->source_vtx[0], 2);
	  tmp_deltar2_bb_cluster += std::pow(tracker_hit->get_y()/CLHEP::m - best_betabeta->source_vtx[1], 2);
	  if (tmp_deltar2_bb_cluster < deltar2_bb_cluster) deltar2_bb_cluster = tmp_deltar2_bb_cluster;

	  // deltar cell - calo1 vertex
	  tmp_deltar2_bb_cluster = 0;
	  tmp_deltar2_bb_cluster += std::pow(tracker_hit->get_x()/CLHEP::m - best_betabeta->calo_vtx1[0], 2);
	  tmp_deltar2_bb_cluster += std::pow(tracker_hit->get_y()/CLHEP::m - best_betabeta->calo_vtx1[1], 2);
	  if (tmp_deltar2_bb_cluster < deltar2_bb_cluster) deltar2_bb_cluster = tmp_deltar2_bb_cluster;

	  // deltar cell - calo2 vertex
	  tmp_deltar2_bb_cluster = 0;
	  tmp_deltar2_bb_cluster += std::pow(tracker_hit->get_x()/CLHEP::m - best_betabeta->calo_vtx2[0], 2);
	  tmp_deltar2_bb_cluster += std::pow(tracker_hit->get_y()/CLHEP::m - best_betabeta->calo_vtx2[1], 2);
	  if (tmp_deltar2_bb_cluster < deltar2_bb_cluster) deltar2_bb_cluster = tmp_deltar2_bb_cluster;

	  // ideally should compare with all betabeta cells..
	}

	const float deltar_bb_cluster = std::sqrt(deltar2_bb_cluster);

	// printf("[%d_%d] extra cluster with deltat = %.3f us and deltar = %.3f m\n", eh_run, eh_event, deltat_bb_cluster, deltar_bb_cluster);
	_output_deltart_bb_cluster_->Fill(deltat_bb_cluster, deltar_bb_cluster);

	// perform betabeta veto if cluster is correlated in time
	// if ((-2.5 < deltat_bb_cluster) && (deltat_bb_cluster < +5.0)) {
	if (std::abs(deltat_bb_cluster) < _deltat_extra_cluster_)
	  has_intime_cluster = true;

	else if (deltat_bb_cluster >= _deltat_extra_cluster_) {
	  if (deltar_bb_cluster < _deltar_delayed_cluster_)
	    has_delayed_cluster = true;
	}

      } // for (cluster_id)

      // iterate over all other calo hits for possible time
      // correlation of betabeta candidate with an extra hit

      const double betabeta_mean_calo_time = 0.5 * (best_betabeta->t1 + best_betabeta->t2);

      bool has_intime_calo = false;

      for (const datatools::handle<snemo::datamodel::calibrated_calorimeter_hit> & calo_hit : CD.calorimeter_hits()) {

	const int om_num = snemo::datamodel::om_num(calo_hit->get_geom_id());

	if (om_num == best_betabeta->om1) continue;
	if (om_num == best_betabeta->om2) continue;

	const float deltat_bb_calo = calo_hit->get_time()/CLHEP::ns - betabeta_mean_calo_time;
	// printf("[%d_%d] extra calo with deltat = %.3f ns\n", eh_run, eh_event, deltat_bb_calo);
	_output_deltat_bb_calo_->Fill(deltat_bb_calo);

	if (std::abs(deltat_bb_calo) < _deltat_extra_calo_)
	  has_intime_calo = true;
      }

      // skip any betabeta candidate with an extra cluster intime
      if (has_intime_cluster)
	continue;

      // best_betabeta->flag already initialed
      // best_betabeta->flag |= (1 << 0); // 0x1 = side1 bit
      // best_betabeta->flag |= (1 << 1); // 0x2 = side2 bit
      // best_betabeta->flag |= (1 << 2); // 0x4 = same side bit

      // if (xxxxxxxx)
      //   best_betabeta->flag |= (1 << 3); // 0x8

      if (has_intime_cluster)
        best_betabeta->flag |= (1 << 4); // 0x10 = intime cluster

      if (has_delayed_cluster)
	best_betabeta->flag |= (1 << 5); // 0x20 = delayed cluster

      if (has_intime_calo)
	best_betabeta->flag |= (1 << 6); // 0x40 = intime calo

      // add waveform flag

      std::memcpy(&_betabeta_data_, best_betabeta, sizeof(betabeta_data));

      // _betabeta_data_.run = best_betabeta->run;
      // _betabeta_data_.event = best_betabeta->event;
      // _betabeta_data_.flag = best_betabeta->flag;
      // _betabeta_data_.cluster1 = best_betabeta->cluster1;
      // _betabeta_data_.cluster2 = best_betabeta->cluster2;
      // _betabeta_data_.particle1 = best_betabeta->particle1;
      // _betabeta_data_.particle2 = best_betabeta->particle2;
      // _betabeta_data_.om1 = best_betabeta->om1;
      // _betabeta_data_.om2 = best_betabeta->om2;
      // _betabeta_data_.e1 = best_betabeta->e1;
      // _betabeta_data_.e2 = best_betabeta->e2;
      // _betabeta_data_.t1 = best_betabeta->t1;
      // _betabeta_data_.t2 = best_betabeta->t2;
      // _betabeta_data_.l1 = best_betabeta->l1;
      // _betabeta_data_.l2 = best_betabeta->l2;

      // for (int xyz=0; xyz<3; xyz++) {
      // 	_betabeta_data_.calo_vtx1[xyz] = best_betabeta->calo_vtx1[xyz];
      // 	_betabeta_data_.calo_vtx2[xyz] = best_betabeta->calo_vtx2[xyz];
      // 	_betabeta_data_.source_vtx[xyz] = best_betabeta->source_vtx[xyz];
      // }

      // _betabeta_data_.source_deltay = best_betabeta->source_deltay;
      // _betabeta_data_.source_deltaz = best_betabeta->source_deltaz;

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
  _output_deltart_bb_cluster_->Write("", TObject::kOverwrite);
  _output_deltat_bb_calo_->Write("", TObject::kOverwrite);
  _output_file_->Close();

  this->_set_initialized(false);
}

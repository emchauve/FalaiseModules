#include <bayeux/dpp/chain_module.h>

#include <falaise/snemo/datamodels/calibrated_data.h>
#include <falaise/snemo/datamodels/tracker_cluster.h>
#include <falaise/snemo/datamodels/tracker_clustering_data.h>
#include <falaise/snemo/datamodels/tracker_clustering_solution.h>
#include <falaise/snemo/datamodels/geomid_utils.h>
#include <falaise/property_set.h>

#include <cstdio>
#include <limits>

////////////////////////////////////////////////////////////////

class cd_killer_module : public dpp::chain_module
{
public:
  cd_killer_module();
  ~cd_killer_module();

  void initialize(const datatools::properties &,
		  datatools::service_manager &,
		  dpp::module_handle_dict_type &);
  
  dpp::chain_module::process_status process(datatools::things &);

  void finalize();

private:
  bool calo_dead_om[712];
  bool tracker_dead_gg[2034];

  double calo_energy_threshold;
  double tracker_max_radius;

  DPP_MODULE_REGISTRATION_INTERFACE(cd_killer_module);
};

DPP_MODULE_REGISTRATION_IMPLEMENT(cd_killer_module, "cd_killer_module");

////////////////////////////////////////////////////////////////

cd_killer_module::cd_killer_module() : dpp::chain_module()
{
  
}

cd_killer_module::~cd_killer_module()
{
  this->finalize();
}

void cd_killer_module::initialize(const datatools::properties & config,
			   datatools::service_manager & ,
			   dpp::module_handle_dict_type & )
{
  const falaise::property_set fps {config};

  memset(calo_dead_om, false, 712*sizeof(bool));
  memset(tracker_dead_gg, false, 2034*sizeof(bool));

  if (config.has_key("calo.dead_om"))
    {
      std::vector<int> dead_om_num_v;
      config.fetch("calo.dead_om", dead_om_num_v);

      printf("-> %zd dead OM: ", dead_om_num_v.size());

      for (const int dead_om_num : dead_om_num_v) {
	calo_dead_om[dead_om_num] = true;
	printf("%d ", dead_om_num);}

      printf("\n");
    }

  if (config.has_key("tracker.dead_gg"))
    {
      std::vector<int> dead_cell_num_v;
      config.fetch("tracker.dead_gg", dead_cell_num_v);

      printf("-> %zd dead CELL: ", dead_cell_num_v.size());

      for (const int dead_cell_num : dead_cell_num_v) {
	tracker_dead_gg[dead_cell_num] = true;
	printf("%d ", dead_cell_num);}

      printf("\n");
    }

  calo_energy_threshold = fps.get<double>("calo.energy_threshold", -1);

  if (calo_energy_threshold > 0)
    printf("-> calo energy threshold = %.2f MeV\n", calo_energy_threshold);

  if (config.has_key("tracker.maxr_mm")) {
    tracker_max_radius = config.fetch_real("tracker.maxr_mm") * CLHEP::mm;
    printf("-> tracker max radius = %.1f mm\n", tracker_max_radius/CLHEP::mm);
  } else {
    tracker_max_radius = std::numeric_limits<double>::max();
  }

  this->_set_initialized(true);
}

dpp::chain_module::process_status cd_killer_module::process(datatools::things &event)
{
  snemo::datamodel::calibrated_data & CD = event.grab<snemo::datamodel::calibrated_data>("CD");

  for (auto tracker_hit_it = CD.tracker_hits().begin(); tracker_hit_it != CD.tracker_hits().end() ; )
    {
      const auto & tracker_hit = tracker_hit_it->get();

      const int a_cell_num = snemo::datamodel::gg_num(tracker_hit.get_geom_id());

      if ((tracker_dead_gg[a_cell_num]) || (tracker_hit.get_r() > tracker_max_radius))
	tracker_hit_it = CD.tracker_hits().erase(tracker_hit_it);

      else
	++tracker_hit_it;
    }


  for (auto calo_hit_it = CD.calorimeter_hits().begin(); calo_hit_it != CD.calorimeter_hits().end() ; )
    {
      const auto & calo_hit = calo_hit_it->get();

      const int an_om_num = snemo::datamodel::om_num(calo_hit.get_geom_id());

      if (calo_dead_om[an_om_num])	
	calo_hit_it = CD.calorimeter_hits().erase(calo_hit_it);

      else if (calo_hit.get_energy() < calo_energy_threshold)
	calo_hit_it = CD.calorimeter_hits().erase(calo_hit_it);

      else
	++calo_hit_it;
    }

  if (event.has("TCD")) {

    snemo::datamodel::tracker_clustering_data & TCD = event.grab<snemo::datamodel::tracker_clustering_data>("TCD");

    for (auto & tcd_solution : TCD.solutions()) {
      for (auto & tcd_cluster : tcd_solution->get_clusters()) {
	for (auto tracker_hit_it = tcd_cluster->hits().begin(); tracker_hit_it != tcd_cluster->hits().end() ; ) {

	  const auto & tracker_hit = tracker_hit_it->get();

	  const int a_cell_num = snemo::datamodel::gg_num(tracker_hit.get_geom_id());

	  if ((tracker_dead_gg[a_cell_num]) || (tracker_hit.get_r() > tracker_max_radius))
	    tracker_hit_it = tcd_cluster->hits().erase(tracker_hit_it);

	  else
	    ++tracker_hit_it;

	} // for (tracker_hit_it)
      } // for (tcd_cluster)
    } // for (tcd_solution)

  } // if has("TCD")

  return dpp::base_module::PROCESS_SUCCESS;
}

void cd_killer_module::finalize()
{
  this->_set_initialized(false);
}

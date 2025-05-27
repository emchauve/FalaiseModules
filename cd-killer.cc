#include <bayeux/dpp/chain_module.h>
#include <bayeux/genbb_help/primary_event.h>
#include <bayeux/mctools/base_step_hit.h>
#include <bayeux/mctools/simulated_data.h>

#include <falaise/snemo/datamodels/event_header.h>
#include <falaise/snemo/datamodels/geomid_utils.h>
#include <falaise/snemo/datamodels/calibrated_data.h>

class cd_killer : public dpp::chain_module
{
public:
  cd_killer();
  virtual ~cd_killer();
  virtual void initialize(const datatools::properties &myConfig,
                          datatools::service_manager& flServices,
                          dpp::module_handle_dict_type& moduleDict);

  virtual dpp::chain_module::process_status process(datatools::things &event);
  
  virtual void finalize ();

private:
  DPP_MODULE_REGISTRATION_INTERFACE(cd_killer);

  bool calo_dead_om[712];
  bool tracker_dead_gg[2034];

  float calo_energy_threshold;
};

DPP_MODULE_REGISTRATION_IMPLEMENT(cd_killer, "cd_killer");

cd_killer::cd_killer() : dpp::chain_module() {}

cd_killer::~cd_killer()
{
  this->finalize();
}

void cd_killer::initialize(const datatools::properties  &myConfig,
			      datatools::service_manager   &/*&flServices*/,
			      dpp::module_handle_dict_type &/*moduleDict*/)
{
  if (myConfig.has_key("calo.dead_om"))
    {
      std::vector<int> dead_om_num_v;
      myConfig.fetch("calo.dead_om", dead_om_num_v);

      memset(calo_dead_om, false, 712*sizeof(bool));

      printf("-> %zd dead OM: ", dead_om_num_v.size());

      for (int dead_om_num : dead_om_num_v) {
	calo_dead_om[dead_om_num] = true;
	printf("%d ", dead_om_num);}

      printf("\n");
    }

  if (myConfig.has_key("tracker.dead_gg"))
    {
      std::vector<int> dead_cell_num_v;
      myConfig.fetch("tracker.dead_gg", dead_cell_num_v);

      memset(tracker_dead_gg, false, 2034*sizeof(bool));

      printf("-> %zd dead CELL: ", dead_cell_num_v.size());

      for (int dead_cell_num : dead_cell_num_v) {
	tracker_dead_gg[dead_cell_num] = true;
	printf("%d ", dead_cell_num);}

      printf("\n");
    }

  if (myConfig.has_key("calo.energy_threshold"))
    {
      calo_energy_threshold = myConfig.fetch_real("calo.energy_threshold");
      printf("-> calo energy threshold: %.1f MeV\n", calo_energy_threshold);
    }

  this->_set_initialized(true);
}

dpp::chain_module::process_status cd_killer::process(datatools::things &event)
{
  snemo::datamodel::calibrated_data & CD = event.grab<snemo::datamodel::calibrated_data>("CD");

  for (auto tracker_hit_it = CD.tracker_hits().begin(); tracker_hit_it != CD.tracker_hits().end() ; )
    {
      const auto & tracker_hit = tracker_hit_it->get();

      const int a_cell_num = snemo::datamodel::gg_num(tracker_hit.get_geom_id());

      if (tracker_dead_gg[a_cell_num])
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

  return PROCESS_SUCCESS;
}

void cd_killer::finalize()
{
  this->_set_initialized(false);
}

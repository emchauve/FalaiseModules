#include "ttd-residual-data.cc"

#include <bayeux/dpp/chain_module.h>

#include <falaise/snemo/datamodels/event_header.h>
#include <falaise/snemo/datamodels/geomid_utils.h>
#include <falaise/snemo/datamodels/line_trajectory_pattern.h>
#include <falaise/snemo/datamodels/precalibrated_data.h>
#include <falaise/snemo/datamodels/precalibrated_tracker_hit.h>
#include <falaise/snemo/datamodels/tracker_clustering_data.h>
#include <falaise/snemo/datamodels/tracker_digitized_hit.h>
#include <falaise/snemo/datamodels/tracker_trajectory_data.h>
#include <falaise/snemo/datamodels/unified_digitized_data.h>
#include <falaise/property_set.h>

#include <cstdio>
#include <cmath>

////////////////////////////////////////////////////////////////

class ttd_residual_module : public dpp::chain_module
{
public:
  ttd_residual_module();
  ~ttd_residual_module();

  void initialize(const datatools::properties &,
		  datatools::service_manager &,
		  dpp::module_handle_dict_type &);
  
  dpp::chain_module::process_status process(datatools::things &);

  void finalize();

private:
  double distance_line_point_2d(double, double, double, double);
  double distance_line_point_2d(double, double, double, double, double);

private:
  std::string _ttd_label_;

  ttd_residual_data _ttd_data_;

  DPP_MODULE_REGISTRATION_INTERFACE(ttd_residual_module);
};

DPP_MODULE_REGISTRATION_IMPLEMENT(ttd_residual_module, "ttd_residual_module");

////////////////////////////////////////////////////////////////

ttd_residual_module::ttd_residual_module() : dpp::chain_module()
{
  
}

ttd_residual_module::~ttd_residual_module()
{
  this->finalize();
}

void ttd_residual_module::initialize(const datatools::properties & config,
			   datatools::service_manager & ,
			   dpp::module_handle_dict_type & )
{
  const falaise::property_set fps {config};

  _ttd_label_ = fps.get<std::string>("TTD_label", "TTD");

  this->_set_initialized(true);
}

dpp::chain_module::process_status ttd_residual_module::process(datatools::things &event)
{
  if (!event.has("TTD"))
    {
      std::cout << "*** no TTD bank ***" << std::endl;
      return dpp::base_module::PROCESS_ERROR;
    }

  const snemo::datamodel::tracker_trajectory_data & TTD = event.get<snemo::datamodel::tracker_trajectory_data>("TTD");

  if (TTD.get_solutions().size() == 0)
    return dpp::base_module::PROCESS_SUCCESS;

  const snemo::datamodel::event_header & EH = event.get<snemo::datamodel::event_header>("EH");
  _ttd_data_.run = EH.get_id().get_run_number();
  _ttd_data_.event = EH.get_id().get_event_number();

  const snemo::datamodel::unified_digitized_data & UDD = event.get<snemo::datamodel::unified_digitized_data>("UDD");
  const snemo::datamodel::precalibrated_data & pCD = event.get<snemo::datamodel::precalibrated_data>("pCD");
  const snemo::datamodel::tracker_clustering_data & TCD = event.get<snemo::datamodel::tracker_clustering_data>("TCD");


  const snemo::datamodel::tracker_trajectory_solution & ttd_solution = TTD.get_default_solution();

  for (const datatools::handle<snemo::datamodel::tracker_trajectory> & ttd_trajectory : ttd_solution.get_trajectories()) {

    if (! ttd_trajectory->get_fit_infos().is_best())
      continue;

    const snemo::datamodel::base_trajectory_pattern & tracker_pattern = ttd_trajectory->get_pattern();

    // process only on line fit
    if (tracker_pattern.get_pattern_id() != snemo::datamodel::line_trajectory_pattern::pattern_id())
      continue;

    const auto & line_trajectory = dynamic_cast<const snemo::datamodel::line_trajectory_pattern &> (tracker_pattern);

    const geomtools::vector_3d & line_start_point = line_trajectory.get_first();
    const geomtools::vector_3d & line_stop_point = line_trajectory.get_last();

    // 2D solution in X/Y plane (top view) [a.x + b.y + c = 0]
    const double line2d_xy_a = line_stop_point.y() - line_start_point.y();
    const double line2d_xy_b = line_start_point.x() - line_stop_point.x();
    const double line2d_xy_c = -line2d_xy_a*line_start_point.x() -line2d_xy_b*line_start_point.y();

    // 2D solution in X/Z plane (side view)
    const double line2d_xz_a = line_stop_point.z() - line_start_point.z();
    const double line2d_xz_b = line_start_point.x() - line_stop_point.x();
    // const double line2d_xz_c = -line2d_xy_a*line_start_point.x() -line2d_xy_b*line_start_point.z();
    const double line2d_xz_slope = (line_stop_point.z()-line_start_point.z())/(line_stop_point.x()-line_start_point.x());
    const double line2d_xz_intercept = line_start_point.z() - line2d_xz_slope*line_start_point.x();

    const double phi_from_fit = -std::atan2(line2d_xy_a, line2d_xy_b)/CLHEP::degree;
    const double theta_from_fit = -std::atan2(line2d_xz_a, line2d_xz_b)/CLHEP::degree;

    printf("[%d_%d:%d]   phi %.1f   theta = %.1f\n", _ttd_data_.run, _ttd_data_.event,
	   ttd_trajectory->get_id(), phi_from_fit, theta_from_fit);

    uint32_t trajectory_side = 0;

    for (const auto & tracker_hit : ttd_trajectory->get_cluster().hits()) {

      trajectory_side = tracker_hit->get_geom_id().get(1);

      // retrieve cell center (anode)
      const double & tracker_hit_x = tracker_hit->get_x();
      const double & tracker_hit_y = tracker_hit->get_y();

      const double radius_from_cd = tracker_hit->get_r();
      const double radius_from_fit = this->distance_line_point_2d(line2d_xy_a, line2d_xy_b, line2d_xy_c, tracker_hit_x, tracker_hit_y); 
      const double radius_residual = radius_from_fit - radius_from_cd;

      const double z_from_cd = tracker_hit->get_z();
      const double z_from_fit = line2d_xz_intercept + line2d_xz_slope*tracker_hit_x;
      const double z_residual = z_from_fit - z_from_cd;
      
      // flag on cell's neighbour configuration
      //  9 1 2
      //  8 X 3
      //  7 6 5 

      printf("%s   deltaR = %5.2f mm", snemo::datamodel::gg_label(tracker_hit->get_geom_id()).c_str(), radius_residual/CLHEP::mm);

      if (datatools::is_valid(z_residual))
	printf("   deltaZ = %5.2f mm", z_residual/CLHEP::mm);

      printf("\n");

    } // for (tracker_hit)

    // Now loop over all unclustered cells to check possible missed cells
    // by the clusterisation algo (cells crossed by the fitted trajectory).

    for (const auto & tracker_hit : TCD.get_default().get_unclustered_hits()) {

      // consider only tracker_hit on same side as the trajectory
      if (tracker_hit->get_geom_id().get(1) != trajectory_side)
	continue;

      // retrieve cell center (anode)
      const double & tracker_hit_x = tracker_hit->get_x();
      const double & tracker_hit_y = tracker_hit->get_y();

      const double radius_from_fit = this->distance_line_point_2d(line2d_xy_a, line2d_xy_b, line2d_xy_c, tracker_hit_x, tracker_hit_y); 

      if (std::fabs(radius_from_fit) > 4.4*std::sqrt(2.0)*CLHEP::cm)
	continue;

      const double radius_from_cd = tracker_hit->get_r();
      const double radius_residual = radius_from_fit - radius_from_cd;

      const double z_from_cd = tracker_hit->get_z();
      const double z_from_fit = line2d_xz_intercept + line2d_xz_slope*tracker_hit_x;
      const double z_residual = z_from_fit - z_from_cd;

      std::cout << "/!\\ " << snemo::datamodel::gg_label(tracker_hit->get_geom_id())
		<< "  deltaR = " << radius_residual/CLHEP::mm
		<< "  deltaZ = " << z_residual/CLHEP::mm << std::endl;

      const datatools::properties & tracker_hit_properties = tracker_hit->get_auxiliaries();

      const int UDD_parent_index = tracker_hit_properties.fetch_integer("UDD.parent");
      const int pCD_parent_index = tracker_hit_properties.fetch_integer("pCD.parent");

      const snemo::datamodel::tracker_digitized_hit & udd_tracker_hit = UDD.get_tracker_hits().at(UDD_parent_index).get();
      const snemo::datamodel::precalibrated_tracker_hit & pcd_tracker_hit = pCD.tracker_hits().at(pCD_parent_index).get();

      std::cout << "|- UDD[" << UDD_parent_index << "]  bottom = " << udd_tracker_hit.get_times().front().get_bottom_cathode_time()
		<< "  top = " << udd_tracker_hit.get_times().front().get_top_cathode_time() << std::endl;

      std::cout << "`- pCD[" << pCD_parent_index << "]  bottom = " << pcd_tracker_hit.get_bottom_cathode_drift_time()/CLHEP::microsecond
		<< "  top = " << pcd_tracker_hit.get_top_cathode_drift_time()/CLHEP::microsecond << std::endl;

    } // for (unclustered tracker_hit)

  } // for (ttd_trajectory)

  return dpp::base_module::PROCESS_SUCCESS;
}

void ttd_residual_module::finalize()
{
  this->_set_initialized(false);
}

double ttd_residual_module::distance_line_point_2d(double line_p0, double line_p1, double point_x, double point_y)
{
  return std::fabs(line_p1*point_x - point_y + line_p0) / std::sqrt(1 + line_p1*line_p1);
}

double ttd_residual_module::distance_line_point_2d(double line_a, double line_b, double line_c, double point_x, double point_y)
{
  return std::fabs(line_a*point_x + line_b*point_y + line_c) / std::sqrt(line_a*line_a + line_b*line_b);
}

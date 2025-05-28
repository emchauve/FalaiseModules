#include "ttd-residual-data.cc"

#include <bayeux/dpp/chain_module.h>

#include <falaise/snemo/datamodels/event_header.h>
#include <falaise/snemo/datamodels/geomid_utils.h>
#include <falaise/snemo/datamodels/line_trajectory_pattern.h>
#include <falaise/snemo/datamodels/precalibrated_tracker_hit.h>
#include <falaise/snemo/datamodels/tracker_clustering_data.h>
#include <falaise/snemo/datamodels/tracker_digitized_hit.h>
#include <falaise/snemo/datamodels/tracker_trajectory_data.h>
#include <falaise/property_set.h>

#include <cstdio>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"

////////////////////////////////////////////////////////////////

void ttd_residual_data_print(const ttd_residual_data & rd)
{
  printf("[%d_%d] with %zd track(s)\n", rd.run, rd.event, rd.track.size());

  for (const track_data & tr : rd.track)
    {
      printf("- nb_gg=%3zd  flag=%02d  length=%6.1f  theta=%6.1f  phi=%5.1f  chi2/ndf=%5.2f\n",
	     tr.cells.size(), tr.flag, tr.length, tr.theta, tr.phi, tr.chi2ndf);
    }
}

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
  std::string _output_filename_;

  ttd_residual_data _ttd_data_;

  TFile *_output_file_;
  TTree *_output_tree_;

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
  _output_filename_ = fps.get<std::string>("output_filename", "output.root");

  // prepare ROOT output file/tree
  _output_file_ = new TFile(_output_filename_.c_str(), "RECREATE");
  _output_tree_ = new TTree("ttd_residual", "");
  _output_tree_->Branch("ttd", &_ttd_data_);
  gROOT->cd();

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

  if (!TTD.has_default_solution())
    return dpp::base_module::PROCESS_SUCCESS;

  const snemo::datamodel::event_header & EH = event.get<snemo::datamodel::event_header>("EH");
  _ttd_data_.run = EH.get_id().get_run_number();
  _ttd_data_.event = EH.get_id().get_event_number();

  const snemo::datamodel::tracker_clustering_data & TCD = event.get<snemo::datamodel::tracker_clustering_data>("TCD");

  const snemo::datamodel::tracker_trajectory_solution & ttd_solution = TTD.get_default_solution();
  const snemo::datamodel::tracker_clustering_solution & tcd_solution = TCD.get_default();

  ////////////////////////////////

  // loop over clusters and prepare the cluster_to_ttd map

  std::vector<std::vector<int>> cluster_to_ttd;

  const int nb_clusters = tcd_solution.get_clusters().size();

  for (int c=0; c<nb_clusters; c++)
    {
      std::vector<int> a_vector;
      cluster_to_ttd.push_back(a_vector);
    }

  ////////////////////////////////

  // iterate over trajectories to identify multiple trajectory for single cluster;

  for (const datatools::handle<snemo::datamodel::tracker_trajectory> & ttd_trajectory : ttd_solution.get_trajectories())
    {
      // process only with best trajectory
      if (! ttd_trajectory->get_fit_infos().is_best())
	continue;

      // process only with line fit trajectory
      const snemo::datamodel::base_trajectory_pattern & tracker_pattern = ttd_trajectory->get_pattern();
      if (tracker_pattern.get_pattern_id() != snemo::datamodel::line_trajectory_pattern::pattern_id())
	continue;

      const snemo::datamodel::tracker_cluster & ttd_cluster = ttd_trajectory->get_cluster();
      cluster_to_ttd[ttd_cluster.get_cluster_id()].push_back(ttd_trajectory->get_id());
    }

  // for (size_t cluster_id=0; cluster_id<cluster_to_ttd.size(); cluster_id++)
  //   {
  //     if (cluster_to_ttd[cluster_id].size() > 1)
  // 	printf("[%d_%d] multiple trajectories for cluster %zd\n", _ttd_data_.run, _ttd_data_.event, cluster_id);

  //     if (cluster_to_ttd[cluster_id].size() > 2)
  // 	printf("[%d_%d] >2 trajectories for cluster %zd !\n", _ttd_data_.run, _ttd_data_.event, cluster_id);
  //   }

  ////////////////////////////////

  for (const datatools::handle<snemo::datamodel::tracker_trajectory> & ttd_trajectory : ttd_solution.get_trajectories()) {

    // process only with best trajectory
    if (! ttd_trajectory->get_fit_infos().is_best())
      continue;

    // process only with line fit trajectory
    const snemo::datamodel::base_trajectory_pattern & tracker_pattern = ttd_trajectory->get_pattern();
    if (tracker_pattern.get_pattern_id() != snemo::datamodel::line_trajectory_pattern::pattern_id())
      continue;

    const auto & line_trajectory = dynamic_cast<const snemo::datamodel::line_trajectory_pattern &> (tracker_pattern);

    const geomtools::vector_3d & line_first_point = line_trajectory.get_first();
    const geomtools::vector_3d & line_last_point = line_trajectory.get_last();

    // 2D solution in X/Y plane (top view) [a.x + b.y + c = 0]
    const double line2d_xy_a = line_last_point.y() - line_first_point.y();
    const double line2d_xy_b = line_first_point.x() - line_last_point.x();
    const double line2d_xy_c = -line2d_xy_a*line_first_point.x() -line2d_xy_b*line_first_point.y();

    // 2D solution in X/Z plane (side view)
    // const double line2d_xz_a = line_last_point.z() - line_first_point.z();
    // const double line2d_xz_b = line_first_point.x() - line_last_point.x();
    // const double line2d_xz_c = -line2d_xy_a*line_first_point.x() -line2d_xy_b*line_first_point.z();
    const double line2d_xz_slope = (line_last_point.z()-line_first_point.z())/(line_last_point.x()-line_first_point.x());
    const double line2d_xz_intercept = line_first_point.z() - line2d_xz_slope*line_first_point.x();

    const geomtools::vector_3d & line_vector = line_last_point - line_first_point;
    const double theta_from_fit = line_vector.theta()/CLHEP::degree;
    double phi_from_fit = line_vector.phi()/CLHEP::degree;
    if (phi_from_fit < 0) phi_from_fit += 180;

    // prepare new track_data entry and fill it
    _ttd_data_.track.push_back(track_data());
    track_data & _track_data_ = _ttd_data_.track.back();

    _track_data_.flag = 0;
    for (int xyz=0; xyz<3; xyz++)
      {
	_track_data_.first[xyz] = line_first_point[xyz];
	_track_data_.last[xyz] = line_last_point[xyz];
	_track_data_.length += std::pow(line_last_point[xyz] - line_first_point[xyz], 2);
      }
    _track_data_.length = std::sqrt(_track_data_.length);
    _track_data_.theta = theta_from_fit;
    _track_data_.phi = phi_from_fit;

    _track_data_.chi2ndf = ttd_trajectory->get_fit_infos().get_chi2();
    _track_data_.chi2ndf /= ttd_trajectory->get_fit_infos().get_ndof();

    // printf("[%d_%d:%d]   phi %.1f   theta = %.1f\n", _ttd_data_.run, _ttd_data_.event,
    // 	   ttd_trajectory->get_id(), phi_from_fit, theta_from_fit);

    uint32_t trajectory_side = 0;

    for (const auto & tracker_hit : ttd_trajectory->get_cluster().hits()) {

      trajectory_side = tracker_hit->get_geom_id().get(1);

      // retrieve cell center (anode)
      const double & tracker_hit_x = tracker_hit->get_x();
      const double & tracker_hit_y = tracker_hit->get_y();

      const double radius_from_cd = tracker_hit->get_r();
      const double radius_from_fit = this->distance_line_point_2d(line2d_xy_a, line2d_xy_b, line2d_xy_c, tracker_hit_x, tracker_hit_y); 
      // const double radius_residual = radius_from_fit - radius_from_cd;

      const double height_from_cd = tracker_hit->get_z();
      const double height_from_fit = line2d_xz_intercept + line2d_xz_slope*tracker_hit_x; // FALSE!!
      // const double height_residual = height_from_fit - height_from_cd;

      const datatools::properties & tracker_hit_properties = tracker_hit->get_auxiliaries();

      // prepare new cell_data entry and fill it
      _track_data_.cells.push_back(cell_data());
      cell_data & _cell_data_ = _track_data_.cells.back();

      _cell_data_.flag = 0;
      _cell_data_.gg_num = snemo::datamodel::gg_num(tracker_hit->get_geom_id());
      _cell_data_.time_anode = tracker_hit->get_anode_time();

      if (tracker_hit_properties.has_key("bottom_drift_time")) {
	_cell_data_.time_bottom_cathode = tracker_hit_properties.fetch_real("bottom_drift_time");
	_cell_data_.flag |= (1 << 0);
      } else _cell_data_.time_bottom_cathode = 0;

      if (tracker_hit_properties.has_key("top_drift_time")) {
	_cell_data_.time_top_cathode = tracker_hit_properties.fetch_real("top_drift_time");
	_cell_data_.flag |= (1 << 1);
      } else _cell_data_.time_top_cathode = 0;

      _cell_data_.r = radius_from_cd;
      _cell_data_.z = height_from_cd;
      _cell_data_.rfit = radius_from_fit;
      _cell_data_.zfit = height_from_fit;

      // printf("[%d_%d] gg %4d : anode = %5.2f us   bottom = %5.1f us   top = %5.1f us   r = %5.2f mm   z = %5.2f cm\n",
      // 	     _ttd_data_.run, _ttd_data_.event, _cell_data_.gg_num, _cell_data_.time_anode/CLHEP::microsecond,
      // 	     _cell_data_.time_bottom_cathode/CLHEP::microsecond, _cell_data_.time_top_cathode/CLHEP::microsecond,
      // 	     _cell_data_.r/CLHEP::mm, _cell_data_.z/CLHEP::cm);

      // flag on cell's neighbour configuration
      //  9 1 2
      //  8 X 3
      //  7 6 5 

      // printf("%s   deltaR = %5.2f mm", snemo::datamodel::gg_label(tracker_hit->get_geom_id()).c_str(), radius_residual/CLHEP::mm);
      // if (datatools::is_valid(height_residual)) printf("   deltaZ = %5.2f mm", height_residual/CLHEP::mm);
      // printf("\n");

    } // for (tracker_hit)

    // flag the track side and convert theta to get 0
    // for horizontal tracks and [0,+90] (or [0,-90])
    // for track going toward top (or bottom) in the
    // point of view source => calo direction

    if (trajectory_side > 0) {
      _track_data_.flag |= (1 << 0);
      _track_data_.theta = 90 - _track_data_.theta;
    } else {
      _track_data_.theta = _track_data_.theta - 90;
    }

    // convert phi to get [-90;0;+90] for track
    // going the calo left/front/right in top view
    // in the point of view source => calo direction

    if (_track_data_.phi > 90)
      _track_data_.phi = 180 - _track_data_.phi;
    else
      _track_data_.phi = -_track_data_.phi;

    if (ttd_trajectory->get_cluster().is_delayed())
      _track_data_.flag |= (2 << 0);

    // Now loop over all unclustered cells to check possible missed cells
    // by the clusterisation algo (cells crossed by the fitted trajectory).

    // for (const auto & tracker_hit : TCD.get_default().get_unclustered_hits()) {

    //   // consider only tracker_hit on same side as the trajectory
    //   if (tracker_hit->get_geom_id().get(1) != trajectory_side)
    // 	continue;

    //   // retrieve cell center (anode)
    //   const double & tracker_hit_x = tracker_hit->get_x();
    //   const double & tracker_hit_y = tracker_hit->get_y();

    //   const double radius_from_fit = this->distance_line_point_2d(line2d_xy_a, line2d_xy_b, line2d_xy_c, tracker_hit_x, tracker_hit_y); 

    //   if (std::fabs(radius_from_fit) > 4.4*std::sqrt(2.0)*CLHEP::cm)
    // 	continue;

    //   const double radius_from_cd = tracker_hit->get_r();
    //   const double radius_residual = radius_from_fit - radius_from_cd;

    //   const double z_from_cd = tracker_hit->get_z();
    //   const double z_from_fit = line2d_xz_intercept + line2d_xz_slope*tracker_hit_x;
    //   const double z_residual = z_from_fit - z_from_cd;

    //   std::cout << "/!\\ " << snemo::datamodel::gg_label(tracker_hit->get_geom_id())
    // 		<< "  deltaR = " << radius_residual/CLHEP::mm
    // 		<< "  deltaZ = " << z_residual/CLHEP::mm << std::endl;

    //   const datatools::properties & tracker_hit_properties = tracker_hit->get_auxiliaries();

    // } // for (unclustered tracker_hit)

  } // for (ttd_trajectory)

  // ttd_residual_data_print(_ttd_data_);

  _output_tree_->Fill();
  _ttd_data_.track.clear();

  return dpp::base_module::PROCESS_SUCCESS;
}

void ttd_residual_module::finalize()
{
  _output_file_->cd();
  _output_tree_->Write("", TObject::kOverwrite);
  _output_file_->Close();

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

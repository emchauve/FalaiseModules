#include "ttd-residual-data.cc"

#include <bayeux/dpp/chain_module.h>
#include <bayeux/datatools/service_manager.h>

#include <falaise/snemo/datamodels/event_header.h>
#include <falaise/snemo/datamodels/geomid_utils.h>
#include <falaise/snemo/datamodels/calibrated_calorimeter_hit.h>
#include <falaise/snemo/datamodels/calibrated_data.h>
#include <falaise/snemo/datamodels/line_trajectory_pattern.h>
#include <falaise/snemo/datamodels/precalibrated_tracker_hit.h>
#include <falaise/snemo/datamodels/tracker_clustering_data.h>
#include <falaise/snemo/datamodels/tracker_digitized_hit.h>
#include <falaise/snemo/datamodels/tracker_trajectory_data.h>
#include <falaise/snemo/services/geometry.h>
#include <falaise/snemo/services/service_handle.h>
#include <falaise/property_set.h>

#include <cstdio>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"

////////////////////////////////////////////////////////////////

void ttd_residual_data_print(const ttd_residual_data & rd)
{
  printf("[%d_%d] with %zd track(s)\n", rd.run, rd.event, rd.tracks.size());

  for (const track_data & tr : rd.tracks)
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
  void nearest_line_point_2d(double, double, double, double, double, double &, double &);
  void nearest_line_point_2d(const geomtools::vector_3d &, const geomtools::vector_3d &, const double &, const double &, double &, double &);
  void nearest_line_point_3d(const geomtools::vector_3d &, const geomtools::vector_3d &, const double &, const double &, geomtools::vector_3d &);

  bool intersect_with_src_plane(const geomtools::vector_3d &, const geomtools::vector_3d &, geomtools::vector_3d &);
  bool intersect_with_mw_plane(const geomtools::vector_3d &, const geomtools::vector_3d &, const int &, geomtools::vector_3d &);
  bool intersect_with_xw0_plane(const geomtools::vector_3d &, const geomtools::vector_3d &, const int &, geomtools::vector_3d &);
  bool intersect_with_xw1_plane(const geomtools::vector_3d &, const geomtools::vector_3d &, const int &, geomtools::vector_3d &);

private:
  snemo::service_handle<snemo::geometry_svc> _geo_manager_ {};

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
			   datatools::service_manager & services,
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

  _geo_manager_ = snemo::service_handle<snemo::geometry_svc>{services};

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

  if ((_ttd_data_.event % 10000) == 0)
    printf("[%d_%d]\n", _ttd_data_.run, _ttd_data_.event);

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
  //     //     if (cluster_to_ttd[cluster_id].size() > 1)
  //     // 	printf("[%d_%d] multiple trajectories for cluster %zd\n", _ttd_data_.run, _ttd_data_.event, cluster_id);

  //     if (cluster_to_ttd[cluster_id].size() > 2)
  // 	printf("[%d_%d] %zd trajectories for cluster %zd !\n", _ttd_data_.run, _ttd_data_.event, cluster_to_ttd[cluster_id].size(), cluster_id);
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
    // const double line2d_xz_slope = (line_last_point.z()-line_first_point.z())/(line_last_point.x()-line_first_point.x());
    // const double line2d_xz_intercept = line_first_point.z() - line2d_xz_slope*line_first_point.x();

    const geomtools::vector_3d & line_vector = line_last_point - line_first_point;

    const double theta_from_fit = line_vector.theta()/CLHEP::degree;
    double phi_from_fit = line_vector.phi()/CLHEP::degree;
    if (phi_from_fit < 0) phi_from_fit += 180;

    // prepare new track_data entry and fill it
    _ttd_data_.tracks.push_back(track_data());
    track_data & _track_data_ = _ttd_data_.tracks.back();

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

    uint32_t cluster_row_min = std::numeric_limits<uint32_t>::max();
    uint32_t cluster_row_max = std::numeric_limits<uint32_t>::min();

    uint32_t cluster_layer_min = std::numeric_limits<uint32_t>::max();
    uint32_t cluster_layer_max = std::numeric_limits<uint32_t>::min();

    bool cluster_has_cell[2034];
    memset(cluster_has_cell, 0, sizeof(cluster_has_cell));

    for (const auto & tracker_hit : ttd_trajectory->get_cluster().hits()) {

      trajectory_side = tracker_hit->get_geom_id().get(1);

      if (tracker_hit->get_geom_id().get(2) < cluster_layer_min)
	cluster_layer_min = tracker_hit->get_geom_id().get(2);
      if (tracker_hit->get_geom_id().get(2) > cluster_layer_max)
	cluster_layer_max = tracker_hit->get_geom_id().get(2);

      if (tracker_hit->get_geom_id().get(3) < cluster_row_min)
	cluster_row_min = tracker_hit->get_geom_id().get(3);
      if (tracker_hit->get_geom_id().get(3) > cluster_row_max)
	cluster_row_max = tracker_hit->get_geom_id().get(3);

      const uint16_t gg_num = snemo::datamodel::gg_num(tracker_hit->get_geom_id());
      cluster_has_cell[gg_num] = true;

      // retrieve cell center (anode)
      const double & tracker_hit_x = tracker_hit->get_x();
      const double & tracker_hit_y = tracker_hit->get_y();

      geomtools::vector_3d nearest_line_point;
      this->nearest_line_point_3d(line_first_point, line_last_point, tracker_hit_x, tracker_hit_y, nearest_line_point);

      const double radius_from_cd = tracker_hit->get_r();
      const double radius_from_fit = std::sqrt(std::pow(nearest_line_point.x()-tracker_hit_x,2)+std::pow(nearest_line_point.y()-tracker_hit_y,2));
      // this->distance_line_point_2d(line2d_xy_a, line2d_xy_b, line2d_xy_c, tracker_hit_x, tracker_hit_y);

      // const double radius_residual = radius_from_fit - radius_from_cd;

      const double height_from_cd = tracker_hit->get_z();
      const double height_from_fit = nearest_line_point.z(); // line2d_xz_intercept + line2d_xz_slope*tracker_hit_x; // FALSE!!
      // const double height_residual = height_from_fit - height_from_cd;

      const datatools::properties & tracker_hit_properties = tracker_hit->get_auxiliaries();

      // prepare new cell_data entry and fill it
      _track_data_.cells.push_back(cell_data());
      cell_data & _cell_data_ = _track_data_.cells.back();

      _cell_data_.flag = 0;
      _cell_data_.gg_num = gg_num;
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
      _cell_data_.xfit = nearest_line_point.x() - tracker_hit_x;
      _cell_data_.yfit = nearest_line_point.y() - tracker_hit_y;
      _cell_data_.rfit = radius_from_fit;
      _cell_data_.zfit = height_from_fit;

      // printf("[%d_%d] gg %4d : anode = %5.2f us   bottom = %5.1f us   top = %5.1f us   r = %5.2f mm   rfit = %5.2f mm   z = %5.2f cm   zfit = %5.2f cm\n",
      // 	     _ttd_data_.run, _ttd_data_.event, _cell_data_.gg_num, _cell_data_.time_anode/CLHEP::microsecond,
      // 	     _cell_data_.time_bottom_cathode/CLHEP::microsecond, _cell_data_.time_top_cathode/CLHEP::microsecond,
      // 	     _cell_data_.r/CLHEP::mm, _cell_data_.rfit/CLHEP::mm, _cell_data_.z/CLHEP::cm, _cell_data_.zfit/CLHEP::cm);

      // flag on cell's neighbour configuration
      //  9 1 2
      //  8 X 3
      //  7 6 5 

      // printf("%s   deltaR = %5.2f mm", snemo::datamodel::gg_label(tracker_hit->get_geom_id()).c_str(), radius_residual/CLHEP::mm);
      // if (datatools::is_valid(height_residual)) printf("   deltaZ = %5.2f mm", height_residual/CLHEP::mm);
      // printf("\n");

    } // for (tracker_hit)

    // scan the cluster box area and look for untriggered cells

    for (uint32_t cell_row=cluster_row_min; cell_row<=cluster_row_max; cell_row++) {

      for (uint32_t cell_layer=cluster_layer_min; cell_layer<=cluster_layer_max; cell_layer++) {

	const uint16_t gg_num = trajectory_side*1017 + cell_row*9 + cell_layer;

	if (cluster_has_cell[gg_num])
	  continue;

	// build the cell geom ID to retrieve its X and Y location
	// => TODO: optimise with a builtin map[2034] at initialisation
	const geomtools::geom_id cell_geom_id (1204, 0, trajectory_side, cell_layer, cell_row);
	const geomtools::mapping & mapping = _geo_manager_->get_mapping();
	const geomtools::geom_info & cell_ginfo = mapping.get_geom_info(cell_geom_id);
	const geomtools::placement & cell_placement = cell_ginfo.get_world_placement();
	const geomtools::vector_3d & cell_pos  = cell_placement.get_translation();
	const double & tracker_hit_x = cell_pos.getX();
	const double & tracker_hit_y = cell_pos.getY();

	// double near_point_x, near_point_y;
	// this->nearest_line_point_2d(line_first_point, line_last_point, tracker_hit_x, tracker_hit_y, near_point_x, near_point_y);

	geomtools::vector_3d nearest_line_point_missed;
	this->nearest_line_point_3d(line_first_point, line_last_point, tracker_hit_x, tracker_hit_y, nearest_line_point_missed);

	const double nearest_delta_x = nearest_line_point_missed.x() - tracker_hit_x;
	const double nearest_delta_y = nearest_line_point_missed.y() - tracker_hit_y;

	const bool within_x = (std::fabs(nearest_delta_x) < 22.0*CLHEP::mm);
	const bool within_y = (std::fabs(nearest_delta_y) < 22.0*CLHEP::mm);

	if (within_x && within_y) {

	  // printf("[%d_%d] missed gg %d/%d/%d\n", _ttd_data_.run, _ttd_data_.event, trajectory_side, cell_row, cell_layer);

	  _track_data_.missed_cells.push_back(cell_data());
	  cell_data & _missed_cell_data_ = _track_data_.missed_cells.back();
	  _missed_cell_data_.gg_num = gg_num;
	  _missed_cell_data_.flag = 0x10; // missed-in
	  _missed_cell_data_.xfit = nearest_line_point_missed.x() - tracker_hit_x;
	  _missed_cell_data_.yfit = nearest_line_point_missed.y() - tracker_hit_y;
	  _missed_cell_data_.rfit = std::sqrt(std::pow(nearest_delta_x,2) + std::pow(nearest_delta_y,2));
	  _missed_cell_data_.zfit = nearest_line_point_missed.z();
	}

	else {
	  // check for cases where the near point is outside the cell  while the track is still crossing it !

	  // case 1
	  const double intersect_ylow_x = (-line2d_xy_c - line2d_xy_b*(tracker_hit_y-22.0*CLHEP::mm))/line2d_xy_a;
	  const double intersect_ylow_deltax = intersect_ylow_x - tracker_hit_x;
	  const bool intersect_ylow = (std::fabs(intersect_ylow_deltax) < 22.0*CLHEP::mm);

	  // case 2
	  const double intersect_yup_x = (-line2d_xy_c - line2d_xy_b*(tracker_hit_y+22.0*CLHEP::mm))/line2d_xy_a;
	  const double intersect_yup_deltax = intersect_yup_x - tracker_hit_x;
	  const bool intersect_yup = (std::fabs(intersect_yup_deltax) < 22.0*CLHEP::mm);

	  // case 3
	  const double intersect_xlow_y = (-line2d_xy_c - line2d_xy_a*(tracker_hit_x-22.0*CLHEP::mm))/line2d_xy_b;
	  const double intersect_xlow_deltay = intersect_xlow_y - tracker_hit_y;
	  const bool intersect_xlow = (std::fabs(intersect_xlow_deltay) < 22.0*CLHEP::mm);

	  // case 4
	  const double intersect_xup_y = (-line2d_xy_c - line2d_xy_a*(tracker_hit_x+22.0*CLHEP::mm))/line2d_xy_b;
	  const double intersect_xup_deltay = intersect_xup_y - tracker_hit_y;
	  const bool intersect_xup = (std::fabs(intersect_xup_deltay) < 22.0*CLHEP::mm);

	  const bool intersect = intersect_xlow || intersect_xup || intersect_ylow || intersect_yup;

	  if (intersect) {

	    if (intersect_xlow && intersect_xup)
	      printf("*** intersect xlow AND xup (this should not happend)\n");

	    if (intersect_ylow && intersect_yup)
	      printf("*** intersect ylow AND yup (this should not happend)\n");

	    int intersect_deltaxy_min_case = 0;
	    double intersect_deltaxy_min = 0;

	    if (intersect_ylow) {
	      intersect_deltaxy_min = intersect_ylow_deltax;
	      intersect_deltaxy_min_case = 1;
	    } else {
	      intersect_deltaxy_min = intersect_yup_deltax;
	      intersect_deltaxy_min_case = 2;
	    }

	    if (intersect_xlow) {
	      if (intersect_xlow_deltay < intersect_deltaxy_min) {
		intersect_deltaxy_min = intersect_xlow_deltay;
		intersect_deltaxy_min_case = 3;
	      }
	    } else {
	      if (intersect_xup_deltay < intersect_deltaxy_min) {
		intersect_deltaxy_min = intersect_xup_deltay;
		intersect_deltaxy_min_case = 4;
	      }
	    }

	    switch (intersect_deltaxy_min_case) {

	    case 1: // ylow
	      // using 2D line equation in X/Y plane (top view)
	      // [a.x + b.y + c = 0]  =>  x = (b.y + c)/a
	      nearest_line_point_missed.setY(tracker_hit_y-22.0*CLHEP::mm);
	      nearest_line_point_missed.setX((line2d_xy_b * nearest_line_point_missed.y() + line2d_xy_c) / line2d_xy_a);
	      break;

	    case 2: // yup
	      nearest_line_point_missed.setY(tracker_hit_y+22.0*CLHEP::mm);
	      nearest_line_point_missed.setX((line2d_xy_b * nearest_line_point_missed.y() + line2d_xy_c) / line2d_xy_a);
	      break;

	    case 3: // xlow
	      // 2D solution in X/Y plane (top view)
	      // [a.x + b.y + c = 0]  =>  y = (a.x + c) / b
	      nearest_line_point_missed.setX(tracker_hit_x-22.0*CLHEP::mm);
	      nearest_line_point_missed.setY((line2d_xy_a * nearest_line_point_missed.x() + line2d_xy_c) / line2d_xy_b);
	      break;

	    case 4: // xup
	      nearest_line_point_missed.setX(tracker_hit_x+22.0*CLHEP::mm);
	      nearest_line_point_missed.setY((line2d_xy_a * nearest_line_point_missed.x() + line2d_xy_c) / line2d_xy_b);
	      break;
	    }

	    const double line_delta_x = line_last_point.x() - line_first_point.x();
	    const double line_delta_z = line_last_point.z() - line_first_point.z();
	    nearest_line_point_missed.setZ(line_first_point.z() + line_delta_z * ((nearest_line_point_missed.x()-line_first_point.x())/line_delta_x));
	    // // nearest_line_point_missed;

	    // printf("[%d_%d] missed gg %d/%d/%d*  case=%d  z = %.1f\n", _ttd_data_.run, _ttd_data_.event, trajectory_side, cell_row, cell_layer, intersect_deltaxy_min_case, nearest_line_point_missed.z());
	    //  store 3D points...
	    // const double intersect_x_deltay = intersect_xlow ? intersect_xlow_deltay : intersect_xup_deltay;
	    // const double intersect_y_deltax = intersect_ylow ? intersect_ylow_deltax : intersect_yup_deltax;
	    // if (intersect_ylow)
	    //   deltar2_intersect_x = std::pow(,2) + std::power()

	    _track_data_.missed_cells.push_back(cell_data());
	    cell_data & _missed_cell_data_ = _track_data_.missed_cells.back();
	    _missed_cell_data_.gg_num = gg_num;
	    _missed_cell_data_.flag = 0x20; // missed-out
	    _missed_cell_data_.xfit = nearest_line_point_missed.x() - tracker_hit_x;
	    _missed_cell_data_.yfit = nearest_line_point_missed.y() - tracker_hit_y;
	    _missed_cell_data_.rfit = std::sqrt(std::pow(nearest_delta_x,2) + std::pow(nearest_delta_y,2));
	    _missed_cell_data_.zfit = nearest_line_point_missed.z();
	  }

	}
      }
    }

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


    // intercept with source plane ?
    geomtools::vector_3d src_vertex;
    if (this->intersect_with_src_plane(line_first_point, line_vector, src_vertex)) {
      _track_data_.flag |= 0x10;
      _track_data_.src_vtx[0] = src_vertex.y();
      _track_data_.src_vtx[1] = src_vertex.z();
    } else {
      _track_data_.src_vtx[0] = 0;
      _track_data_.src_vtx[1] = 0;
    }

    // intercept with main-wall plane ?
    geomtools::vector_3d mw_vertex;
    if (this->intersect_with_mw_plane(line_first_point, line_vector, trajectory_side, mw_vertex)) {
      _track_data_.flag |= 0x20;
      _track_data_.mw_vtx[0] = mw_vertex.y();
      _track_data_.mw_vtx[1] = mw_vertex.z();
    } else {
      _track_data_.mw_vtx[0] = mw_vertex.y();
      _track_data_.mw_vtx[1] = mw_vertex.z();
    }

    // intercept with x-wall 0 plane ?
    geomtools::vector_3d xw0_vertex;
    if (this->intersect_with_xw0_plane(line_first_point, line_vector, trajectory_side, xw0_vertex)) {
      _track_data_.flag |= 0x40;
      _track_data_.xw0_vtx[0] = xw0_vertex.x();
      _track_data_.xw0_vtx[1] = xw0_vertex.z();
    } else {
      _track_data_.xw0_vtx[0] = xw0_vertex.x();
      _track_data_.xw0_vtx[1] = xw0_vertex.z();
    }

    // intercept with x-wall 0 plane ?
    geomtools::vector_3d xw1_vertex;
    if (this->intersect_with_xw1_plane(line_first_point, line_vector, trajectory_side, xw1_vertex)) {
      _track_data_.flag |= 0x80;
      _track_data_.xw1_vtx[0] = xw1_vertex.x();
      _track_data_.xw1_vtx[1] = xw1_vertex.z();
    } else {
      _track_data_.xw1_vtx[0] = xw1_vertex.x();
      _track_data_.xw1_vtx[1] = xw1_vertex.z();
    }

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


  ////////////////////////////////

  const snemo::datamodel::calibrated_data & CD = event.get<snemo::datamodel::calibrated_data>("CD");

  for (const datatools::handle<snemo::datamodel::calibrated_calorimeter_hit> & calo_hit : CD.calorimeter_hits())
    {
      _ttd_data_.calos.push_back(om_data());
      om_data & _om_data_ = _ttd_data_.calos.back();

      _om_data_.om_num = snemo::datamodel::om_num(calo_hit->get_geom_id());
      _om_data_.energy = calo_hit->get_energy()/CLHEP::MeV;
      _om_data_.time = calo_hit->get_time()/CLHEP::ns;
    }

  // ttd_residual_data_print(_ttd_data_);

  _output_tree_->Fill();
  _ttd_data_.tracks.clear();
  _ttd_data_.calos.clear();

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

void ttd_residual_module::nearest_line_point_2d(double line_a, double line_b, double line_c,
						  double point_x, double point_y, double & near_point_x, double & near_point_y)
{
  near_point_x = ((line_b*point_x+line_a*point_y)/line_a - (line_c/line_b)) / (line_a/line_b + line_b/line_a);
  near_point_y = (-line_c - line_a*near_point_x)/line_b;
  return;
}

void ttd_residual_module::nearest_line_point_2d(const geomtools::vector_3d & line_first_point, const geomtools::vector_3d & line_last_point,
						const double & point_x, const double & point_y, double & near_point_x, double & near_point_y)
{
  //                        _o  <- last
  //  point   -> o       _-'
  //              ',  _-'
  //               _X'
  //            _-'       X = near point
  // first -> o'

  const double line_delta_x = line_last_point.x() - line_first_point.x();
  const double line_delta_y = line_last_point.y() - line_first_point.y();
  const double k = ((line_delta_y) * (point_x-line_first_point.x()) - (line_delta_x) * (point_y-line_first_point.y())) / (std::pow(line_delta_y,2) + std::pow(line_delta_x,2));
  near_point_x = point_x - k * (line_delta_y);
  near_point_y = point_y + k * (line_delta_x);
}

void ttd_residual_module::nearest_line_point_3d(const geomtools::vector_3d & line_first_point, const geomtools::vector_3d & line_last_point,
						const double & point_x, const double & point_y, geomtools::vector_3d & nearest_point)
{
  //                        _o  <- last
  //  point   -> o       _-'
  //              ',  _-'
  //               _X'
  //            _-'       X = nearest point
  // first -> o'

  const double line_delta_x = line_last_point.x() - line_first_point.x();
  const double line_delta_y = line_last_point.y() - line_first_point.y();
  const double line_delta_z = line_last_point.z() - line_first_point.z();

  const double k = ((line_delta_y) * (point_x-line_first_point.x()) - (line_delta_x) * (point_y-line_first_point.y())) / (std::pow(line_delta_y,2) + std::pow(line_delta_x,2));

  nearest_point.setX(point_x - k * (line_delta_y));
  nearest_point.setY(point_y + k * (line_delta_x));
  nearest_point.setZ(line_first_point.z() + line_delta_z * ((nearest_point.x()-line_first_point.x())/line_delta_x));
}

bool ttd_residual_module::intersect_with_src_plane(const geomtools::vector_3d & line_point, const geomtools::vector_3d & line_vector, geomtools::vector_3d & intersect)
{
  const geomtools::vector_3d src_point (0, 0, 0);
  const geomtools::vector_3d src_vector (1, 0, 0);

  const geomtools::vector_3d src_line_vector (line_point - src_point);
  const double src_line_factor1 = -src_vector.dot(src_line_vector);
  const double src_line_factor2 = src_vector.dot(line_vector);

  if (src_line_factor2 == 0)
    return false;

  intersect = (line_point + (src_line_factor1/src_line_factor2) * line_vector);

  if (std::fabs(intersect.z()) > 1750*CLHEP::mm)
    return false;

  if (std::fabs(intersect.y()) > 2750*CLHEP::mm)
    return false;

  return true;
}

bool ttd_residual_module::intersect_with_mw_plane(const geomtools::vector_3d & line_point, const geomtools::vector_3d & line_vector, const int & side, geomtools::vector_3d & intersect)
{
  const int x_sign = (side > 0) ? +1 : -1;

  const geomtools::vector_3d mw_point (x_sign*435*CLHEP::mm, 0, 0);
  const geomtools::vector_3d mw_vector (1, 0, 0);

  const geomtools::vector_3d mw_line_vector (line_point - mw_point);
  const double mw_line_factor1 = -mw_vector.dot(mw_line_vector);
  const double mw_line_factor2 = mw_vector.dot(line_vector);

  if (mw_line_factor2 == 0)
    return false;

  intersect = (line_point + (mw_line_factor1/mw_line_factor2) * line_vector);

  if (std::fabs(intersect.z()) > 1750*CLHEP::mm)
    return false;

  if (std::fabs(intersect.y()) > 2750*CLHEP::mm)
    return false;

  return true;
}


bool ttd_residual_module::intersect_with_xw0_plane(const geomtools::vector_3d & line_point, const geomtools::vector_3d & line_vector, const int & side, geomtools::vector_3d & intersect)
{
  const geomtools::vector_3d xw0_point (0, -2505.5*CLHEP::mm, 0);
  const geomtools::vector_3d xw0_vector (0, 1, 0);

  const geomtools::vector_3d xw0_line_vector (line_point - xw0_point);
  const double xw0_line_factor1 = -xw0_vector.dot(xw0_line_vector);
  const double xw0_line_factor2 = xw0_vector.dot(line_vector);

  if (xw0_line_factor2 == 0)
    return false;

  intersect = (line_point + (xw0_line_factor1/xw0_line_factor2) * line_vector);

  if (side > 0) {
    if ((intersect.x() < 0) || (intersect.x() > 435*CLHEP::mm))
      return false;
  } else {
    if ((intersect.x() > 0) || (intersect.x() < -435*CLHEP::mm))
      return false;
  }

  if (std::fabs(intersect.z()) > 1750*CLHEP::mm)
    return false;

  return true;
}

bool ttd_residual_module::intersect_with_xw1_plane(const geomtools::vector_3d & line_point, const geomtools::vector_3d & line_vector, const int & side, geomtools::vector_3d & intersect)
{
  const geomtools::vector_3d xw1_point (0, +2505.5*CLHEP::mm, 0);
  const geomtools::vector_3d xw1_vector (0, 1, 0);

  const geomtools::vector_3d xw1_line_vector (line_point - xw1_point);
  const double xw1_line_factor1 = -xw1_vector.dot(xw1_line_vector);
  const double xw1_line_factor2 = xw1_vector.dot(line_vector);

  if (xw1_line_factor2 == 0)
    return false;

  intersect = (line_point + (xw1_line_factor1/xw1_line_factor2) * line_vector);

  if (side > 0) {
    if ((intersect.x() < 0) || (intersect.x() > 435*CLHEP::mm))
      return false;
  } else {
    if ((intersect.x() > 0) || (intersect.x() < -435*CLHEP::mm))
      return false;
  }

  if (std::fabs(intersect.z()) > 1750*CLHEP::mm)
    return false;

  return true;
}

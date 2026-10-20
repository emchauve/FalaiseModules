#include <cstdio>
#include <stdint.h>
#include <vector>

////////////////////////////////

struct cell_data
{
  uint16_t gg_num;

  uint16_t flag;
  // 0x01 = has bottom cathode
  // 0x02 = has top cathode
  // 0x04
  // 0x08
  // 0x10 = missed-in
  // 0x20 = missed-out

  float time_anode;
  float time_bottom_cathode;
  float time_top_cathode;

  float r; // CD radius
  float z; // CD height

  float xfit;
  float yfit;
  float rfit; // FIT radius
  float zfit; // FIT height
};

////////////////////////////////

struct om_data
{
  uint16_t om_num;

  float amplitude;
  float charge;
  float energy;
  float time;

  double raw_time;

  // uint16_t flag;
  // waveform ?
};

////////////////////////////////

struct track_data
{
  // uint8_t tcd_id;
  // uint8_t ttd_id;

  uint16_t flag;
  // 0x001 = side
  // 0x002 = delayed
  // 0x004
  // 0x008
  // 0x010 = has source vertex
  // 0x020 = has MW vertex
  // 0x040 = has XW0 vertex
  // 0x080 = has XW1 vertex
  // 0x100 = has GV0 vertex
  // 0x200 = has GV1 vertex
  uint16_t nb_cells;
  // uint16_t nb_cells_withz;

  uint8_t row_min;
  uint8_t row_max;
  uint8_t layer_min;
  uint8_t layer_max;
  float z_min;
  float z_max;

  float first[3];
  float last[3];

  float src_vtx[2];
  float mw_vtx[2];
  float xw0_vtx[2];
  float xw1_vtx[2];
  float gv0_vtx[2];
  float gv1_vtx[2];

  float length;
  float theta;
  float phi;

  float chi2ndf;

  // double raw_time;
  double first_raw_time;
  double mean_raw_time;

  std::vector<cell_data> cells;
  std::vector<cell_data> missed_cells;
};

////////////////////////////////

// cluster data contains cells infos
// with all solution of fitted tracks

// struct cluster_data
// {
//   std::vector<cell_data> cells;
//   std::vector<track_data> tracks;
// };

////////////////////////////////

struct ttd_residual_data
{
  uint32_t run;
  uint32_t event;

  std::vector<track_data> tracks;
  // std::vector<cluster_data> clusters;

  std::vector<om_data> calos;
};

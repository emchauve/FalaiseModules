#include <cstdio>
#include <stdint.h>
#include <vector>

////////////////////////////////

struct cell_data
{
  uint16_t gg_num;

  uint16_t flag;
  // 0x1 = has bottom cathode
  // 0x2 = has top cathode

  float time_anode;
  float time_bottom_cathode;
  float time_top_cathode;

  float r; // CD radius
  float z; // CD height

  float rfit; // FIT radius
  float zfit; // FIT height
};

////////////////////////////////

struct om_data
{
  uint16_t om_num;

  // uint16_t flag;
};

////////////////////////////////

struct track_data
{
  // uint8_t tcd_id;
  // uint8_t ttd_id;

  uint32_t flag;
  // 0x01 = side
  // 0x02 = delayed
  // 0x04
  // 0x08
  // 0x10 = has source vertex
  // 0x20 = has MW vertex
  // 0x40 = has XW0 vertex
  // 0x80 = has XW1 vertex

  float first[3];
  float last[3];

  float src_vtx[2];
  float mw_vtx[2];
  float xw0_vtx[2];
  float xw1_vtx[2];

  float length;
  float theta;
  float phi;

  float chi2ndf;

  std::vector<cell_data> cells;
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

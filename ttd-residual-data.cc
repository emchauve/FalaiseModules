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

struct track_data
{
  // uint8_t tcd_id;
  // uint8_t ttd_id;

  uint32_t flag;
  // 0x1 = side
  // 0x2 = delayed

  float first[3];
  float last[3];

  float length;
  float theta;
  float phi;

  float chi2ndf;

  std::vector<cell_data> cells;
};

////////////////////////////////

struct ttd_residual_data
{
  uint32_t run;
  uint32_t event;

  std::vector<track_data> track;
};

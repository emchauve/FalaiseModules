#include <stdint.h>
#include <vector>

////////////////////////////////

struct cell_data
{
  uint16_t gg_num;

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
  uint32_t flag;

  float first[3];
  float last[3];

  float length;
  float theta;
  float phi;

  std::vector<cell_data> gg;
};

////////////////////////////////

struct ttd_residual_data
{
  uint32_t run;
  uint32_t event;

  std::vector<track_data> track;
};

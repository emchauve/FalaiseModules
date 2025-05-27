#include <cstdio>
#include <stdint.h>
#include <vector>

////////////////////////////////

struct cell_data
{
  uint16_t gg_num;
  uint16_t flag;

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

void ttd_residual_data_print(const ttd_residual_data & rd)
{
  printf("[%d_%d] with %zd track(s)\n", rd.run, rd.event, rd.track.size());

  for (const track_data & tr : rd.track)
    {
      printf("- nb_gg=%3zd  flag=%02d  length=%6.1f  theta=%6.1f  phi=%5.1f  chi2/ndf=%5.2f\n",
	     tr.cells.size(), tr.flag, tr.length, tr.theta, tr.phi, tr.chi2ndf);
    }
}

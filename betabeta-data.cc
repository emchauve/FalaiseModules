#include <stdint.h>
#include <vector>
struct
alpha_data
{
  // float time;
  // float length;
  double mean_anode;
  // float vtx0[3];
  // float vtx1[3];
  // uint16_t nb_cells;
};

struct beta_data
{
  float energy;
  float time;

  float source_vtx[3];
  float calo_vtx[3];

  float length;
  double mean_anode;

  uint16_t om_num;
  uint16_t nb_cells;
  uint16_t cluster_id;
  uint16_t particle_id;
};

struct gamma_data
{
  float energy;
  float time;

  uint16_t om_num;
};

struct betabeta_data
{
  uint32_t run;
  uint32_t event;

  // beta_data beta[2];
  // float source_vtx[3];
  // float deltay;
  // float deltaz;
  // float deltar;

  // uint16_t flag1, flag2;
  uint16_t cluster1, cluster2;
  uint16_t particle1, particle2;

  uint16_t om1, om2;
  float e1, e2;
  float t1, t2;
  float l1, l2;

  float calo_vtx1[3];
  float calo_vtx2[3];

  // float source_vtx1[3];
  // float source_vtx2[3];
  float source_vtx[3];
  float source_deltay;
  float source_deltaz;

  // std::vector<alpha_data> alphas;
  // std::vector<gamma_data> gammas;

};


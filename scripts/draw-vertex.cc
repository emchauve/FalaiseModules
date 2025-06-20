
TH2F src_vtx_histo ("src_vtx_histo", "SOURCE VERTEX;Y (m); Z(m)", 3000, -3, 3, 2000, -2, 2);
TH2F mw0_vtx_histo ("mw0_vtx_histo", "MW0 VERTEX;Y (m); Z(m)", 600, -3, 3, 400, -2, 2);
TH2F mw1_vtx_histo ("mw1_vtx_histo", "MW1 VERTEX;Y (m); Z(m)", 600, -3, 3, 400, -2, 2);
TH2F xw0_vtx_histo ("xw0_vtx_histo", "XW0 VERTEX;Y (m); X(m)", 200, -1, 1, 400, -2, 2);
TH2F xw1_vtx_histo ("xw1_vtx_histo", "XW1 VERTEX;Y (m); X(m)", 200, -1, 1, 400, -2, 2);

// ttd_residual->Project("src_vtx_histo", "tracks[].src_vtx[1]/1E3 : tracks[].src_vtx[0]/1E3", "(tracks.flag & 0x10)!=0");
// ttd_residual->Project("src_vtx_histo2", "tracks[].src_vtx[1]/1E3 : tracks[].src_vtx[0]/1E3", "(tracks.flag & 0x10)!=0 && @tracks.size()==1");

const float mc_source_y[6] = {-2.090, -1.254, -0.418, +0.418, +1.254, +2.090};
const float mc_source_z[7] = {-1.275, -0.847, -0.425, +0.000, +0.429, +0.847, +1.279}; // MC
// const float mc_source_y[6] = {-2.088, -1.254, -0.419, +0.417, +1.254, +2.086}; // DATA
// const float mc_source_z[7] = {-1.318, -0.887, -0.449, -0.008, +0.438, +0.881, +1.315}; // DATA

const float mc_source_sigmay = 0.100;
const float mc_source_sigmaz = 0.100;

std::vector<TH2F*> om_vtx_histos;

float interpolatex(TH1D *histo, int bin_min, int bin_max, float threshold)
{
  const float xmin = histo->GetBinCenter(bin_min);
  const float xmax = histo->GetBinCenter(bin_max);

  const float ymin = histo->GetBinContent(bin_min);
  const float ymax = histo->GetBinContent(bin_max);
  
  const float slope = (ymax-ymin)/(xmax-xmin);
  const float intercept = ymin - slope*xmin;

  return (threshold - intercept) / slope;
}

bool get_mc_source(const float y, const float z, int & source_index, float & source_distance)
{
  // a optimiser!!!
  for (int source_index_y=0; source_index_y<6; source_index_y++) {
    for (int source_index_z=0; source_index_z<7; source_index_z++) {
      const double deltay2 = std::pow((mc_source_y[source_index_y] - y)/mc_source_sigmay, 2);
      const double deltaz2 = std::pow((mc_source_z[source_index_z] - z)/mc_source_sigmaz, 2);
      // const double deltaz2 = std::pow((mc_source_z[source_index_z] + 1.573 - z)/mc_source_sigmaz, 2);
      const double deltar2 = deltay2 + deltaz2;
      if (deltar2 < 1) {
	source_index = 7*source_index_y + source_index_z;
	source_distance = std::sqrt(deltar2);
	return true;
      }
    }
  }

  return false;
}

void draw_vertex ()
{
  gStyle->SetOptTitle(1);

  src_vtx_histo.Reset();
  mw0_vtx_histo.Reset(); mw1_vtx_histo.Reset();
  xw0_vtx_histo.Reset(); xw1_vtx_histo.Reset();

  for (int om=0; om<520; om++) {
    TH2F *om_vtx_histo = new TH2F(Form("om_vtx_histo_%03d",om), "OM VERTEX;Y (m); Z(m)", 300, -3, 3, 200, -2, 2);
    om_vtx_histos.push_back(om_vtx_histo);
  }

  ttd_residual_data *ttd_data = new ttd_residual_data;
  ttd_residual->SetBranchAddress("ttd", &ttd_data);

  const int ttd_entries = ttd_residual->GetEntries();
  const int ttd_modulo_entries = ttd_entries/100;

  for (int entry=0; entry<ttd_entries; entry++) {

    if ((entry % ttd_modulo_entries) == 0) {
      printf("\r %3d %% ", entry*100/ttd_entries);
      fflush(stdout);}

    ttd_residual->GetEntry(entry);

    for (const track_data & track_d : ttd_data->tracks) {

      // if (track_d.cells.size() < 8)
      // 	continue;

      if ((track_d.flag & 0x10) == 0)
	continue;

      src_vtx_histo.Fill(track_d.src_vtx[0]/1E3, track_d.src_vtx[1]/1E3);

      if ((track_d.flag & 0x21) == 0x20)
	mw0_vtx_histo.Fill(-track_d.mw_vtx[0]/1E3, track_d.mw_vtx[1]/1E3);
    
      if ((track_d.flag & 0x21) == 0x21)
	mw1_vtx_histo.Fill(track_d.mw_vtx[0]/1E3, track_d.mw_vtx[1]/1E3);

      if ((track_d.flag & 0x40) !=0)
	xw0_vtx_histo.Fill(track_d.xw0_vtx[0]/1E3, track_d.xw0_vtx[1]/1E3);

      if ((track_d.flag & 0x80) !=0)
	xw1_vtx_histo.Fill(track_d.xw1_vtx[0]/1E3, track_d.xw1_vtx[1]/1E3);
    }

  }

  printf("\n");

  new TCanvas; src_vtx_histo.Draw("colz");
  new TCanvas; mw0_vtx_histo.Draw("colz");
  new TCanvas; mw1_vtx_histo.Draw("colz");
  TCanvas *C_xw = new TCanvas; C_xw->Divide(2,1);
  C_xw->cd(1); xw0_vtx_histo.Draw("colz");
  C_xw->cd(2); xw1_vtx_histo.Draw("colz");
}

void draw_vertex_bi ()
{
  gStyle->SetOptTitle(1);

  for (int om=0; om<520; om++) {
    TH2F *om_vtx_histo = new TH2F(Form("om_vtx_histo_%03d",om), "OM VERTEX;Y (m); Z(m)", 300, -3, 3, 200, -2, 2);
    om_vtx_histos.push_back(om_vtx_histo);
  }
      
  src_vtx_histo.Reset();
  mw0_vtx_histo.Reset(); mw1_vtx_histo.Reset();
  xw0_vtx_histo.Reset(); xw1_vtx_histo.Reset();

  ttd_residual_data *ttd_data = new ttd_residual_data;
  ttd_residual->SetBranchAddress("ttd", &ttd_data);

  const int ttd_entries = ttd_residual->GetEntries();
  const int ttd_modulo_entries = ttd_entries/100;

  for (int entry=0; entry<ttd_entries; entry++) {

    if ((entry % ttd_modulo_entries) == 0) {
      printf("\r %3d %% ", entry*100/ttd_entries);
      fflush(stdout);}
 
   ttd_residual->GetEntry(entry);

    size_t track_index = 0;

    float source_distance_min[42];
    int track_index_min[42];

    for (int source=0; source<42; source++) {
      source_distance_min[source] = std::numeric_limits<float>::max();
      track_index_min[source] = -1;
    }

    for (size_t track_index=0; track_index<ttd_data->tracks.size(); track_index++) {

      const track_data & track_d = ttd_data->tracks[track_index];

      if ((track_d.flag & 0x10) == 0)
	continue;

      int source_index; float source_distance;
      if (get_mc_source(track_d.src_vtx[0]/1E3, track_d.src_vtx[1]/1E3, source_index, source_distance)) {

	if (source_distance < source_distance_min[source_index]) {
	  source_distance_min[source_index] = source_distance;
	  track_index_min[source_index] = track_index;
	}
      }

    } // for (track_index

    for (int source=0; source<42; source++) {

      if (track_index_min[source] == -1)
	continue;

      const track_data & track_d = ttd_data->tracks[track_index_min[source]];

      src_vtx_histo.Fill(track_d.src_vtx[0]/1E3, track_d.src_vtx[1]/1E3);

      if ((track_d.flag & 0x21) == 0x20) {
	mw0_vtx_histo.Fill(-track_d.mw_vtx[0]/1E3, track_d.mw_vtx[1]/1E3);
	for (const om_data & om_d : ttd_data->calos)
	  if (om_d.om_num < 260)
	    om_vtx_histos[om_d.om_num]->Fill(-track_d.mw_vtx[0]/1E3, track_d.mw_vtx[1]/1E3);
      }

      if ((track_d.flag & 0x21) == 0x21) {
	mw1_vtx_histo.Fill(track_d.mw_vtx[0]/1E3, track_d.mw_vtx[1]/1E3);
	for (const om_data & om_d : ttd_data->calos)
	  if ((om_d.om_num >= 260) && (om_d.om_num < 520))
	    om_vtx_histos[om_d.om_num]->Fill(-track_d.mw_vtx[0]/1E3, track_d.mw_vtx[1]/1E3);
      }

      if ((track_d.flag & 0x40) !=0)
	xw0_vtx_histo.Fill(track_d.xw0_vtx[0]/1E3, track_d.xw0_vtx[1]/1E3);

      if ((track_d.flag & 0x80) !=0)
	xw1_vtx_histo.Fill(track_d.xw1_vtx[0]/1E3, track_d.xw1_vtx[1]/1E3);
    }

  }

  printf("\n");

  // new TCanvas; src_vtx_histo.Draw("colz");
  // new TCanvas; mw0_vtx_histo.Draw("colz");
  // new TCanvas; mw1_vtx_histo.Draw("colz");
  // TCanvas *C_xw = new TCanvas; C_xw->Divide(2,1);
  // C_xw->cd(1); xw0_vtx_histo.Draw("colz");
  // C_xw->cd(2); xw1_vtx_histo.Draw("colz");

  TH1F *entries_histo = new TH1F ("entries_histo", "", 500, 0, 10000);

  float it_row_sum[13], it_row_stat[13];
  float fr_row_sum[13], fr_row_stat[13];
  memset(it_row_sum, 0, sizeof(it_row_sum));
  memset(it_row_stat, 0, sizeof(it_row_stat));
  memset(fr_row_sum, 0, sizeof(fr_row_sum));
  memset(fr_row_stat, 0, sizeof(fr_row_stat));

  for (int om=0; om<520; om++) {

    entries_histo->Fill(om_vtx_histos[om]->GetEntries());
    // if (om_vtx_histos[om]->GetEntries() < 2500) printf("%d\n", om);

    if (om==0) {new TCanvas; om_vtx_histos[0]->Draw("colz");}
    if (om==260) {new TCanvas; om_vtx_histos[260]->Draw("colz");}

    if (om_vtx_histos[om]->GetEntries() < 2500) continue;

    TH1D *om_vtx_histoy = om_vtx_histos[om]->ProjectionX("om_vtx_histoy");
    const float thresholdy = om_vtx_histoy->GetMaximum()/2.0;

    int biny_min = om_vtx_histoy->GetMaximumBin();
    for (biny_min--; biny_min>1; biny_min--) {
      if (om_vtx_histoy->GetBinContent(biny_min) < thresholdy) break;
    }

    int biny_max = om_vtx_histoy->GetMaximumBin();
    for (biny_max++; biny_max<om_vtx_histoy->GetNbinsX(); biny_max++) {
      if (om_vtx_histoy->GetBinContent(biny_max) < thresholdy) break;
    }

    const float histo_ymin = interpolatex(om_vtx_histoy, biny_min, biny_min+1, thresholdy);
    const float histo_ymax = interpolatex(om_vtx_histoy, biny_max-1, biny_min, thresholdy);
    const float histo_y = 0.5 * (histo_ymin + histo_ymax);

    delete om_vtx_histoy;

    //

    TH1D *om_vtx_histoz = om_vtx_histos[om]->ProjectionY("om_vtx_histoz");
    const float thresholdz = om_vtx_histoz->GetMaximum()/2.0;

    int binz_min = om_vtx_histoy->GetMaximumBin();
    for (binz_min--; binz_min>1; binz_min--) {
      if (om_vtx_histoz->GetBinContent(binz_min) < thresholdz) break;
    }

    int binz_max = om_vtx_histoy->GetMaximumBin();
    for (binz_max++; binz_max<om_vtx_histoy->GetNbinsX(); binz_max++) {
      if (om_vtx_histoz->GetBinContent(binz_max) < thresholdz) break;
    }

    const float histo_zmin = interpolatex(om_vtx_histoy, binz_min, binz_min+1, thresholdz);
    const float histo_zmax = interpolatex(om_vtx_histoy, binz_max-1, binz_max, thresholdz);
    const float histo_z = 0.5 * (histo_zmin + histo_zmax);

    delete om_vtx_histoz;

    //

    (new TMarker (histo_y, histo_z, 20))->Draw();

    const int row = (om % 13);

    if (om < 260) {
      it_row_sum[row] += histo_z;
      it_row_stat[row]++;
    } else {
      fr_row_sum[row] += histo_z;
      fr_row_stat[row]++;
    }

  }


  for (int row=0; row<13; row++) {
    if (row == 0) continue;
    if (row == 12) continue;

    if (it_row_sum[row] != 0)
      printf("IT %d %.4f\n", row, it_row_sum[row]/it_row_stat[row]);

    if (fr_row_sum[row] != 0)
      printf("FR %d %.4f\n", row, fr_row_sum[row]/fr_row_stat[row]);

  }

}

// TH1F deltar_src_1d ("deltar_src_1d", "", 1000, 0, 1000);
TH2F deltar_src_2d ("deltar_src_2d", "", 250, -500, 500, 250, -500, 500);

void draw_vertex_bb ()
{
  gStyle->SetOptTitle(1);

  src_vtx_histo.Reset();
  mw0_vtx_histo.Reset(); mw1_vtx_histo.Reset();
  xw0_vtx_histo.Reset(); xw1_vtx_histo.Reset();

  // deltar_src_1d.Reset();
  deltar_src_2d.Reset();

  for (int om=0; om<520; om++) {
    TH2F *om_vtx_histo = new TH2F(Form("om_vtx_histo_%03d",om), "OM VERTEX;Y (m); Z(m)", 150, -3, 3, 100, -2, 2);
    om_vtx_histos.push_back(om_vtx_histo);
  }

  ttd_residual_data *ttd_data = new ttd_residual_data;
  ttd_residual->SetBranchAddress("ttd", &ttd_data);

  const int ttd_entries = ttd_residual->GetEntries();
  const int ttd_modulo_entries = ttd_entries/100;

  for (int entry=0; entry<ttd_entries; entry++) {

    if ((entry % ttd_modulo_entries) == 0) {
      printf("\r %3d %% ", entry*100/ttd_entries);
      fflush(stdout);}

    ttd_residual->GetEntry(entry);

    std::vector<size_t> selected_tracks_indexes;

    std::vector<double> deltar2_src_vtx;
    std::vector<int> closest_track_index;
    deltar2_src_vtx.reserve(ttd_data->tracks.size());
    closest_track_index.reserve(ttd_data->tracks.size());

    for (size_t ti=0; ti<ttd_data->tracks.size(); ti++) {

      const track_data & track_i = ttd_data->tracks[ti];

      if (track_i.cells.size() < 5)
	continue; // skip short track

      if ((track_i.flag & 0x10) == 0)
	continue; // skip if no vertex on source foil

      deltar2_src_vtx.push_back(std::numeric_limits<double>::max());
      closest_track_index.push_back(-1);

      for (size_t tj=ti+1; tj<ttd_data->tracks.size(); tj++) {

	const track_data & track_j = ttd_data->tracks[tj];

	if (track_i.cells[0].gg_num == track_j.cells[0].gg_num)
	  continue; // skip fit from same cluster

	if (track_i.cells.size() < 5)
	  continue; // skip short track

	if ((track_j.flag & 0x10) == 0)
	  continue; // skip if no vertex on source foil

	const float deltar2 = std::pow((track_j.src_vtx[0]-track_i.src_vtx[0]),2) + std::pow((track_j.src_vtx[1]-track_i.src_vtx[1]),2);

	if (deltar2 < deltar2_src_vtx.back()) {
	  // if (closest_track_index.back() != -1)
	  //   printf("[%d] better track (%zd,%zd)\n", entry, ti, tj);
	  deltar2_src_vtx.back() = deltar2;
	  closest_track_index.back() = tj;
	}

      } // for (tj)

      if (closest_track_index.back() >= 0) {
      // 	printf("[%d] (%zd,%d) = %.1f mm\n", entry, ti, closest_track_index.back(), std::sqrt(deltar2_src_vtx.back()));

	const track_data & track_j = ttd_data->tracks[closest_track_index.back()];

	const float deltax2 = std::pow(track_j.src_vtx[0]-track_i.src_vtx[0],2);
	const float deltay2 = std::pow(track_j.src_vtx[1]-track_i.src_vtx[1],2);

	if ((deltax2/(300.*300) + deltay2/(300.*300)) < 1) {
	  // deltar_src_1d.Fill(std::sqrt(deltar2_src_vtx.back()));
	  deltar_src_2d.Fill(track_j.src_vtx[0]-track_i.src_vtx[0], track_j.src_vtx[1]-track_i.src_vtx[1]);

	  src_vtx_histo.Fill(0.5*(track_i.src_vtx[0]+track_j.src_vtx[0])/1E3, 0.5*(track_i.src_vtx[1]+track_j.src_vtx[1])/1E3);

	  if ((track_i.flag & 0x21) == 0x20) {
	    mw0_vtx_histo.Fill(-track_i.mw_vtx[0]/1E3, track_i.mw_vtx[1]/1E3);
	    for (const om_data & om_d : ttd_data->calos)
	      if (om_d.om_num < 260)
		om_vtx_histos[om_d.om_num]->Fill(-track_i.mw_vtx[0]/1E3, track_i.mw_vtx[1]/1E3);
	  }

	  if ((track_i.flag & 0x21) == 0x21) {
	    mw1_vtx_histo.Fill(track_i.mw_vtx[0]/1E3, track_i.mw_vtx[1]/1E3);
	    for (const om_data & om_d : ttd_data->calos)
	      if ((om_d.om_num >= 260) && (om_d.om_num < 520))
		om_vtx_histos[om_d.om_num]->Fill(-track_i.mw_vtx[0]/1E3, track_i.mw_vtx[1]/1E3);
	  }

	  if ((track_i.flag & 0x40) !=0)
	    xw0_vtx_histo.Fill(track_i.xw0_vtx[0]/1E3, track_i.xw0_vtx[1]/1E3);

	  if ((track_i.flag & 0x80) !=0)
	    xw1_vtx_histo.Fill(track_i.xw1_vtx[0]/1E3, track_i.xw1_vtx[1]/1E3);

	  
	  if ((track_j.flag & 0x21) == 0x20) {
	    mw0_vtx_histo.Fill(-track_j.mw_vtx[0]/1E3, track_j.mw_vtx[1]/1E3);
	    for (const om_data & om_d : ttd_data->calos)
	      if (om_d.om_num < 260)
		om_vtx_histos[om_d.om_num]->Fill(-track_j.mw_vtx[0]/1E3, track_j.mw_vtx[1]/1E3);
	  }

	  if ((track_j.flag & 0x21) == 0x21) {
	    mw1_vtx_histo.Fill(track_j.mw_vtx[0]/1E3, track_j.mw_vtx[1]/1E3);
	    for (const om_data & om_d : ttd_data->calos)
	      if ((om_d.om_num >= 260) && (om_d.om_num < 520))
		om_vtx_histos[om_d.om_num]->Fill(-track_j.mw_vtx[0]/1E3, track_j.mw_vtx[1]/1E3);
	  }

	  if ((track_j.flag & 0x40) !=0)
	    xw0_vtx_histo.Fill(track_j.xw0_vtx[0]/1E3, track_j.xw0_vtx[1]/1E3);

	  if ((track_j.flag & 0x80) !=0)
	    xw1_vtx_histo.Fill(track_j.xw1_vtx[0]/1E3, track_j.xw1_vtx[1]/1E3);

	}

      }

    } // for (ti)

  } // for (entry)

  printf("\n");
  
  // (new TCanvas)->SetLogy();
  // deltar_src_1d.Draw();

  (new TCanvas)->SetLogz();
  deltar_src_2d.Draw("colz");

  // src_vtx_histo.Fill(track_d.src_vtx[0]/1E3, track_d.src_vtx[1]/1E3);

  // // if (!has_odd_om_it && (track_d.flag & 0x21) == 0x20)
  // if ((track_d.flag & 0x21) == 0x20)
  //   mw0_vtx_histo.Fill(-track_d.mw_vtx[0]/1E3, track_d.mw_vtx[1]/1E3);
    
  // // if (!has_odd_om_fr && (track_d.flag & 0x21) == 0x21)
  // if ((track_d.flag & 0x21) == 0x21)
  //   mw1_vtx_histo.Fill(track_d.mw_vtx[0]/1E3, track_d.mw_vtx[1]/1E3);

  // if ((track_d.flag & 0x40) !=0)
  //   xw0_vtx_histo.Fill(track_d.xw0_vtx[0]/1E3, track_d.xw0_vtx[1]/1E3);

  // if ((track_d.flag & 0x80) !=0)
  //   xw1_vtx_histo.Fill(track_d.xw1_vtx[0]/1E3, track_d.xw1_vtx[1]/1E3);

  new TCanvas; src_vtx_histo.Draw("colz");
  new TCanvas; mw0_vtx_histo.Draw("colz");
  new TCanvas; mw1_vtx_histo.Draw("colz");
  TCanvas *C_xw = new TCanvas; C_xw->Divide(2,1);
  C_xw->cd(1); xw0_vtx_histo.Draw("colz");
  C_xw->cd(2); xw1_vtx_histo.Draw("colz");

  TH1F *entries_histo = new TH1F ("entries_histo", "", 500, 0, 10000);

  float it_row_sum[13], it_row_stat[13];
  float fr_row_sum[13], fr_row_stat[13];
  memset(it_row_sum, 0, sizeof(it_row_sum));
  memset(it_row_stat, 0, sizeof(it_row_stat));
  memset(fr_row_sum, 0, sizeof(fr_row_sum));
  memset(fr_row_stat, 0, sizeof(fr_row_stat));

  for (int om=0; om<520; om++) {

    entries_histo->Fill(om_vtx_histos[om]->GetEntries());
    // if (om_vtx_histos[om]->GetEntries() < 2500) printf("%d\n", om);

    if (om==0) {new TCanvas; om_vtx_histos[0]->Draw("colz");}
    if (om==260) {new TCanvas; om_vtx_histos[260]->Draw("colz");}

    if (om_vtx_histos[om]->GetEntries() < 1500) continue;

    TH1D *om_vtx_histoy = om_vtx_histos[om]->ProjectionX("om_vtx_histoy");
    const float thresholdy = om_vtx_histoy->GetMaximum()/2.0;

    int biny_min = om_vtx_histoy->GetMaximumBin();
    for (biny_min--; biny_min>1; biny_min--) {
      if (om_vtx_histoy->GetBinContent(biny_min) < thresholdy) break;
    }

    int biny_max = om_vtx_histoy->GetMaximumBin();
    for (biny_max++; biny_max<om_vtx_histoy->GetNbinsX(); biny_max++) {
      if (om_vtx_histoy->GetBinContent(biny_max) < thresholdy) break;
    }

    const float histo_ymin = interpolatex(om_vtx_histoy, biny_min, biny_min+1, thresholdy);
    const float histo_ymax = interpolatex(om_vtx_histoy, biny_max-1, biny_min, thresholdy);
    const float histo_y = 0.5 * (histo_ymin + histo_ymax);

    delete om_vtx_histoy;

    //

    TH1D *om_vtx_histoz = om_vtx_histos[om]->ProjectionY("om_vtx_histoz");
    const float thresholdz = om_vtx_histoz->GetMaximum()/2.0;

    int binz_min = om_vtx_histoy->GetMaximumBin();
    for (binz_min--; binz_min>1; binz_min--) {
      if (om_vtx_histoz->GetBinContent(binz_min) < thresholdz) break;
    }

    int binz_max = om_vtx_histoy->GetMaximumBin();
    for (binz_max++; binz_max<om_vtx_histoy->GetNbinsX(); binz_max++) {
      if (om_vtx_histoz->GetBinContent(binz_max) < thresholdz) break;
    }

    const float histo_zmin = interpolatex(om_vtx_histoy, binz_min, binz_min+1, thresholdz);
    const float histo_zmax = interpolatex(om_vtx_histoy, binz_max-1, binz_max, thresholdz);
    const float histo_z = 0.5 * (histo_zmin + histo_zmax);

    delete om_vtx_histoz;

    //

    (new TMarker (histo_y, histo_z, 20))->Draw();

    const int row = (om % 13);

    if (om < 260) {
      it_row_sum[row] += histo_z;
      it_row_stat[row]++;
    } else {
      fr_row_sum[row] += histo_z;
      fr_row_stat[row]++;
    }
    
  }

  for (int row=0; row<13; row++) {
    if (row == 0) continue;
    if (row == 12) continue;

    if (it_row_sum[row] != 0)
      printf("IT %d %.4f\n", row, it_row_sum[row]/it_row_stat[row]);

    if (fr_row_sum[row] != 0)
      printf("FR %d %.4f\n", row, fr_row_sum[row]/fr_row_stat[row]);

  }

}

double gate1d_func (double *x, double *p)
{
  const double X = x[0];

  const double X_norm  = p[0];
  const double X_mean  = p[1];
  const double X_width = p[2]/2.0;
  const double X_sigma = p[3];
  const double X_bg    = p[4];

  const int nb_steps = 100;
  const double X_min = X - 5*X_sigma;
  const double X_max = X + 5*X_sigma;
  const double X_step = (X_max - X_min) / nb_steps;

  double X_sum = 0;

  for (int step=1; step<=nb_steps/2; step++) {

    const double XX1 = X_min + (step - 0.5) * X_step;

    if (TMath::Abs(XX1 - X_mean) < X_width)
      X_sum += X_norm * TMath::Gaus(X, XX1, X_sigma);

    const double XX2 = X_max - (step - 0.5) * X_step;

    if (TMath::Abs(XX2 - X_mean) < X_width)
      X_sum += X_norm * TMath::Gaus(X, XX2, X_sigma);

  }

  return X_bg + X_sum/nb_steps;
}

gStyle->SetOptFit(1);

// TH2F mw_vtx_om ("mw_vtx_om", ";Y (m);Z (m)", 230, -3000, 3000, 156, -2000, 2000); // ~10 bins per OM
// TH2F mw_vtx_om ("mw_vtx_om", ";Y (m);Z (m)", 469, -3000, 3000, 312, -2000, 2000); // ~20 bins per OM
// TH2F mw_vtx_om ("mw_vtx_om", ";Y (m);Z (m)", 1172, -3000, 3000, 781, -2000, 2000); // ~50 bins per OM
TH2F mw_vtx_om ("mw_vtx_om", "MW VERTEX;Y (m); Z (m)", 600, -3, 3, 400, -2, 2);
TH1F mw_vtx_om_y ("mw_vtx_om_y", "MW VERTEX;Y (m);", 1200, -3, 3);
TH1F mw_vtx_om_z ("mw_vtx_om_z", "MW VERTEX;Z (m);", 800, -2, 2);
TF1 mw_vtx_om_yfit ("mw_vtx_om_yfit", gate1d_func, -3, 3, 5);
TF1 mw_vtx_om_zfit ("mw_vtx_om_zfit", gate1d_func, -2, 2, 5);

void fit_vertex_om (const int om_num)
{

  ttd_residual_data *ttd_data = new ttd_residual_data;
  ttd_residual->SetBranchAddress("ttd", &ttd_data);

  const int ttd_entries = ttd_residual->GetEntries();
  const int mw_flag = (om_num < 260) ? 0x20 : 0x21;
  const int vtx_sgn = (om_num < 260) ? -1 : +1;
  mw_vtx_om.Reset();
  mw_vtx_om_y.Reset();
  mw_vtx_om_z.Reset();

  for (int entry=0; entry<ttd_entries; entry++) {

    ttd_residual->GetEntry(entry);

    bool has_om = false;
    for (const om_data & om_d : ttd_data->calos)
      if (om_d.om_num == om_num) {
	has_om = true; break;
      }

    if (!has_om) continue;

    for (const track_data & track_d : ttd_data->tracks) {
      if ((track_d.flag & mw_flag) != mw_flag) continue;
      if (track_d.cells.size() < 8) continue;
      mw_vtx_om.Fill(vtx_sgn*track_d.mw_vtx[0]/1E3, track_d.mw_vtx[1]/1E3);
      mw_vtx_om_y.Fill(vtx_sgn*track_d.mw_vtx[0]/1E3);
      mw_vtx_om_z.Fill(track_d.mw_vtx[1]/1E3);
    }

  } // for (entry)

  TCanvas *C0 = (TCanvas *)(gROOT->FindObject("C_draw_vertex_om_2d"));
  if (C0 == nullptr)  C0 = new TCanvas ("C_draw_vertex_om_2d", "");
  C0->cd(); mw_vtx_om.Draw("colz");

  TCanvas *C = (TCanvas *)(gROOT->FindObject("C_draw_vertex_om"));

  if (C == nullptr) {
    C = new TCanvas ("C_draw_vertex_om", "", 1400, 500);
    C->Divide(2,1);
  }

  C->cd(1);
  mw_vtx_om_y.Draw("e");
  mw_vtx_om_yfit.SetParNames("NORM", "Y0", "Ywidth", "Ysigma", "BG");
  mw_vtx_om_yfit.SetLineColor(kRed+1);
  mw_vtx_om_yfit.SetLineWidth(2);
  mw_vtx_om_yfit.SetRange(mw_vtx_om_y.GetMean()-2*0.256, mw_vtx_om_y.GetMean()+2*0.256);
  mw_vtx_om_yfit.SetParameter(0, mw_vtx_om_y.GetMaximum()*3);
  mw_vtx_om_yfit.SetParameter(1, mw_vtx_om_y.GetMean());
  mw_vtx_om_yfit.SetParameter(2, 0.375);
  mw_vtx_om_yfit.SetParameter(3, 0.025);
  mw_vtx_om_yfit.SetParameter(4, 10);
  mw_vtx_om_yfit.SetNpx(100);
  mw_vtx_om_yfit.Draw("lsame");
  mw_vtx_om_y.Fit("mw_vtx_om_yfit", "0RQ");
  mw_vtx_om_yfit.SetRange(mw_vtx_om_yfit.GetParameter(1)-1.25*0.256, mw_vtx_om_yfit.GetParameter(1)+1.25*0.256);
  mw_vtx_om_y.Fit("mw_vtx_om_yfit", "0R");

  C->cd(2);
  mw_vtx_om_z.Draw("e");
  mw_vtx_om_zfit.SetParNames("NORM", "Z0", "Zwidth", "Zsigma", "BG");
  mw_vtx_om_zfit.SetLineColor(kRed+1);
  mw_vtx_om_zfit.SetLineWidth(2);
  mw_vtx_om_zfit.SetRange(mw_vtx_om_z.GetMean()-2*0.256, mw_vtx_om_z.GetMean()+2*0.256);
  mw_vtx_om_zfit.SetParameter(0, mw_vtx_om_z.GetMaximum()*3);
  mw_vtx_om_zfit.SetParameter(1, mw_vtx_om_z.GetMean());
  mw_vtx_om_zfit.SetParameter(2, 0.275);
  mw_vtx_om_zfit.SetParameter(3, 0.025);
  mw_vtx_om_zfit.SetParameter(4, 10);
  mw_vtx_om_zfit.SetNpx(100);
  mw_vtx_om_zfit.Draw("lsame");
  mw_vtx_om_z.Fit("mw_vtx_om_zfit", "0RQ");
  mw_vtx_om_zfit.SetRange(mw_vtx_om_zfit.GetParameter(1)-1.25*0.256, mw_vtx_om_zfit.GetParameter(1)+1.25*0.256);
  mw_vtx_om_z.Fit("mw_vtx_om_zfit", "0R");

}



TH1F theta0 ("theta0", "", 1800, -180, 180);
TH1F theta1 ("theta1", "", 1800, -180, 180);

TH1F phi0 ("phi0", "", 1800, -180, 180);
TH1F phi1 ("phi1", "", 1800, -180, 180);

void draw_vertex_source (const float srcy, const float srcz, const float srcy_r=30, const float srcz_r=50)
{
  src_vtx_histo.Reset();
  mw0_vtx_histo.Reset(); mw1_vtx_histo.Reset();
  xw0_vtx_histo.Reset(); xw1_vtx_histo.Reset();

  theta0.Reset(); theta1.Reset();
  phi0.Reset(); phi1.Reset();

  ttd_residual_data *ttd_data = new ttd_residual_data;
  ttd_residual->SetBranchAddress("ttd", &ttd_data);

  const int ttd_entries = ttd_residual->GetEntries();

  for (int entry=0; entry<ttd_entries; entry++) {

    ttd_residual->GetEntry(entry);

    for (const track_data & track_d : ttd_data->tracks) {

      if ((track_d.flag & 0x10) == 0)
	continue; // no source vertex

      if ((std::pow((track_d.src_vtx[0] - srcy)/srcy_r, 2) +
	   std::pow((track_d.src_vtx[1] - srcz)/srcz_r, 2)) > 1)
	continue;

      src_vtx_histo.Fill(track_d.src_vtx[0]/1E3, track_d.src_vtx[1]/1E3);

      if ((track_d.flag & 0x01) == 0) {
	mw0_vtx_histo.Fill(-track_d.mw_vtx[0]/1E3, track_d.mw_vtx[1]/1E3);
	theta0.Fill(track_d.theta);
	phi0.Fill(track_d.phi);
      } else {
	mw1_vtx_histo.Fill(track_d.mw_vtx[0]/1E3, track_d.mw_vtx[1]/1E3);
	theta1.Fill(track_d.theta);
	phi1.Fill(track_d.phi);
      }
    }

  } // for (entry)
    
  new TCanvas;
  src_vtx_histo.Draw("colz");

  new TCanvas;
  mw0_vtx_histo.Draw("colz");

  new TCanvas;
  mw1_vtx_histo.Draw("colz");

}

  

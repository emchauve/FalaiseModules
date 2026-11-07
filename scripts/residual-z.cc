
double residual_sigma_func (double *x, double *p)
{
  // [0]
  // [0]+[1] = max
  // [0] = min
  // min+(max-min*(1-*std::sqrt(1-std::pow(x/[1],2))

  const double X = x[0];

  const double XRANGE = p[0];
  const double SIGMA0 = p[1];
  const double SIGMA_LIMIT = p[2];

  const double XRATIO = X/XRANGE;
  double SIGMA = 0;

  if (std::fabs(XRATIO) < 1)
    SIGMA = SIGMA0 * std::sqrt(1-std::pow(X/XRANGE,2));

  return std::max(SIGMA, SIGMA_LIMIT);
}

double residual_bottom_sigma_func (double *x, double *p)
{
  const double X = x[0];

  const double XRANGE = p[0];
  const double SIGMA_A = p[1];
  const double SIGMA_B = p[2];
  const double SIGMA_LIMIT = p[3];

  const double XRATIO = X/XRANGE;
  double SIGMA = 0;

  if (std::fabs(XRATIO) < 1)
    SIGMA = std::sqrt(std::pow(SIGMA_A,2)*std::pow(1+XRATIO,2) + std::pow(SIGMA_B,2)*(1+XRATIO));

  return std::max(SIGMA, SIGMA_LIMIT);
}

double residual_top_sigma_func (double *x, double *p)
{
  const double X = x[0];

  const double XRANGE = p[0];
  const double SIGMA_A = p[1];
  const double SIGMA_B = p[2];
  const double SIGMA_LIMIT = p[3];

  const double XRATIO = X/XRANGE;
  double SIGMA = 0;

  if (std::fabs(XRATIO) < 1)
    SIGMA = std::sqrt(std::pow(SIGMA_A,2)*std::pow(1-XRATIO,2) + std::pow(SIGMA_B,2)*(1-XRATIO));

  return std::max(SIGMA, SIGMA_LIMIT);
}


// TH2F ppt_with_both_cathodes ("ppt_with_both_cathodes", ";Z (mm);PPT (us)", 300, -1500, 1500, 200, 0, 100);
// TH2F ppt_with_bottom_cathode ("ppt_with_bottom_cathode", ";Z (mm);PPT (us)", 300, -1500, 1500, 200, 0, 100);
// TH2F ppt_with_top_cathode ("ppt_with_top_cathode", ";Z (mm);PPT (us)", 300, -1500, 1500, 200, 0, 100);
// TH2F ppt_all_cases ("ppt_all_cases", ";Z (mm);PPT (us)", 300, -1500, 1500, 200, 0, 100);
TH2F ppt_with_both_cathodes[2034];
TH2F ppt_with_bottom_cathode[2034];
TH2F ppt_with_top_cathode[2034];
TH2F ppt_all_cases[2034];
TH1F ppt_mean_1d ("ppt_mean_1d", "", 200, 0, 100);

TLine znorm_low_line (-1380*0.9, 0, -1380*0.9, 100);
TLine znorm_up_line (1380*0.9, 0, 1380*0.9, 100);

void init_ppt_histos()
{
  for (int gg=0; gg<2034; gg++) {
    ppt_with_both_cathodes[gg] = TH2F(Form("ppt_with_both_cathodes_%04d",gg), ";Z (mm);PPT (us)", 300, -1500, 1500, 200, 0, 100);
    ppt_with_bottom_cathode[gg] = TH2F(Form("ppt_with_bottom_cathode_%04d",gg), ";Z (mm);PPT (us)", 300, -1500, 1500, 200, 0, 100);
    ppt_with_top_cathode[gg] = TH2F(Form("ppt_with_top_cathode_%04d",gg), ";Z (mm);PPT (us)", 300, -1500, 1500, 200, 0, 100);
    ppt_all_cases[gg] = TH2F(Form("ppt_all_cases_%04d",gg), ";Z (mm);PPT (us)", 300, -1500, 1500, 200, 0, 100);
  }

}

void ppt_measurement(const int gg_num_to_draw)
{
  init_ppt_histos();
  // ppt_with_both_cathodes.Reset();
  // ppt_with_bottom_cathode.Reset();
  // ppt_with_top_cathode.Reset();
  // ppt_all_cases.Reset();
  ppt_mean_1d.Reset();

  TFile *ttd_residual_file = new TFile ("snemo_run-2683_ttd-trackfit_residual.root", "READ");
  TTree *ttd_residual = (TTree*)(ttd_residual_file->Get("ttd_residual"));

  ttd_residual_data *ttd_data = new ttd_residual_data;
  ttd_residual->SetBranchAddress("ttd", &ttd_data);

  const long ttd_entries = ttd_residual->GetEntries();
  const long ttd_modulo_entries = ttd_entries/100 - 1;

  for (long entry=0; entry<ttd_entries; entry++) {

    if ((entry % ttd_modulo_entries) == 0) {
      printf("\r %3ld %% ", entry*100/ttd_entries);
      fflush(stdout);}

    ttd_residual->GetEntry(entry);

    for (const track_data & track_d : ttd_data->tracks) {

      // track with at least 5 cells
      if (track_d.nb_cells < 5) continue;

      // skip "delayed" tracks
      if (track_d.flag & 0x2) continue;

      // skip track of single row width (left/right ambiguity cases)
      if (track_d.row_min == track_d.row_max) continue;

      // skip track of single layer width (left/right ambiguity cases)
      if (track_d.layer_min == track_d.layer_max) continue;

      // // skip horizontal tracks
      // if (std::fabs(track_d.theta) < 22.5) continue;

      // // skip bad fit
      // if (track_d.chi2ndf > 10) continue;

      for (const cell_data & cell_d : track_d.cells) {

	// const float z_residual = cell_d.z - cell_d.zfit;
	const float z_norm = cell_d.zfit/1380.;

	float cell_ppt = 0;
        switch (cell_d.flag & 0x3) {

	case 0x1: // bottom cathode only
	  cell_ppt = 2*cell_d.time_bottom_cathode / (1 + z_norm);
	  ppt_with_bottom_cathode[cell_d.gg_num].Fill(cell_d.zfit, cell_ppt/1E3);
	  ppt_all_cases[cell_d.gg_num].Fill(cell_d.zfit, cell_ppt/1E3);
	  break;

	case 0x2: // top cathode only
	  cell_ppt = 2*cell_d.time_top_cathode / (1 - z_norm);
	  ppt_with_top_cathode[cell_d.gg_num].Fill(cell_d.zfit, cell_ppt/1E3);
	  ppt_all_cases[cell_d.gg_num].Fill(cell_d.zfit, cell_ppt/1E3);
	  break;

	case 0x3: // bottom + top cathodes
	  cell_ppt = cell_d.time_bottom_cathode + cell_d.time_top_cathode;
	  ppt_with_both_cathodes[cell_d.gg_num].Fill(cell_d.zfit, cell_ppt/1E3);
	  ppt_all_cases[cell_d.gg_num].Fill(cell_d.zfit, cell_ppt/1E3);
	  break;

	} // switch
      
      } // for (cell_data)

    } // for (track_data)
      
  } // for (entry)

  printf("\n");

  std::ofstream ppt_file ("ppt_file.txt");

  for (int gg=0; gg<2034; gg++) {
    ppt_with_both_cathodes[gg].GetXaxis()->SetRangeUser(-1380*0.9, +1380*0.9);
    const float mean_ppt_with_both_cathodes = ppt_with_both_cathodes[gg].GetMean(2);
    const float mean_ppt_with_both_cathodes_error = ppt_with_both_cathodes[gg].GetMeanError(2);
    const float mean_ppt_with_both_cathodes_rms = ppt_with_both_cathodes[gg].GetRMS(2);
    ppt_with_both_cathodes[gg].GetXaxis()->SetRangeUser(-1500, 1500);
    // printf("<ppt> B+T = %4.1f +- %4.1f (ERR) +- %4.1f (RMS)\n", mean_ppt_with_both_cathodes, mean_ppt_with_both_cathodes_error, mean_ppt_with_both_cathodes_rms);

    ppt_with_bottom_cathode[gg].GetXaxis()->SetRangeUser(-1380*0.9, +1380*0.9);
    const float mean_ppt_with_bottom_cathode = ppt_with_bottom_cathode[gg].GetMean(2);
    const float mean_ppt_with_bottom_cathode_error = ppt_with_bottom_cathode[gg].GetMeanError(2);
    const float mean_ppt_with_bottom_cathode_rms = ppt_with_bottom_cathode[gg].GetRMS(2);
    ppt_with_bottom_cathode[gg].GetXaxis()->SetRangeUser(-1500, 1500);
    // printf("<ppt>  B  = %4.1f +- %4.1f (ERR) +- %4.1f (RMS)\n", mean_ppt_with_bottom_cathode, mean_ppt_with_bottom_cathode_error, mean_ppt_with_bottom_cathode_rms);

    ppt_with_top_cathode[gg].GetXaxis()->SetRangeUser(-1380*0.9, +1380*0.9);
    const float mean_ppt_with_top_cathode = ppt_with_top_cathode[gg].GetMean(2);
    const float mean_ppt_with_top_cathode_error = ppt_with_top_cathode[gg].GetMeanError(2);
    const float mean_ppt_with_top_cathode_rms = ppt_with_top_cathode[gg].GetRMS(2);
    ppt_with_top_cathode[gg].GetXaxis()->SetRangeUser(-1500, 1500);
    // printf("<ppt>  T  = %4.1f +- %4.1f (ERR) +- %4.1f (RMS)\n", mean_ppt_with_top_cathode, mean_ppt_with_top_cathode_error, mean_ppt_with_top_cathode_rms);

    ppt_all_cases[gg].GetXaxis()->SetRangeUser(-1380*0.9, +1380*0.9);
    const float mean_ppt_all_cases = ppt_all_cases[gg].GetMean(2);
    const float mean_ppt_all_cases_error = ppt_all_cases[gg].GetMeanError(2);
    const float mean_ppt_all_cases_rms = ppt_all_cases[gg].GetRMS(2);
    ppt_all_cases[gg].GetXaxis()->SetRangeUser(-1500, 1500);
    // printf("<ppt> ALL = %4.1f +- %4.1f (ERR) +- %4.1f (RMS)\n", mean_ppt_all_cases, mean_ppt_all_cases_error, mean_ppt_all_cases_rms);

    ppt_mean_1d.Fill(mean_ppt_all_cases);

    if (mean_ppt_all_cases > 40)
      ppt_file << gg << "\t" << mean_ppt_all_cases << "\t" << mean_ppt_all_cases_error << "\t" << mean_ppt_all_cases_rms << std::endl;

  }

  ppt_file.close();

  //

  TCanvas *ppt_canvas = (TCanvas*)(gROOT->FindObject("ppt_canvas"));
  if (ppt_canvas == nullptr)
    ppt_canvas = new TCanvas ("ppt_canvas", "");
  else ppt_canvas->cd();

  ppt_with_both_cathodes[gg_num_to_draw].SetMarkerColor(kRed+1);
  ppt_with_both_cathodes[gg_num_to_draw].SetMarkerStyle(21);
  ppt_with_both_cathodes[gg_num_to_draw].SetMarkerSize(0.5);
  ppt_with_both_cathodes[gg_num_to_draw].Draw("scat");

  ppt_with_bottom_cathode[gg_num_to_draw].SetMarkerColor(kBlue);
  ppt_with_bottom_cathode[gg_num_to_draw].SetMarkerStyle(23);
  ppt_with_bottom_cathode[gg_num_to_draw].SetMarkerSize(0.5);
  ppt_with_bottom_cathode[gg_num_to_draw].Draw("scatsame");

  ppt_with_top_cathode[gg_num_to_draw].SetMarkerColor(kGreen+2);
  ppt_with_top_cathode[gg_num_to_draw].SetMarkerStyle(22);
  ppt_with_top_cathode[gg_num_to_draw].SetMarkerSize(0.5);
  ppt_with_top_cathode[gg_num_to_draw].Draw("scatsame");

  znorm_low_line.SetLineStyle(kDotted);
  znorm_low_line.Draw();

  znorm_up_line.SetLineStyle(kDotted);
  znorm_up_line.Draw();

}


std::vector<float> zmean;

std::vector<TH1F*> residual_histos;
std::vector<TH1F*> residual_bottom_histos;
std::vector<TH1F*> residual_top_histos;

TH2F residual_histo_2d ("residual_histo_2d", "", 300, -1500, 1500, 200, -100, 100);
TH2F residual_bottom_histo_2d ("residual_bottom_histo_2d", "", 300, -1500, 1500, 200, -100, 100);
TH2F residual_top_histo_2d ("residual_top_histo_2d", "", 300, -1500, 1500, 200, -100, 100);

void residual_z()
{
  float ppt_table[2034];
  memset(ppt_table, 0, sizeof(ppt_table));

  {
    std::ifstream ppt_file ("data/snemo_run-2683_ppt.txt");

    int gg; float ppt;
    while (ppt_file >> gg >> ppt)
      ppt_table[gg] = ppt*1E3;

    ppt_file.close();
  }

  const float zstep = 25;

  for (float z=-1350; z<1351; z+=zstep) {
    zmean.push_back(z);
    residual_histos.push_back(new TH1F(Form("residual_histo_%02zu",residual_histos.size()), "", 200, -100, 100));
    residual_bottom_histos.push_back(new TH1F(Form("residual_bottom_histo_%02zu",residual_histos.size()), "", 100, -200, 200));
    residual_top_histos.push_back(new TH1F(Form("residual_top_histo_%02zu",residual_histos.size()), "", 100, -200, 200));
  }

  TFile *ttd_residual_file = new TFile ("snemo_run-2683_ttd-trackfit_residual.root", "READ");
  TTree *ttd_residual = (TTree*)(ttd_residual_file->Get("ttd_residual"));

  ttd_residual_data *ttd_data = new ttd_residual_data;
  ttd_residual->SetBranchAddress("ttd", &ttd_data);

  const long ttd_entries = ttd_residual->GetEntries();
  const long ttd_modulo_entries = ttd_entries/100 - 1;

  for (long entry=0; entry<ttd_entries; entry++) {

    if ((entry % ttd_modulo_entries) == 0) {
      printf("\r %3ld %% ", entry*100/ttd_entries);
      fflush(stdout);}

    ttd_residual->GetEntry(entry);

    for (const track_data & track_d : ttd_data->tracks) {

      // track with at least 5 cells
      if (track_d.nb_cells < 5) continue;

      // skip "delayed" tracks
      if (track_d.flag & 0x2) continue;

      // skip track of single row width (left/right ambiguity cases)
      if (track_d.row_min == track_d.row_max) continue;

      // skip track of single layer width (left/right ambiguity cases)
      if (track_d.layer_min == track_d.layer_max) continue;

      int nb_cells_with_full_ppt = 0;

      for (const cell_data & cell_d : track_d.cells) {

	// const float z_norm = cell_d.zfit/1380.;
	float cell_znorm=0, cell_ppt=0, cell_z = 0, z_residual = 0;

        switch (cell_d.flag & 0x3) {
	  
	case 0x1: // bottom cathode
	  if (ppt_table[cell_d.gg_num] == 0) continue;
	  cell_ppt = ppt_table[cell_d.gg_num];
	  cell_znorm = (2*cell_d.time_bottom_cathode - cell_ppt) / cell_ppt;
	  cell_z = 1380 * (cell_znorm + 0.005 * std::pow(1+cell_znorm,2));
	  // cell_z = 1380 * (2*cell_d.time_bottom_cathode - ppt_table[cell_d.gg_num]) / ppt_table[cell_d.gg_num];
	  z_residual = cell_z - cell_d.zfit;
	  for (size_t i=0; i<zmean.size(); i++)
	    if (std::fabs(cell_d.zfit-zmean[i]) < zstep/2) {
	      residual_bottom_histos[i]->Fill(z_residual);
	      break;}
	  residual_bottom_histo_2d.Fill(cell_d.zfit, z_residual);
	  break;

	case 0x2: // top cathode
	  if (ppt_table[cell_d.gg_num] == 0) continue;
	  cell_ppt = ppt_table[cell_d.gg_num];
	  cell_znorm = (cell_ppt - 2*cell_d.time_top_cathode) / cell_ppt;
	  cell_z = 1380 * (cell_znorm - 0.005 * std::pow(1-cell_znorm,2));
	  z_residual = cell_z - cell_d.zfit;
	  for (size_t i=0; i<zmean.size(); i++)
	    if (std::fabs(cell_d.zfit-zmean[i]) < zstep/2) {
	      residual_top_histos[i]->Fill(z_residual);
	      break;}
	  residual_top_histo_2d.Fill(cell_d.zfit, z_residual);
	  break;
	  
	case 0x3: // bottom + top cathodes
	  cell_z = cell_d.z;
	  z_residual = cell_d.z - cell_d.zfit;
	  for (size_t i=0; i<zmean.size(); i++)
	    if (std::fabs(cell_d.zfit-zmean[i]) < zstep/2) {
	      residual_histos[i]->Fill(z_residual);
	      break;}
	  residual_histo_2d.Fill(cell_d.zfit, z_residual);
	  nb_cells_with_full_ppt++;
	  break;

	default:
	  continue;

	} // switch (cell_flag)

	// printf("[%ld] %04d-%d = %7.1f mm\n", entry, cell_d.gg_num, (cell_d.flag&0x3), cell_z);

      } // for (cell_data)

      // if ((nb_cells_with_full_ppt < (track_d.nb_cells - nb_cells_with_full_ppt)) &&
      // 	  (nb_cells_with_full_ppt == 2) && (track_d.nb_cells > 10))
      // 	printf("[%ld] track with %d/%d full ppt\n", entry, nb_cells_with_full_ppt, track_d.nb_cells);

    } // for (track_data)
      
  } // for (entry)

  printf("\n");
  // return;

  TGraphErrors *residual_mean_graph = new TGraphErrors;
  TGraphErrors *residual_sigma_graph = new TGraphErrors;

  TGraphErrors *residual_bottom_mean_graph = new TGraphErrors;
  TGraphErrors *residual_bottom_sigma_graph = new TGraphErrors;

  TGraphErrors *residual_top_mean_graph = new TGraphErrors;
  TGraphErrors *residual_top_sigma_graph = new TGraphErrors;
  
  for (size_t i=0; i<zmean.size(); i++) {

    TF1 residual_fit ("residual_fit", "gaus", -100, 100);
    residual_histos[i]->Fit("residual_fit", "0QR");

    const float sigma_range = 1.5;

    //

    float residual_mean = residual_fit.GetParameter(1);
    float residual_sigma = residual_fit.GetParameter(2);
    residual_fit.SetRange(residual_mean-sigma_range*residual_sigma, residual_mean+sigma_range*residual_sigma);
    residual_histos[i]->Fit("residual_fit", "0QR");

    residual_mean = residual_fit.GetParameter(1);
    residual_sigma = residual_fit.GetParameter(2);
    residual_fit.SetRange(residual_mean-sigma_range*residual_sigma, residual_mean+sigma_range*residual_sigma);
    residual_histos[i]->Fit("residual_fit", "0QR");

    residual_mean = residual_fit.GetParameter(1);
    residual_sigma = residual_fit.GetParameter(2);

    const int npoint = residual_mean_graph->GetN();
    residual_mean_graph->SetPoint(npoint, zmean[i], residual_mean);
    residual_mean_graph->SetPointError(npoint, 0, residual_fit.GetParError(1));

    residual_sigma_graph->SetPoint(npoint, zmean[i], residual_sigma);
    residual_sigma_graph->SetPointError(npoint, 0, residual_fit.GetParError(2));

    //

    residual_mean = residual_fit.GetParameter(1);
    residual_sigma = residual_fit.GetParameter(2);
    residual_fit.SetRange(residual_mean-sigma_range*residual_sigma, residual_mean+sigma_range*residual_sigma);
    residual_bottom_histos[i]->Fit("residual_fit", "0QR");

    residual_mean = residual_fit.GetParameter(1);
    residual_sigma = residual_fit.GetParameter(2);
    residual_fit.SetRange(residual_mean-sigma_range*residual_sigma, residual_mean+sigma_range*residual_sigma);
    residual_bottom_histos[i]->Fit("residual_fit", "0QR");

    residual_mean = residual_fit.GetParameter(1);
    residual_sigma = residual_fit.GetParameter(2);

    residual_bottom_mean_graph->SetPoint(npoint, zmean[i], residual_mean);
    residual_bottom_mean_graph->SetPointError(npoint, 0, residual_fit.GetParError(1));

    residual_bottom_sigma_graph->SetPoint(npoint, zmean[i], residual_sigma);
    residual_bottom_sigma_graph->SetPointError(npoint, 0, residual_fit.GetParError(2));

    //

    residual_mean = residual_fit.GetParameter(1);
    residual_sigma = residual_fit.GetParameter(2);
    residual_fit.SetRange(residual_mean-sigma_range*residual_sigma, residual_mean+sigma_range*residual_sigma);
    residual_top_histos[i]->Fit("residual_fit", "0QR");

    residual_mean = residual_fit.GetParameter(1);
    residual_sigma = residual_fit.GetParameter(2);
    residual_fit.SetRange(residual_mean-sigma_range*residual_sigma, residual_mean+sigma_range*residual_sigma);
    residual_top_histos[i]->Fit("residual_fit", "0QR");

    residual_mean = residual_fit.GetParameter(1);
    residual_sigma = residual_fit.GetParameter(2);

    residual_top_mean_graph->SetPoint(npoint, zmean[i], residual_mean);
    residual_top_mean_graph->SetPointError(npoint, 0, residual_fit.GetParError(1));

    residual_top_sigma_graph->SetPoint(npoint, zmean[i], residual_sigma);
    residual_top_sigma_graph->SetPointError(npoint, 0, residual_fit.GetParError(2));

  }

  // (new TCanvas)->SetLogy();
  // residual_histos[6]->SetLineColor(kBlue);
  // residual_histos[6]->Draw("e");
  // residual_histos[4]->SetLineColor(kBlue+1);
  // residual_histos[4]->Draw("esame");
  // residual_histos[2]->SetLineColor(kBlue+2);
  // residual_histos[2]->Draw("esame");
  // residual_histos[0]->SetLineColor(kBlue+3);
  // residual_histos[0]->Draw("esame");
  
  TCanvas *residual_mean_canvas = new TCanvas;
  residual_mean_canvas->SetGridx();
  residual_mean_canvas->SetGridy();
  residual_mean_graph->SetTitle(";Z (mm);Z residual (mm)");
  residual_mean_graph->GetXaxis()->SetNdivisions(508);
  residual_mean_graph->GetYaxis()->SetNdivisions(508);
  residual_mean_graph->GetYaxis()->SetRangeUser(-2.5, 2.5);
  residual_mean_graph->Draw("apl");
  residual_bottom_mean_graph->SetLineColor(kBlue);
  residual_bottom_mean_graph->SetMarkerColor(kBlue);
  residual_bottom_mean_graph->Draw("plsame");
  residual_top_mean_graph->SetLineColor(kGreen+2);
  residual_top_mean_graph->SetMarkerColor(kGreen+2);
  residual_top_mean_graph->Draw("plsame");

  TCanvas *residual_sigma_canvas = new TCanvas;
  residual_sigma_canvas->SetGridx();
  residual_sigma_canvas->SetGridy();
  residual_sigma_graph->SetTitle(";Z (mm);Z error (mm)");
  residual_sigma_graph->GetXaxis()->SetNdivisions(508);
  residual_sigma_graph->GetYaxis()->SetNdivisions(508);
  residual_sigma_graph->GetYaxis()->SetRangeUser(0, 50);
  residual_sigma_graph->Draw("apl");
  residual_bottom_sigma_graph->SetLineColor(kBlue);
  residual_bottom_sigma_graph->SetMarkerColor(kBlue);
  residual_bottom_sigma_graph->Draw("plsame");
  residual_top_sigma_graph->SetLineColor(kGreen+2);
  residual_top_sigma_graph->SetMarkerColor(kGreen+2);
  residual_top_sigma_graph->Draw("plsame");

  // TF1 *residual_sigma_fit = new TF1 ("residual_sigma_fit", "[0]*std::sqrt(1-std::pow(x/[1],2))", -1350, 1350);
  TF1 *residual_sigma_fit = new TF1 ("residual_sigma_fit", residual_sigma_func, -1400, 1400, 3);
  residual_sigma_fit->SetParNames("H0", "#sigma_{0}", "#sigma_{lim}");
  residual_sigma_fit->SetParameters(1380, 12.5, 54.2);
  residual_sigma_fit->FixParameter(0, 1380);
  residual_sigma_fit->SetLineColor(kRed+1);
  residual_sigma_fit->SetLineWidth(2);
  residual_sigma_fit->SetNpx(100);
  residual_sigma_graph->Fit("residual_sigma_fit", "0", "", -1400, 1400);
  residual_sigma_fit->Draw("lsame");

  TF1 *residual_bottom_sigma_fit = new TF1 ("residual_bottom_sigma_fit", residual_bottom_sigma_func, -1400, 1400, 4);
  residual_bottom_sigma_fit->SetParNames("H0", "#sigma_{a}", "#sigma_{b}", "#sigma_{lim}");
  residual_bottom_sigma_fit->SetParameters(1380, 12.5, 30.0, 4.2, 5.0);
  residual_bottom_sigma_fit->FixParameter(0, 1380);
  residual_bottom_sigma_fit->FixParameter(1, residual_sigma_fit->GetParameter(1));
  residual_bottom_sigma_fit->FixParameter(3, residual_sigma_fit->GetParameter(2));
  residual_bottom_sigma_fit->SetLineColor(kBlue);
  residual_bottom_sigma_fit->SetLineWidth(2);
  residual_bottom_sigma_fit->SetNpx(100);
  residual_bottom_sigma_graph->Fit("residual_bottom_sigma_fit", "0", "", -1400, 1400);
  residual_bottom_sigma_fit->Draw("lsame");

  TF1 *residual_top_sigma_fit = new TF1 ("residual_top_sigma_fit", residual_top_sigma_func, -1400, 1400, 4);
  residual_top_sigma_fit->SetParNames("H0", "#sigma_{a}", "#sigma_{b}", "#sigma_{lim}");
  residual_top_sigma_fit->SetParameters(1380, 12.5, 30.0, 4.2, 5.0);
  residual_top_sigma_fit->FixParameter(0, 1380);
  residual_top_sigma_fit->FixParameter(1, residual_sigma_fit->GetParameter(1));
  residual_top_sigma_fit->FixParameter(3, residual_sigma_fit->GetParameter(2));
  residual_top_sigma_fit->SetLineColor(kGreen);
  residual_top_sigma_fit->SetLineWidth(2);
  residual_top_sigma_fit->SetNpx(100);
  residual_top_sigma_graph->Fit("residual_top_sigma_fit", "0", "", -1400, 1400);
  residual_top_sigma_fit->Draw("lsame");

  //                -1.5/+1.5  -1.0/+1.0
  // chi2ndf <20 : 12.6 - 4.7  11.4 - 4.4
  // chi2ndf <10 : 11.9 - 4.4  10.7 - 4.2
  // chi2ndf  <5 : 11.2 - 4.0   9.8 - 3.7
}

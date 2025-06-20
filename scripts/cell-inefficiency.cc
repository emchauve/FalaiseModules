// efficiency_cell_map_polar.Draw
// [0] * (1 + TMath::Erf( (x-[1]) / [2] )
// 0.5, 23, 3

TEllipse cell_ellipse (0, 0, 22, 22);
TBox cell_box (-22, -22, 22, 22);

TH2F cell_map ("cell_map", "Cell Hit;X (mm); Y (mm)", 125, -25, 25, 125, -25, 25);
TH2F missed_cell_map ("missed_cell_map", "missed Cell Hits;X (mm); Y (mm)", 125, -25, 25, 125, -25, 25);
TH2F total_cell_map ("total_cell_map", ";X (mm); Y (mm)", 125, -25, 25, 125, -25, 25);
TH2F efficiency_cell_map ("efficiency_cell_map", "Cell inefficiency;X (mm); Y (mm)", 125, -25, 25, 125, -25, 25);

TH2F cell_map_polar ("cell_map_polar", "Cell Hit;#Theta (degrees); R (mm)", 360/10, -180, 180, 350/10, 0, 35);
TH2F missed_cell_map_polar ("missed_cell_map_polar", "missed ;#Theta (degrees); R (mm)", 360/10, -180, 180, 350/10, 0, 35);
TH2F total_cell_map_polar ("total_cell_map_polar", ";#Theta (degrees); R (mm)", 360/10, -180, 180, 350/10, 0, 35);
TH2F efficiency_cell_map_polar ("efficiency_cell_map_polar", "Cell inefficiency ;#Theta (degrees); R (mm)", 360/10, -180, 180, 350/10, 0, 35);

TH2F cell_histo ("cell_histo", "", 113+2, 0-1, 113+1, 18+3, 0, 18+3);
TH2F missed_cell_histo ("missed_cell_histo", "", 113+2, 0-1, 113+1, 18+3, 0, 18+3);
TH2F total_cell_histo ("total_cell_histo", "", 113+2, 0-1, 113+1, 18+3, 0, 18+3);
TH2F efficiency_cell_histo ("efficiency_cell_histo", "Cell inefficiency", 113+2, 0-1, 113+1, 18+3, 0, 18+3);

std::vector<uint16_t> dead_cells = {10, 19, 24, 27, 39, 81, 99, 189, 223, 507, 508, 509, 513, 575, 630, 675, 676, 677, 678, 679, 680, 681, 682, 683, 701, 756, 782, 802, 912, 1305, 1448, 1525, 1526, 1616, 1718, 1739, 1778, 1916, 1921}; // 1718 after run XXX

void cell_inefficiency()
{
  gStyle->SetOptTitle(1);

  cell_box.SetFillStyle(0);
  cell_box.SetLineStyle(kSolid);
  cell_box.SetLineWidth(5);

  cell_ellipse.SetFillStyle(0);
  cell_ellipse.SetLineStyle(kDotted);
  cell_ellipse.SetLineWidth(5);

  bool is_dead[2034];
  memset(is_dead, 0, sizeof(is_dead));
  for (const uint16_t dead_cell : dead_cells)
    is_dead[dead_cell] = true;
  
  ttd_residual_data *ttd_data = new ttd_residual_data;
  ttd_residual->SetBranchAddress("ttd", &ttd_data);

  const int ttd_entries = ttd_residual->GetEntries();

  for (int entry=0; entry<ttd_entries; entry++) {

    ttd_residual->GetEntry(entry);

    for (const track_data & track_d : ttd_data->tracks) {

      if (track_d.cells.size() < 5)
	continue;

      if (track_d.chi2ndf > 2)
	continue;

      if ((track_d.flag & 0x02) != 0)
	continue; // not delayed

      for (const cell_data & cell_d : track_d.cells) {

	if (is_dead[cell_d.gg_num])
	  continue;

	// if (cell_d.gg_num != 936)
	//   continue;

	cell_map.Fill(cell_d.xfit, cell_d.yfit);

	const float cell_d_r_fit = std::sqrt(std::pow(cell_d.xfit,2)+std::pow(cell_d.yfit,2));
	const float cell_d_theta_fit = std::atan2(cell_d.yfit,cell_d.xfit) * (180.0/TMath::Pi());
	cell_map_polar.Fill(cell_d_theta_fit, cell_d_r_fit);

	const uint16_t gg_side = (cell_d.gg_num < 1017) ? 0 : 1;
	const uint16_t gg_row = (cell_d.gg_num - gg_side*1017) / 9;
	const uint16_t gg_layer = cell_d.gg_num % 9;
	if (gg_side == 0) cell_histo.Fill(0.5+gg_row, 11.5+gg_layer);
	else cell_histo.Fill(0.5+gg_row, 9.5-gg_layer);	
      }

      for (const cell_data & cell_d : track_d.missed_cells) {

	if (is_dead[cell_d.gg_num])
	  continue;

	// if (cell_d.gg_num != 936)
	//   continue;

	missed_cell_map.Fill(cell_d.xfit, cell_d.yfit);

	const float cell_d_r_fit = std::sqrt(std::pow(cell_d.xfit,2)+std::pow(cell_d.yfit,2));
	const float cell_d_theta_fit = std::atan2(cell_d.yfit,cell_d.xfit) * (180.0/TMath::Pi());
	missed_cell_map_polar.Fill(cell_d_theta_fit, cell_d_r_fit);

	const uint16_t gg_side = (cell_d.gg_num < 1017) ? 0 : 1;
	const uint16_t gg_row = (cell_d.gg_num - gg_side*1017) / 9;
	const uint16_t gg_layer = cell_d.gg_num % 9;
	if (gg_side == 0) missed_cell_histo.Fill(0.5+gg_row, 11.5+gg_layer);
	else missed_cell_histo.Fill(0.5+gg_row, 9.5-gg_layer);	
      }
    } // for (track_d)

  } // for (entry)

  total_cell_map.Add(&cell_map, &missed_cell_map);
  efficiency_cell_map.Divide(&missed_cell_map, &total_cell_map);

  total_cell_map_polar.Add(&cell_map_polar, &missed_cell_map_polar);
  efficiency_cell_map_polar.Divide(&missed_cell_map_polar, &total_cell_map_polar);

  total_cell_histo.Add(&cell_histo, &missed_cell_histo);
  efficiency_cell_histo.Divide(&missed_cell_histo, &total_cell_histo);
    
  TCanvas *C = new TCanvas ("C", "", 1000, 500);
  C->Divide(2,1);
  C->cd(1)->SetLogz();
  cell_map.Draw("colz");
  cell_box.Draw();
  C->cd(2)->SetLogz();
  missed_cell_map.Draw("colz"); 
  cell_box.Draw();

  (new TCanvas ("C2", "", 500, 500))->SetLogz();
  efficiency_cell_map.Draw("colz");
  cell_box.Draw();

  (new TCanvas ("C3", "", 1500, 500))->SetLogz();
  efficiency_cell_histo.SetTitle(";Cell Row;Cell Layer");
  efficiency_cell_histo.GetYaxis()->SetBinLabel(1, "FR");
  efficiency_cell_histo.GetYaxis()->SetTickSize(0.005);
  for (int fr_layer=8; fr_layer>=0; fr_layer-=2)
    efficiency_cell_histo.GetYaxis()->SetBinLabel(2+(8-fr_layer), Form("%d", fr_layer));
  for (int it_layer=0; it_layer<9; it_layer+=2)
    efficiency_cell_histo.GetYaxis()->SetBinLabel(12+it_layer, Form("%d", it_layer));
  efficiency_cell_histo.GetYaxis()->SetBinLabel(21, "IT");
  for (int row=0; row<113; row+=2)
    efficiency_cell_histo.GetXaxis()->SetBinLabel(2+row, Form("%d", row));
  efficiency_cell_histo.GetZaxis()->SetRangeUser(0, 1);
  efficiency_cell_histo.SetContour(100);
  efficiency_cell_histo.Draw("colz");

  // for (int gg_side=0; gg_side<2; gg_side++)

  // int gg_side = 0;
  for (int gg_row=0; gg_row<113; gg_row++)
    for (int gg_layer=0; gg_layer<9; gg_layer++)
      if (efficiency_cell_histo.GetBinContent(2+gg_row, 12+gg_layer) > 0.5)
	printf("%d, ", 0 + gg_row*9 + gg_layer);

  // int gg_side = 1;
  for (int gg_row=0; gg_row<113; gg_row++)
    for (int gg_layer=0; gg_layer<9; gg_layer++)
      if (efficiency_cell_histo.GetBinContent(2+gg_row, 10-gg_layer) > 0.5)
	printf("%d, ", 1017 + gg_row*9 + gg_layer);

  printf("\n");



  std::vector<float> radiuses;
  radiuses.reserve(360);

  // TGraph *efficiency_vs_radius = new TGraph;

  for (int bin=1; bin<=efficiency_cell_map_polar.GetNbinsX(); bin++) {

    TH1D *histo1d = efficiency_cell_map_polar.ProjectionY("histo1d", 1+bin, 1+bin);

    int binx=1;

    for (binx=1; binx<=histo1d->GetNbinsX(); binx++) {
      if (histo1d->GetBinContent(binx) > 0.5) break;
    }

    if (binx > histo1d->GetNbinsX())
      continue;

    // printf("binx = %d\n", binx);

    const float x1 = histo1d->GetBinCenter(binx-1);
    const float x2 = histo1d->GetBinCenter(binx);
    const float y1 = histo1d->GetBinContent(binx-1);
    const float y2 = histo1d->GetBinContent(binx);

    const float slope = (y2 - y1) / (x2 - x1);
    const float intercept = y1 - slope * x1;

    // efficiency_vs_radius->SetPoint(efficiency_vs_radius->GetN(),
    // 				   efficiency_cell_map_polar.GetXaxis()->GetBinCenter(bin),
    // 				   (0.5 - intercept)/slope);

    radiuses.push_back((0.5-intercept)/slope);
  
    delete histo1d;
  }

  // const float radius_mean = ;
  const float radius_mean = TMath::Mean(radiuses.begin(), radiuses.end());
  const float radius_rms = TMath::RMS(radiuses.begin(), radiuses.end());
  printf("%.2f %.2f\n", radius_mean, radius_rms);

  // new TCanvas;
  // efficiency_vs_radius->Draw("al");

}




void cell_inefficiency_per_run()
{
  gStyle->SetOptTitle(1);

  bool is_dead[2034];
  memset(is_dead, 0, sizeof(is_dead));
  for (const uint16_t dead_cell : dead_cells)
    is_dead[dead_cell] = true;

  TGraphErrors *radius_vs_run = new TGraphErrors;
  
  ttd_residual_data *ttd_data = new ttd_residual_data;
  ttd_residual->SetBranchAddress("ttd", &ttd_data);

  const int ttd_entries = ttd_residual->GetEntries();

  uint32_t last_run = 0;

  for (int entry=0; entry<ttd_entries; entry++) {

    ttd_residual->GetEntry(entry);

    if (last_run != ttd_data->run) {

      if (last_run > 0) { 

	total_cell_map.Add(&cell_map, &missed_cell_map);
	efficiency_cell_map.Divide(&missed_cell_map, &total_cell_map);

	total_cell_map_polar.Add(&cell_map_polar, &missed_cell_map_polar);
	efficiency_cell_map_polar.Divide(&missed_cell_map_polar, &total_cell_map_polar);

	total_cell_histo.Add(&cell_histo, &missed_cell_histo);
	efficiency_cell_histo.Divide(&missed_cell_histo, &total_cell_histo);

	std::vector<float> radiuses;
	radiuses.reserve(360);

	for (int bin=1; bin<=efficiency_cell_map_polar.GetNbinsX(); bin++) {

	  TH1D *histo1d = efficiency_cell_map_polar.ProjectionY("histo1d", 1+bin, 1+bin);

	  int binx=1;

	  for (binx=1; binx<=histo1d->GetNbinsX(); binx++) {
	    if (histo1d->GetBinContent(binx) > 0.5) break;
	  }

	  if (binx > histo1d->GetNbinsX())
	    continue;

	  const float x1 = histo1d->GetBinCenter(binx-1);
	  const float x2 = histo1d->GetBinCenter(binx);
	  const float y1 = histo1d->GetBinContent(binx-1);
	  const float y2 = histo1d->GetBinContent(binx);
	  const float slope = (y2 - y1) / (x2 - x1);
	  const float intercept = y1 - slope * x1;

	  radiuses.push_back((0.5-intercept)/slope);
  
	  delete histo1d;
	}

	const float radius_mean = TMath::Mean(radiuses.begin(), radiuses.end());
	const float radius_rms = TMath::RMS(radiuses.begin(), radiuses.end());
	printf("%d %.2f %.2f\n", ttd_data->run, radius_mean, radius_rms);

	const int point= radius_vs_run->GetN();
	radius_vs_run->SetPoint(point, ttd_data->run, radius_mean);
	radius_vs_run->SetPointError(point, 0, radius_rms);
      }

      cell_map.Reset(); missed_cell_map.Reset(); total_cell_map.Reset(); efficiency_cell_map.Reset();
      cell_map_polar.Reset(); missed_cell_map_polar.Reset(); total_cell_map_polar.Reset(); efficiency_cell_map.Reset();
      cell_histo.Reset(); missed_cell_histo.Reset(); total_cell_histo.Reset(); efficiency_cell_histo.Reset();

      last_run = ttd_data->run;
    }

    for (const track_data & track_d : ttd_data->tracks) {

      if (track_d.cells.size() < 5)
	continue;

      if (track_d.chi2ndf > 2)
	continue;

      if ((track_d.flag & 0x02) != 0)
	continue; // not delayed

      for (const cell_data & cell_d : track_d.cells) {

	if (is_dead[cell_d.gg_num])
	  continue;

	cell_map.Fill(cell_d.xfit, cell_d.yfit);

	const float cell_d_r_fit = std::sqrt(std::pow(cell_d.xfit,2)+std::pow(cell_d.yfit,2));
	const float cell_d_theta_fit = std::atan2(cell_d.yfit,cell_d.xfit) * (180.0/TMath::Pi());
	cell_map_polar.Fill(cell_d_theta_fit, cell_d_r_fit);

	const uint16_t gg_side = (cell_d.gg_num < 1017) ? 0 : 1;
	const uint16_t gg_row = (cell_d.gg_num - gg_side*1017) / 9;
	const uint16_t gg_layer = cell_d.gg_num % 9;
	if (gg_side == 0) cell_histo.Fill(0.5+gg_row, 11.5+gg_layer);
	else cell_histo.Fill(0.5+gg_row, 9.5-gg_layer);	
      }

      for (const cell_data & cell_d : track_d.missed_cells) {

	if (is_dead[cell_d.gg_num])
	  continue;

	missed_cell_map.Fill(cell_d.xfit, cell_d.yfit);

	const float cell_d_r_fit = std::sqrt(std::pow(cell_d.xfit,2)+std::pow(cell_d.yfit,2));
	const float cell_d_theta_fit = std::atan2(cell_d.yfit,cell_d.xfit) * (180.0/TMath::Pi());
	missed_cell_map_polar.Fill(cell_d_theta_fit, cell_d_r_fit);

	const uint16_t gg_side = (cell_d.gg_num < 1017) ? 0 : 1;
	const uint16_t gg_row = (cell_d.gg_num - gg_side*1017) / 9;
	const uint16_t gg_layer = cell_d.gg_num % 9;
	if (gg_side == 0) missed_cell_histo.Fill(0.5+gg_row, 11.5+gg_layer);
	else missed_cell_histo.Fill(0.5+gg_row, 9.5-gg_layer);	
      }
    } // for (track_d)

  } // for (entry)

  new TCanvas;
  radius_vs_run->GetYaxis()->SetRangeUser(0, 30);
  radius_vs_run->Draw("apl");

}

  

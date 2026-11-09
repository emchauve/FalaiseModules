TLine fit_om_lines[4];

void fit_om_vertex (TH2F *om_vtx_histo)
{
  gStyle->SetOptTitle(1);

  TH1D *om_vtx_histo_px = om_vtx_histo->ProjectionX("om_vtx_histo_px"); // Y
  om_vtx_histo_px->SetTitle(Form("%s Y", om_vtx_histo->GetTitle()));

  TH1D *om_vtx_histo_px_deriv = new TH1D(*om_vtx_histo_px);
  om_vtx_histo_px_deriv->SetTitle("Y derivative"); // Form("%s Y derivative", om_vtx_histo->GetTitle()));
  om_vtx_histo_px_deriv->Reset();

  for (int bin=1+3; bin<=om_vtx_histo_px->GetNbinsX()-3; bin++) {
    float content = 0;
    content += om_vtx_histo_px->GetBinContent(bin+3);
    content += om_vtx_histo_px->GetBinContent(bin+2);
    content += om_vtx_histo_px->GetBinContent(bin+1);
    content -= om_vtx_histo_px->GetBinContent(bin-1);
    content -= om_vtx_histo_px->GetBinContent(bin-2);
    content -= om_vtx_histo_px->GetBinContent(bin-3);
    om_vtx_histo_px_deriv->SetBinContent(bin, content);
  }

  const int px_min_bin = om_vtx_histo_px_deriv->GetMaximumBin();
  const int px_max_bin = om_vtx_histo_px_deriv->GetMinimumBin();

  //

  TH1D *om_vtx_histo_py = om_vtx_histo->ProjectionY("om_vtx_histo_py"); // Z
  om_vtx_histo_py->SetTitle(Form("%s Z", om_vtx_histo->GetTitle()));

  TH1D *om_vtx_histo_py_deriv = new TH1D(*om_vtx_histo_py);
  om_vtx_histo_py_deriv->SetTitle("Z derivative"); // Form("%s Z derivative", om_vtx_histo->GetTitle()));
  om_vtx_histo_py_deriv->Reset();

  for (int bin=1+3; bin<=om_vtx_histo_py->GetNbinsX()-3; bin++) {
    float content = 0;
    content += om_vtx_histo_py->GetBinContent(bin+3);
    content += om_vtx_histo_py->GetBinContent(bin+2);
    content += om_vtx_histo_py->GetBinContent(bin+1);
    content -= om_vtx_histo_py->GetBinContent(bin-1);
    content -= om_vtx_histo_py->GetBinContent(bin-2);
    content -= om_vtx_histo_py->GetBinContent(bin-3);
    om_vtx_histo_py_deriv->SetBinContent(bin, content);
  }

  const int py_min_bin = om_vtx_histo_py_deriv->GetMaximumBin();
  const int py_max_bin = om_vtx_histo_py_deriv->GetMinimumBin();

  // build px2 histo and its derivative

  TH1D *om_vtx_histo_px2 = om_vtx_histo->ProjectionX("om_vtx_histo_px2", py_min_bin, py_max_bin); // Y
  om_vtx_histo_px2->SetTitle(Form("%s Y", om_vtx_histo->GetTitle()));

  TH1D *om_vtx_histo_px2_deriv = new TH1D(*om_vtx_histo_px);
  om_vtx_histo_px2_deriv->SetTitle("Y derivative"); // Form("%s Y derivative", om_vtx_histo->GetTitle()));
  om_vtx_histo_px2_deriv->Reset();

  for (int bin=1+3; bin<=om_vtx_histo_px2->GetNbinsX()-3; bin++) {
    float content = 0;
    content += om_vtx_histo_px2->GetBinContent(bin+3);
    content += om_vtx_histo_px2->GetBinContent(bin+2);
    content += om_vtx_histo_px2->GetBinContent(bin+1);
    content -= om_vtx_histo_px2->GetBinContent(bin-1);
    content -= om_vtx_histo_px2->GetBinContent(bin-2);
    content -= om_vtx_histo_px2->GetBinContent(bin-3);
    om_vtx_histo_px2_deriv->SetBinContent(bin, content);
  }

  const int px2_min_bin = om_vtx_histo_px2_deriv->GetMaximumBin();
  const int px2_max_bin = om_vtx_histo_px2_deriv->GetMinimumBin();

  // build py2 histo and its derivative

  TH1D *om_vtx_histo_py2 = om_vtx_histo->ProjectionY("om_vtx_histo_py2", px_min_bin, px_max_bin); // Z
  om_vtx_histo_py2->SetTitle(Form("%s Z", om_vtx_histo->GetTitle()));

  TH1D *om_vtx_histo_py2_deriv = new TH1D(*om_vtx_histo_py);
  om_vtx_histo_py2_deriv->SetTitle("Y derivative"); // Form("%s Y derivative", om_vtx_histo->GetTitle()));
  om_vtx_histo_py2_deriv->Reset();

  for (int bin=1+3; bin<=om_vtx_histo_py2->GetNbinsX()-3; bin++) {
    float content = 0;
    content += om_vtx_histo_py2->GetBinContent(bin+3);
    content += om_vtx_histo_py2->GetBinContent(bin+2);
    content += om_vtx_histo_py2->GetBinContent(bin+1);
    content -= om_vtx_histo_py2->GetBinContent(bin-1);
    content -= om_vtx_histo_py2->GetBinContent(bin-2);
    content -= om_vtx_histo_py2->GetBinContent(bin-3);
    om_vtx_histo_py2_deriv->SetBinContent(bin, content);
  }

  const int py2_min_bin = om_vtx_histo_py2_deriv->GetMaximumBin();
  const int py2_max_bin = om_vtx_histo_py2_deriv->GetMinimumBin();

  const float px2_min = om_vtx_histo_px2_deriv->GetBinCenter(px2_min_bin);
  const float px2_max = om_vtx_histo_px2_deriv->GetBinCenter(px2_max_bin);

  const float py2_min = om_vtx_histo_py2_deriv->GetBinCenter(py2_min_bin);
  const float py2_max = om_vtx_histo_py2_deriv->GetBinCenter(py2_max_bin);

  // printf("Y width = %4.1f cm   Z-width = %4.1f\n", 100*(px2_max-px2_min), 100*(py2_max-py2_min));

  //

  TCanvas *c1d = (TCanvas *)(gROOT->FindObject("c1d"));

  if (c1d == nullptr) {
    c1d = new TCanvas ("c1d", "", 700+350, 500+250);
    c1d->Divide(2, 2);
  }

  c1d->cd(1);
  om_vtx_histo_px->SetLineColor(kBlack);
  om_vtx_histo_px->SetLineWidth(2);
  om_vtx_histo_px->Draw();
  om_vtx_histo_px2->SetLineColor(kGreen+2);
  om_vtx_histo_px2->SetLineWidth(2);
  om_vtx_histo_px2->Draw("same");

  c1d->cd(2);
  om_vtx_histo_py->SetLineColor(kBlack);
  om_vtx_histo_py->SetLineWidth(2);
  om_vtx_histo_py->Draw();
  om_vtx_histo_py2->SetLineColor(kGreen+2);
  om_vtx_histo_py2->SetLineWidth(2);
  om_vtx_histo_py2->Draw("same");

  c1d->cd(3);
  om_vtx_histo_px_deriv->SetLineColor(kBlack);
  om_vtx_histo_px_deriv->SetLineWidth(2);
  om_vtx_histo_px_deriv->Draw();
  om_vtx_histo_px2_deriv->SetLineColor(kRed+1);
  om_vtx_histo_px2_deriv->SetLineWidth(2);
  om_vtx_histo_px2_deriv->Draw("same");

  c1d->cd(4);
  om_vtx_histo_py_deriv->SetLineColor(kBlack);
  om_vtx_histo_py_deriv->SetLineWidth(2);
  om_vtx_histo_py_deriv->Draw();
  om_vtx_histo_py2_deriv->SetLineColor(kRed+1);
  om_vtx_histo_py2_deriv->SetLineWidth(2);
  om_vtx_histo_py2_deriv->Draw("same");


  c1d->cd(3);

  TF1 *px2_fit = new TF1("px2_fit", "gaus", -3, 3);
  px2_fit->SetNpx(600);

  px2_fit->SetParameter(0, om_vtx_histo_py2_deriv->GetMaximum());
  px2_fit->SetParameter(1, px2_min);
  px2_fit->SetParameter(2, 0.025);
  px2_fit->SetRange(px2_min-0.025, px2_min+0.025);
  om_vtx_histo_px2_deriv->Fit("px2_fit", "0QR");
  float px3_min = px2_fit->GetParameter(1);
  px2_fit->SetRange(px3_min-0.025, px3_min+0.025);
  om_vtx_histo_px2_deriv->Fit("px2_fit", "0QR");
  px3_min = px2_fit->GetParameter(1);

  px2_fit->SetParameter(0, om_vtx_histo_py2_deriv->GetMinimum());
  px2_fit->SetParameter(1, px2_max);
  px2_fit->SetParameter(2, 0.025);
  px2_fit->SetRange(px2_max-0.025, px2_max+0.025);
  om_vtx_histo_px2_deriv->Fit("px2_fit", "0QR");
  float px3_max = px2_fit->GetParameter(1);
  px2_fit->SetRange(px3_max-0.025, px3_max+0.025);
  om_vtx_histo_px2_deriv->Fit("px2_fit", "0QR");
  px3_max = px2_fit->GetParameter(1);

  c1d->cd(4);

  TF1 *py2_fit = new TF1("py2_fit", "gaus", -3, 3);
  py2_fit->SetNpx(600);

  py2_fit->SetParameter(0, om_vtx_histo_py2_deriv->GetMaximum());
  py2_fit->SetParameter(1, py2_min);
  py2_fit->SetParameter(2, 0.025);
  py2_fit->SetRange(py2_min-0.025, py2_min+0.025);
  om_vtx_histo_py2_deriv->Fit("py2_fit", "0QR");
  float py3_min = py2_fit->GetParameter(1);
  py2_fit->SetRange(py3_min-0.025, py3_min+0.025);
  om_vtx_histo_py2_deriv->Fit("py2_fit", "0QR");
  py3_min = py2_fit->GetParameter(1);

  py2_fit->SetParameter(0, om_vtx_histo_py2_deriv->GetMinimum());
  py2_fit->SetParameter(1, py2_max);
  py2_fit->SetParameter(2, 0.025);
  py2_fit->SetRange(py2_max-0.025, py2_max+0.025);
  om_vtx_histo_py2_deriv->Fit("py2_fit", "0QR");
  float py3_max = py2_fit->GetParameter(1);
  py2_fit->SetRange(py3_max-0.025, py3_max+0.025);
  om_vtx_histo_py2_deriv->Fit("py2_fit", "0QR");
  py3_max = py2_fit->GetParameter(1);

  printf("Y width = %4.1f cm   Z-width = %4.1f\n", 100*(px3_max-px3_min), 100*(py3_max-py3_min));

  //

  TCanvas *c2d = (TCanvas *)(gROOT->FindObject("c2d"));

  if (c2d == nullptr)
    c2d = new TCanvas ("c2d", "", 700, 500);

  c2d->cd();
  om_vtx_histo->Draw("colz");

  fit_om_lines[0].SetLineWidth(2);
  fit_om_lines[0].SetX1(px3_min);
  fit_om_lines[0].SetX2(px3_min);
  fit_om_lines[0].SetY1(py3_min - (py3_max-py3_min));
  fit_om_lines[0].SetY2(py3_max + (py3_max-py3_min));
  fit_om_lines[0].Draw();

  fit_om_lines[1].SetLineWidth(2);
  fit_om_lines[1].SetX1(px3_max);
  fit_om_lines[1].SetX2(px3_max);
  fit_om_lines[1].SetY1(py3_min - (py3_max-py3_min));
  fit_om_lines[1].SetY2(py3_max + (py3_max-py3_min));
  fit_om_lines[1].Draw();

  fit_om_lines[2].SetLineWidth(2);
  fit_om_lines[2].SetX1(px3_min - (px3_max-px3_min));
  fit_om_lines[2].SetX2(px3_max + (px3_max-px3_min));
  fit_om_lines[2].SetY1(py3_min);
  fit_om_lines[2].SetY2(py3_min);
  fit_om_lines[2].Draw();

  fit_om_lines[3].SetLineWidth(2);
  fit_om_lines[3].SetX1(px3_min - (px3_max-px3_min));
  fit_om_lines[3].SetX2(px3_max + (px3_max-px3_min));
  fit_om_lines[3].SetY1(py3_max);
  fit_om_lines[3].SetY2(py3_max);
  fit_om_lines[3].Draw();

  om_vtx_histo->GetName();
  const char *om_char = &(om_vtx_histo->GetName()[13]);
  const int om_num = std::atoi(om_char);

  // std::ofstream fit_output ("fit-vertex.txt", std::ios::app);
  // fit_output << om_num << "\t"
  // 	     << 0.5*(px3_min+px3_max) << "\t"
  // 	     << 0.5*(py3_min+py3_max) << "\t"
  // 	     << px3_max-px3_min << "\t"
  // 	     << py3_max-py3_min << "\t"
  // 	     << std::endl;
  // fit_output.close();

}


void fit_all_om_vertex()
{
  for (int om=0; om<520; om++) {

    if ((om % 13) == 0) continue;
    if ((om % 13) == 12) continue;

    TH2F *om_vtx_histo = (TH2F*)(gROOT->FindObject(Form("om_vtx_histo_%03d",om)));
    if (om_vtx_histo->GetEntries() < 10000) continue;
    fit_om_vertex(om_vtx_histo);
  }

}

TH2F mw0_vtx_histo ("mw0_vtx_histo", "MW0 VERTEX;Y (m); Z (m)", 600, -3, 3, 400, -2, 2);
TH2F mw1_vtx_histo ("mw1_vtx_histo", "MW1 VERTEX;Y (m); Z (m)", 600, -3, 3, 400, -2, 2);

TH2F mw0_widthx_histo ("mw0_widthx_histo", "", 20, 0, 20, 13, 0, 13);
TH2F mw0_widthy_histo ("mw0_widthy_histo", "", 20, 0, 20, 13, 0, 13);
TH2F mw1_widthx_histo ("mw1_widthx_histo", "", 20, 0, 20, 13, 0, 13);
TH2F mw1_widthy_histo ("mw1_widthy_histo", "", 20, 0, 20, 13, 0, 13);

const float mc_mw_posx[20] = {-2460.5, -2201.5, -1942.5, -1683.5, -1424.5, -1165.5, -906.5, -647.5, -388.5, -129.5,
			   +129.5, +388.5, +647.5, +906.5, +1165.5, +1424.5, +1683.5, +1942.5, +2201.5, +2460.5};

const float mc_mw_posy[13] = {-1554, -1295, -1036, -777, -518, -259, 0, 259, 518, 777, 1036, 1295, 1554};

void draw_om_vertex()
{
  TCanvas *mw0_vtx_canvas = new TCanvas("mw0_vtx_canvas", "", 700, 500);
  mw0_vtx_histo.Draw("colz");

  TCanvas *mw1_vtx_canvas = new TCanvas("mw1_vtx_canvas", "", 700, 500);
  mw1_vtx_histo.Draw("colz");

  std::ifstream fit_output ("fit-vertex.txt");

  int om_num;
  float data_om_posx, data_om_posy;
  float data_om_widthx, data_om_widthy;

  float data_it_mw_posx[20];
  float data_it_mw_posx_stat[20];
  memset(data_it_mw_posx, 0, sizeof(data_it_mw_posx));
  memset(data_it_mw_posx_stat, 0, sizeof(data_it_mw_posx_stat));

  float data_it_mw_posy[13];
  float data_it_mw_posy_stat[13];
  memset(data_it_mw_posy, 0, sizeof(data_it_mw_posy));
  memset(data_it_mw_posy_stat, 0, sizeof(data_it_mw_posy_stat));

  float data_fr_mw_posx[20];
  float data_fr_mw_posx_stat[20];
  memset(data_fr_mw_posx, 0, sizeof(data_fr_mw_posx));
  memset(data_fr_mw_posx_stat, 0, sizeof(data_fr_mw_posx_stat));

  float data_fr_mw_posy[13];
  float data_fr_mw_posy_stat[13];
  memset(data_fr_mw_posy, 0, sizeof(data_fr_mw_posy));
  memset(data_fr_mw_posy_stat, 0, sizeof(data_fr_mw_posy_stat));

  while (fit_output >> om_num) {

    fit_output >> data_om_posx >> data_om_posy >> data_om_widthx >> data_om_widthy;

    if (om_num < 260) {

      mw0_vtx_canvas->cd();

      const int om_column = om_num/13;
      const int om_row = om_num%13;
      mw0_widthx_histo.SetBinContent(1+om_column, 1+om_row, data_om_widthx);
      mw0_widthy_histo.SetBinContent(1+om_column, 1+om_row, data_om_widthy);

      TMarker *mc_marker = new TMarker(mc_mw_posx[om_column]/1E3, mc_mw_posy[om_row]/1E3, 20);
      mc_marker->SetMarkerColor(kRed+1);
      mc_marker->Draw();

      TMarker *data_marker = new TMarker(data_om_posx, data_om_posy, 2);
      data_marker->SetMarkerSize(2);
      data_marker->Draw();

      if ((om_column > 0) && (om_column < 19)) {
	data_it_mw_posx[om_column] += data_om_posx;
	data_it_mw_posx_stat[om_column]++;
      }

      if ((om_row > 1) && (om_row < 11)) {
	data_it_mw_posy[om_row] += data_om_posy;
	data_it_mw_posy_stat[om_row]++;
      }

    }

    else {
      mw1_vtx_canvas->cd();

      const int om_column = (om_num-260)/13;
      const int om_row = om_num%13;
      mw1_widthx_histo.SetBinContent(1+om_column, 1+om_row, data_om_widthx);
      mw1_widthy_histo.SetBinContent(1+om_column, 1+om_row, data_om_widthy);

      TMarker *mc_marker = new TMarker(mc_mw_posx[om_column]/1E3, mc_mw_posy[om_row]/1E3, 20);
      mc_marker->SetMarkerColor(kRed+1);
      mc_marker->Draw();

      TMarker *data_marker = new TMarker(data_om_posx, data_om_posy, 2);
      data_marker->SetMarkerSize(2);
      data_marker->Draw();

      if ((om_column > 0) && (om_column < 19)) {
	data_fr_mw_posx[om_column] += data_om_posx;
	data_fr_mw_posx_stat[om_column]++;
      }

      if ((om_row > 1) && (om_row < 11)) {
	data_fr_mw_posy[om_row] += data_om_posy;
	data_fr_mw_posy_stat[om_row]++;
      }

    }


    TBox *box = new TBox(data_om_posx-data_om_widthx/2,
			 data_om_posy-data_om_widthy/2,
			 data_om_posx+data_om_widthx/2,
			 data_om_posy+data_om_widthy/2);
    box->SetFillStyle(0);
    // box->SetLineColor(kBlack);
    // box->SetLineWidth(2);
    box->Draw();

  }

  fit_output.close();


  for (int om_column=0; om_column<19; om_column++) {

    mw0_vtx_canvas->cd();

    if (data_it_mw_posx_stat[om_column] != 0) {
      data_it_mw_posx[om_column] /= data_it_mw_posx_stat[om_column];
      TLine *line = new TLine (data_it_mw_posx[om_column], -1.5, data_it_mw_posx[om_column], 1.5);
      line->SetLineStyle(kDotted); line->Draw();
    }

    mw1_vtx_canvas->cd();

    if (data_fr_mw_posx_stat[om_column] != 0) {
      data_fr_mw_posx[om_column] /= data_fr_mw_posx_stat[om_column];
      TLine *line = new TLine (data_fr_mw_posx[om_column], -1.5, data_fr_mw_posx[om_column], 1.5);
      line->SetLineStyle(kDotted); line->Draw();
    }

  }

  for (int om_row=0; om_row<13; om_row++) {

    mw0_vtx_canvas->cd();

    if (data_it_mw_posy_stat[om_row] != 0) {
      data_it_mw_posy[om_row] /= data_it_mw_posy_stat[om_row];
      TLine *line = new TLine (-2.75, data_it_mw_posy[om_row], 2.75, data_it_mw_posy[om_row]);
      line->SetLineStyle(kDotted); line->Draw();
    }

    mw1_vtx_canvas->cd();

    if (data_fr_mw_posy_stat[om_row] != 0) {
      data_fr_mw_posy[om_row] /= data_fr_mw_posy_stat[om_row];
      TLine *line = new TLine (-2.75, data_fr_mw_posy[om_row], 2.75, data_fr_mw_posy[om_row]);
      line->SetLineStyle(kDotted); line->Draw();
    }

  }

  printf("COLUMN :     MC    DATA-IT  DATA-FR   <DATA>\n");
  // printf("  [XX] :  XXXX.X  XXXX.X  XXXX.X");
  for (int om_column=0; om_column<20; om_column++) {
    printf("  [%02d] :  %7.1f  %7.1f  %7.1f  %7.1f\n",
	   om_column, mc_mw_posx[om_column],
	   1000*data_it_mw_posx[om_column], 1000*data_fr_mw_posx[om_column],
	   1000*(data_it_mw_posx[om_column]+data_fr_mw_posx[om_column])/2);
  }

  printf("\n");

  printf(" ROW   :     MC    DATA-IT  DATA-FR   <DATA>\n");
  // printf("  [XX] :  XXXX.X  XXXX.X  XXXX.X");
  for (int om_row=0; om_row<13; om_row++) {
    printf("  [%02d] :  %7.1f  %7.1f  %7.1f  %7.1f\n",
	   om_row, mc_mw_posy[om_row],
	   1000*data_it_mw_posy[om_row], 1000*data_fr_mw_posy[om_row],
	   1000*(data_it_mw_posy[om_row]+data_fr_mw_posy[om_row])/2);
    if ((om_row > 1) && (om_row < 10))
      printf("                                            %+7.1f\n",
	     1000*(data_it_mw_posy[om_row+1]+data_fr_mw_posy[om_row+1])/2-1000*(data_it_mw_posy[om_row]+data_fr_mw_posy[om_row])/2);
  }


}

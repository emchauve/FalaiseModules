TLine fit_om_lines[4];

void fit_vertex (TH2F *om_vtx_histo)
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

  //

  TCanvas *c2d = (TCanvas *)(gROOT->FindObject("c2d"));

  if (c2d == nullptr)
    c2d = new TCanvas ("c2d", "", 700, 500);

  c2d->cd();
  om_vtx_histo->Draw("colz");

  const float px2_min = om_vtx_histo_px2_deriv->GetBinCenter(px2_min_bin);
  const float px2_max = om_vtx_histo_px2_deriv->GetBinCenter(px2_max_bin);

  const float py2_min = om_vtx_histo_py2_deriv->GetBinCenter(py2_min_bin);
  const float py2_max = om_vtx_histo_py2_deriv->GetBinCenter(py2_max_bin);

  printf("Y width = %4.1f cm   Z-width = %4.1f\n", 100*(px2_max-px2_min), 100*(py2_max-py2_min));

  fit_om_lines[0].SetLineWidth(2);
  fit_om_lines[0].SetX1(px2_min);
  fit_om_lines[0].SetX2(px2_min);
  fit_om_lines[0].SetY1(py2_min - (py2_max-py2_min));
  fit_om_lines[0].SetY2(py2_max + (py2_max-py2_min));
  fit_om_lines[0].Draw();

  fit_om_lines[1].SetLineWidth(2);
  fit_om_lines[1].SetX1(px2_max);
  fit_om_lines[1].SetX2(px2_max);
  fit_om_lines[1].SetY1(py2_min - (py2_max-py2_min));
  fit_om_lines[1].SetY2(py2_max + (py2_max-py2_min));
  fit_om_lines[1].Draw();

  fit_om_lines[2].SetLineWidth(2);
  fit_om_lines[2].SetX1(px2_min - (px2_max-px2_min));
  fit_om_lines[2].SetX2(px2_max + (px2_max-px2_min));
  fit_om_lines[2].SetY1(py2_min);
  fit_om_lines[2].SetY2(py2_min);
  fit_om_lines[2].Draw();

  fit_om_lines[3].SetLineWidth(2);
  fit_om_lines[3].SetX1(px2_min - (px2_max-px2_min));
  fit_om_lines[3].SetX2(px2_max + (px2_max-px2_min));
  fit_om_lines[3].SetY1(py2_max);
  fit_om_lines[3].SetY2(py2_max);
  fit_om_lines[3].Draw();

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
}

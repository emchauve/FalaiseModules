double _tracker_drift_model_manu_params_[5] = {0.263223, -0.030965, -0.571594, 6.01392e-02, 1.13142e+03};

double _tracker_drift_model_manu_ (double time_usec)
{
  if (time_usec < -1)
    return 0;

  else if (time_usec < 0)
    return 1E-9;

  const double f2 = _tracker_drift_model_manu_params_[3] * TMath::Log(1 + time_usec * _tracker_drift_model_manu_params_[4]);

  if (time_usec > 10)
    return f2;

  const double f1 = _tracker_drift_model_manu_params_[0] * TMath::Exp(_tracker_drift_model_manu_params_[1]*time_usec)/TMath::Power(time_usec,_tracker_drift_model_manu_params_[2]);

  return TMath::Min(f1, f2);
}

double _tracker_drift_model_manu_func_ (double *x, double *)
{
  return _tracker_drift_model_manu_ (x[0]) * 44;
}

double _tracker_drift_model_error_manu_params_[4] = {3.41215e-02, 2.87786e-04, -1.75261e+01, 6.33008e-01};

double _tracker_drift_model_error_manu_ (double time_usec)
{
  if (time_usec < 0)
    time_usec = 0;

  double value = _tracker_drift_model_error_manu_params_[0];
  value += _tracker_drift_model_error_manu_params_[1] * time_usec;
  value += TMath::Exp(_tracker_drift_model_error_manu_params_[2] + _tracker_drift_model_error_manu_params_[3]*time_usec);

  return value;
}

double _tracker_drift_model_error_manu_func_ (double *x, double *)
{
  return _tracker_drift_model_error_manu_ (x[0]);
}

TH2F r_drift_histo ("r_drift_histo", ";Drift time (us);Rfit (mm)", 3000, 0, 15, 1500, 0, 30);
TH2F r_residual_histo ("r_residual_histo", ";Rfit (mm);Rcd - Rfit (mm)", 300, 0, 30, 200, -10, 10);
TH2F z_residual_histo ("z_residual_histo", ";Zfit (m);Zcd - Zfit (cm)", 300, -1.5, 1.5, 250, -5, 5);

TF1 drift_model_manu ("drift_model_manu", _tracker_drift_model_manu_func_, 0, 100, 0);
drift_model_manu.SetLineColor(kBlack);
drift_model_manu.SetLineWidth(2);
drift_model_manu.SetNpx(1000);

TF1 drift_model_error_manu ("drift_model_error_manu", _tracker_drift_model_error_manu_func_, 0, 100, 0);
drift_model_error_manu.SetLineColor(kBlack);
drift_model_error_manu.SetLineWidth(2);
drift_model_error_manu.SetNpx(1000);

void draw_residual()
{
  gStyle->SetOptTitle(1);

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

      // if (std::fabs(track_d.theta) > 1)
      // 	continue; // horizontal tracks only

      // if (std::fabs(track_d.phi) > 1)
      // 	continue; // horizontal tracks only

      if (track_d.cells.size() < 8)
	continue;

      if (track_d.chi2ndf > 2)
	continue;

      if ((track_d.flag & 0x02) != 0)
	continue; // skip delayed

      for (const cell_data & cell_d : track_d.cells) {

	if ((cell_d.flag & 0x3) != 0x3)
	  continue; // missing cathode(s)

	// printf("R = %.1f %.1f\n", cell_d.rfit, cell_d.r - cell_d.rfit);
	// printf("Z = %.1f %.1f\n", cell_d.zfit, cell_d.z - cell_d.zfit);
	r_drift_histo.Fill(cell_d.time_anode/1E3, cell_d.rfit);
	r_residual_histo.Fill(cell_d.rfit, (cell_d.r-cell_d.rfit));
	z_residual_histo.Fill(cell_d.zfit/1E3, (cell_d.z-cell_d.zfit)/1E1);

      } // for (cell_d)

    } // for (track_d)

  } // for (entry)

  printf("\n");

  (new TCanvas)->SetLogz();
  r_residual_histo.Draw("colz");

  (new TCanvas)->SetLogz();
  z_residual_histo.Draw("colz");

  (new TCanvas)->SetLogz();
  r_drift_histo.Draw("colz");

}

  

double GaussGauss(double *x, double *par)
{
  // par[0] = mean of 1st gaussian
  // par[1] = stddev of 1st gaussian
  // par[2] = mean of 2nd gaussian
  // par[3] = stddev of 2nd gaussian
  // par[4] = weight of 1st gaussian
  // par[5] = weight of 2nd gaussian
  double gauss1=exp(-0.5*pow((x[0]-par[0])/par[1], 2));
  double gauss2=exp(-0.5*pow((x[0]-par[2])/par[3], 2));
  double result=par[4]*gauss1+par[5]*gauss2;
  return result;
}

double PedGauss(double *x, double *par)
{
  // par[0] = mean of gaussian
  // par[1] = stddev of gaussian
  // par[2] = weight of gaussian
  double result=par[2]*exp(-0.5*pow((x[0]-par[0])/par[1], 2));
  return result;
}

struct FitParameters
{
  string pT_lo;
  string pT_hi;
  string eta_lo;
  string eta_hi;
  int pTbin_lo;
  int pTbin_hi;
  int etabin_lo;
  int etabin_hi;
  
  float fitRangeUS_lo, fitRangeUS_hi;
  float fitRangeSS_lo, fitRangeSS_hi;
  float parLimitUS_lo[6], parLimitUS_hi[6];
  float parLimitSS_lo[6], parLimitSS_hi[6];
};

void fit(TH3F *hOS, TH3F *hSS, FitParameters *fitParameters, std::ofstream &outfile)
{
  hOS->GetXaxis()->SetRange(fitParameters->pTbin_lo, fitParameters->pTbin_hi);
  hOS->GetYaxis()->SetRange(fitParameters->etabin_lo, fitParameters->etabin_hi);
  TH1F *h_mll_OS=(TH1F*)hOS->Project3D("z")->Clone("h_mll_OS");
  
  hSS->GetXaxis()->SetRange(fitParameters->pTbin_lo, fitParameters->pTbin_hi);
  hSS->GetYaxis()->SetRange(fitParameters->etabin_lo, fitParameters->etabin_hi);
  TH1F *h_mll_SS=(TH1F*)hSS->Project3D("z")->Clone("h_mll_SS");
  
  TH1F *h_mll_US=(TH1F*)h_mll_OS->Clone("h_mll_US");
  h_mll_US->Add(h_mll_SS);
  
  h_mll_US->Sumw2();
  h_mll_SS->Sumw2();
  
  h_mll_US->Rebin(2);
  h_mll_SS->Rebin(2);
  
  h_mll_US->GetXaxis()->SetRangeUser(fitParameters->fitRangeUS_lo, fitParameters->fitRangeUS_hi);
  h_mll_SS->GetXaxis()->SetRangeUser(fitParameters->fitRangeSS_lo, fitParameters->fitRangeSS_hi);
  
  TF1 *f_gaussgauss_US=new TF1("f_gaussgauss_US", GaussGauss, fitParameters->fitRangeUS_lo, fitParameters->fitRangeUS_hi, 6);
  // h_mll_US->Fit("gaus");
  f_gaussgauss_US->SetParLimits(0, fitParameters->parLimitUS_lo[0], fitParameters->parLimitUS_hi[0]);
  f_gaussgauss_US->SetParLimits(1, fitParameters->parLimitUS_lo[1], fitParameters->parLimitUS_hi[1]);
  f_gaussgauss_US->SetParLimits(2, fitParameters->parLimitUS_lo[2], fitParameters->parLimitUS_hi[2]);
  f_gaussgauss_US->SetParLimits(3, fitParameters->parLimitUS_lo[3], fitParameters->parLimitUS_hi[3]);
  f_gaussgauss_US->SetParLimits(4, fitParameters->parLimitUS_lo[4], fitParameters->parLimitUS_hi[4]);
  f_gaussgauss_US->SetParLimits(5, fitParameters->parLimitUS_lo[5], fitParameters->parLimitUS_hi[5]);
  h_mll_US->Fit(f_gaussgauss_US, "RL");
  
  TF1 *f_gaussgauss_US_ped=new TF1("f_gaussgauss_US_ped", PedGauss, fitParameters->fitRangeUS_lo, fitParameters->fitRangeUS_hi, 6);
  f_gaussgauss_US_ped->SetParameter(0, f_gaussgauss_US->GetParameter(2));
  f_gaussgauss_US_ped->SetParameter(1, f_gaussgauss_US->GetParameter(3));
  f_gaussgauss_US_ped->SetParameter(2, f_gaussgauss_US->GetParameter(5));
  f_gaussgauss_US_ped->SetLineColor(kGreen);
  
  TF1 *f_gaussgauss_SS=new TF1("f_gaussgauss_SS", GaussGauss, fitParameters->fitRangeSS_lo, fitParameters->fitRangeSS_hi, 6);
  // h_mll_SS->Fit("gaus");
  f_gaussgauss_SS->SetParLimits(0, fitParameters->parLimitSS_lo[0], fitParameters->parLimitSS_hi[0]);
  f_gaussgauss_SS->SetParLimits(1, fitParameters->parLimitSS_lo[1], fitParameters->parLimitSS_hi[1]);
  f_gaussgauss_SS->SetParLimits(2, fitParameters->parLimitSS_lo[2], fitParameters->parLimitSS_hi[2]);
  f_gaussgauss_SS->SetParLimits(3, fitParameters->parLimitSS_lo[3], fitParameters->parLimitSS_hi[3]);
  f_gaussgauss_SS->SetParLimits(4, fitParameters->parLimitSS_lo[4], fitParameters->parLimitSS_hi[4]);
  f_gaussgauss_SS->SetParLimits(5, fitParameters->parLimitSS_lo[5], fitParameters->parLimitSS_hi[5]);
  h_mll_SS->Fit(f_gaussgauss_SS, "RL");
  
  TF1 *f_gaussgauss_SS_ped=new TF1("f_gaussgauss_SS_ped", PedGauss, fitParameters->fitRangeSS_lo, fitParameters->fitRangeSS_hi, 6);
  f_gaussgauss_SS_ped->SetParameter(0, f_gaussgauss_SS->GetParameter(2));
  f_gaussgauss_SS_ped->SetParameter(1, f_gaussgauss_SS->GetParameter(3));
  f_gaussgauss_SS_ped->SetParameter(2, f_gaussgauss_SS->GetParameter(5));
  f_gaussgauss_SS_ped->SetLineColor(kGreen);
  
  double integral_US=f_gaussgauss_US->Integral(fitParameters->fitRangeUS_lo, fitParameters->fitRangeUS_hi);
  double integral_US_ped=f_gaussgauss_US_ped->Integral(fitParameters->fitRangeUS_lo, fitParameters->fitRangeUS_hi);
  double integral_SS=f_gaussgauss_SS->Integral(fitParameters->fitRangeSS_lo, fitParameters->fitRangeSS_hi);
  double integral_SS_ped=f_gaussgauss_SS_ped->Integral(fitParameters->fitRangeSS_lo, fitParameters->fitRangeSS_hi);
  double qFR=0.5*(integral_SS-integral_SS_ped)/(integral_US-integral_US_ped);
  double qFR_err=qFR*f_gaussgauss_SS->GetParError(4)/f_gaussgauss_SS->GetParameter(4);
  
  std::cout<<"integral_US = "<<integral_US<<std::endl;
  std::cout<<"integral_US_ped = "<<integral_US_ped<<std::endl;
  std::cout<<"integral_SS = "<<integral_SS<<std::endl;
  std::cout<<"integral_SS_ped = "<<integral_SS_ped<<std::endl;
  std::cout<<"charge flip rate = "<<qFR<<" +- "<<qFR_err<<std::endl;
  
  string pT_lo=fitParameters->pT_lo;
  string pT_hi=fitParameters->pT_hi;
  string eta_lo=fitParameters->eta_lo;
  string eta_hi=fitParameters->eta_hi;
  
  h_mll_US->SetTitle(("OS+SS ee, "+pT_lo+"<p_{T}<"+pT_hi+" GeV, "+eta_lo+"<|#eta|<"+eta_hi+"; m_{ll} (GeV)").c_str());
  h_mll_SS->SetTitle(("SS ee, "+pT_lo+"<p_{T}<"+pT_hi+" GeV, "+eta_lo+"<|#eta|<"+eta_hi+"; m_{ll} (GeV)").c_str());
  
  TCanvas *c_mll_US=new TCanvas("c_mll_US", "c_mll_US", 700, 700);
  h_mll_US->Draw();
  f_gaussgauss_US_ped->Draw("same");
  c_mll_US->SaveAs(("qFRDashboard/c_mll_US_pT_"+pT_lo+"_"+pT_hi+"_eta_"+eta_lo+"_"+eta_hi+".png").c_str());
  
  TCanvas *c_mll_SS=new TCanvas("c_mll_SS", "c_mll_SS", 700, 700);
  h_mll_SS->Draw();
  f_gaussgauss_SS_ped->Draw("same");
  c_mll_SS->SaveAs(("qFRDashboard/c_mll_SS_pT_"+pT_lo+"_"+pT_hi+"_eta_"+eta_lo+"_"+eta_hi+".png").c_str());
  
  // Write to HTML file
  outfile<<"<tr>"<<std::endl;
  outfile<<" <td>"<<std::endl;
  outfile<<"  <img src=\"c_mll_US_pT_"<<pT_lo<<"_"<<pT_hi<<"_eta_"<<eta_lo<<"_"<<eta_hi<<".png\"/> <br/>"<<std::endl;
  outfile<<"  <b> Peak Gaussian fit parameters: </b> <br/>"<<std::endl;
  outfile<<"   <blockquote> "<<std::endl;
  outfile<<"    Mean = "<<f_gaussgauss_US->GetParameter(0)<<" &plusmn "<<f_gaussgauss_US->GetParError(0)<<"</br>"<<std::endl;
  outfile<<"    Width = "<<f_gaussgauss_US->GetParameter(1)<<" &plusmn "<<f_gaussgauss_US->GetParError(1)<<"</br>"<<std::endl;
  outfile<<"    Normalization = "<<f_gaussgauss_US->GetParameter(4)<<" &plusmn "<<f_gaussgauss_US->GetParError(4)<<"<br/>"<<std::endl;
  outfile<<"   </blockquote>"<<std::endl;
  outfile<<"  <b> Pedestal Gaussian fit parameters: </b> <br/>"<<std::endl;
  outfile<<"   <blockquote> "<<std::endl;
  outfile<<"    Mean = "<<f_gaussgauss_US->GetParameter(2)<<" &plusmn "<<f_gaussgauss_US->GetParError(2)<<"</br>"<<std::endl;
  outfile<<"    Width = "<<f_gaussgauss_US->GetParameter(3)<<" &plusmn "<<f_gaussgauss_US->GetParError(3)<<"</br>"<<std::endl;
  outfile<<"    Normalization = "<<f_gaussgauss_US->GetParameter(5)<<" &plusmn "<<f_gaussgauss_US->GetParError(5)<<"<br/>"<<std::endl;
  outfile<<"   </blockquote>"<<std::endl;
  outfile<<" </td>"<<std::endl;
  outfile<<" <td>"<<std::endl;
  outfile<<"  <img src=\"c_mll_SS_pT_"<<pT_lo<<"_"<<pT_hi<<"_eta_"<<eta_lo<<"_"<<eta_hi<<".png\"/> <br/>"<<std::endl;
  outfile<<"  <b> Peak Gaussian fit parameters: </b> <br/>"<<std::endl;
  outfile<<"   <blockquote> "<<std::endl;
  outfile<<"    Mean = "<<f_gaussgauss_SS->GetParameter(0)<<" &plusmn "<<f_gaussgauss_SS->GetParError(0)<<"</br>"<<std::endl;
  outfile<<"    Width = "<<f_gaussgauss_SS->GetParameter(1)<<" &plusmn "<<f_gaussgauss_SS->GetParError(1)<<"</br>"<<std::endl;
  outfile<<"    Normalization = "<<f_gaussgauss_SS->GetParameter(4)<<" &plusmn "<<f_gaussgauss_SS->GetParError(4)<<"<br/>"<<std::endl;
  outfile<<"   </blockquote>"<<std::endl;
  outfile<<"  <b> Pedestal Gaussian fit parameters: </b> <br/>"<<std::endl;
  outfile<<"   <blockquote> "<<std::endl;
  outfile<<"    Mean = "<<f_gaussgauss_SS->GetParameter(2)<<" &plusmn "<<f_gaussgauss_SS->GetParError(2)<<"</br>"<<std::endl;
  outfile<<"    Width = "<<f_gaussgauss_SS->GetParameter(3)<<" &plusmn "<<f_gaussgauss_SS->GetParError(3)<<"</br>"<<std::endl;
  outfile<<"    Normalization = "<<f_gaussgauss_SS->GetParameter(5)<<" &plusmn "<<f_gaussgauss_SS->GetParError(5)<<"<br/>"<<std::endl;
  outfile<<"   </blockquote>"<<std::endl;
  outfile<<" </td>"<<std::endl;
  outfile<<" <td>"<<std::endl;
  outfile<<"  <b> qFR calculation </b> <br/> <br/>"<<std::endl;
  outfile<<"  <blockquote>"<<std::endl;
  outfile<<"  Unknown sign, integral under peak = "<<integral_US<<"<br/>"<<std::endl;
  outfile<<"  Unknown sign, integral under pedestal = "<<integral_US_ped<<"<br/> <br/>"<<std::endl;
  outfile<<"  Same sign, integral under peak = "<<integral_SS<<"<br/>"<<std::endl;
  outfile<<"  Same sign, integral under pedestal = "<<integral_SS_ped<<"<br/> <br/>"<<std::endl;
  outfile<<"  Relative fit uncertainty from SS peak normalization = "<<f_gaussgauss_SS->GetParError(4)/f_gaussgauss_SS->GetParameter(4)<<"</br><br/>"<<std::endl;
  outfile<<"  </blockquote>"<<std::endl;
  outfile<<"  <b> qFR = "<<qFR<<" &plusmn "<<qFR_err<<"</b><br/>"<<std::endl;
  outfile<<" </td>"<<std::endl;
  outfile<<"</tr>"<<std::endl;
  
}

void VLQ_ChargeFlipRate_Fit_ee()
{
  TFile *f_in=TFile::Open("Data_templates.root", "read");
  f_in->cd();
  
  TH3F *h_PromptRate_ee_TT_OS=(TH3F*)f_in->Get("H_PromptRate_ee_TT_OS");
  TH3F *h_PromptRate_ee_TT_SS=(TH3F*)f_in->Get("H_PromptRate_ee_TT_SS");
  
  // x is pT, 12 bins
  // y is eta, 5 bins
  // Draw the eta and pT bins
  // TH2F *h_pT_eta=(TH2F*)h_PromptRate_ee_TT_OS->Project3D("xy")->Clone("h_pT_eta");
  // TCanvas *c=new TCanvas("c", "c", 700, 700);
  // h_pT_eta->Draw("colz");
  // c->SaveAs("c.png");
  
  std::vector<FitParameters> v_fitParameters(12);
  
  v_fitParameters.at(0).pT_lo="30";
  v_fitParameters.at(0).pT_hi="40";
  v_fitParameters.at(0).eta_lo="0";
  v_fitParameters.at(0).eta_hi="0.8";
  v_fitParameters.at(0).pTbin_lo=3;
  v_fitParameters.at(0).pTbin_hi=4;
  v_fitParameters.at(0).etabin_lo=1;
  v_fitParameters.at(0).etabin_hi=2;
  v_fitParameters.at(0).fitRangeUS_lo=60; v_fitParameters.at(0).fitRangeUS_hi=130; 
  v_fitParameters.at(0).fitRangeSS_lo=20; v_fitParameters.at(0).fitRangeSS_hi=200; 
  v_fitParameters.at(0).parLimitUS_lo[0]=80;  v_fitParameters.at(0).parLimitUS_hi[0]=110;
  v_fitParameters.at(0).parLimitUS_lo[1]=1;  v_fitParameters.at(0).parLimitUS_hi[1]=10;
  v_fitParameters.at(0).parLimitUS_lo[2]=60;  v_fitParameters.at(0).parLimitUS_hi[2]=120;
  v_fitParameters.at(0).parLimitUS_lo[3]=10;  v_fitParameters.at(0).parLimitUS_hi[3]=100;
  v_fitParameters.at(0).parLimitUS_lo[4]=5000; v_fitParameters.at(0).parLimitUS_hi[4]=15000;
  v_fitParameters.at(0).parLimitUS_lo[5]=100;  v_fitParameters.at(0).parLimitUS_hi[5]=1000;
  v_fitParameters.at(0).parLimitSS_lo[0]=80;  v_fitParameters.at(0).parLimitSS_hi[0]=110;
  v_fitParameters.at(0).parLimitSS_lo[1]=5;  v_fitParameters.at(0).parLimitSS_hi[1]=40;
  v_fitParameters.at(0).parLimitSS_lo[2]=30;  v_fitParameters.at(0).parLimitSS_hi[2]=110;
  v_fitParameters.at(0).parLimitSS_lo[3]=20;  v_fitParameters.at(0).parLimitSS_hi[3]=100;
  v_fitParameters.at(0).parLimitSS_lo[4]=5; v_fitParameters.at(0).parLimitSS_hi[4]=20;
  v_fitParameters.at(0).parLimitSS_lo[5]=2;  v_fitParameters.at(0).parLimitSS_hi[5]=10;
  
  v_fitParameters.at(1).pT_lo="40";
  v_fitParameters.at(1).pT_hi="60";
  v_fitParameters.at(1).eta_lo="0";
  v_fitParameters.at(1).eta_hi="0.8";
  v_fitParameters.at(1).pTbin_lo=5;
  v_fitParameters.at(1).pTbin_hi=7;
  v_fitParameters.at(1).etabin_lo=1;
  v_fitParameters.at(1).etabin_hi=2;
  v_fitParameters.at(1).fitRangeUS_lo=60; v_fitParameters.at(1).fitRangeUS_hi=130; 
  v_fitParameters.at(1).fitRangeSS_lo=20; v_fitParameters.at(1).fitRangeSS_hi=200; 
  v_fitParameters.at(1).parLimitUS_lo[0]=80;  v_fitParameters.at(1).parLimitUS_hi[0]=110;
  v_fitParameters.at(1).parLimitUS_lo[1]=1;  v_fitParameters.at(1).parLimitUS_hi[1]=10;
  v_fitParameters.at(1).parLimitUS_lo[2]=60;  v_fitParameters.at(1).parLimitUS_hi[2]=120;
  v_fitParameters.at(1).parLimitUS_lo[3]=10;  v_fitParameters.at(1).parLimitUS_hi[3]=100;
  v_fitParameters.at(1).parLimitUS_lo[4]=15000; v_fitParameters.at(1).parLimitUS_hi[4]=30000;
  v_fitParameters.at(1).parLimitUS_lo[5]=100;  v_fitParameters.at(1).parLimitUS_hi[5]=2000;
  v_fitParameters.at(1).parLimitSS_lo[0]=80;  v_fitParameters.at(1).parLimitSS_hi[0]=110;
  v_fitParameters.at(1).parLimitSS_lo[1]=5;  v_fitParameters.at(1).parLimitSS_hi[1]=40;
  v_fitParameters.at(1).parLimitSS_lo[2]=30;  v_fitParameters.at(1).parLimitSS_hi[2]=170;
  v_fitParameters.at(1).parLimitSS_lo[3]=20;  v_fitParameters.at(1).parLimitSS_hi[3]=200;
  v_fitParameters.at(1).parLimitSS_lo[4]=5; v_fitParameters.at(1).parLimitSS_hi[4]=20;
  v_fitParameters.at(1).parLimitSS_lo[5]=2;  v_fitParameters.at(1).parLimitSS_hi[5]=10;
  
  v_fitParameters.at(2).pT_lo="60";
  v_fitParameters.at(2).pT_hi="80";
  v_fitParameters.at(2).eta_lo="0";
  v_fitParameters.at(2).eta_hi="0.8";
  v_fitParameters.at(2).pTbin_lo=8;
  v_fitParameters.at(2).pTbin_hi=9;
  v_fitParameters.at(2).etabin_lo=1;
  v_fitParameters.at(2).etabin_hi=2;
  v_fitParameters.at(2).fitRangeUS_lo=60; v_fitParameters.at(2).fitRangeUS_hi=130; 
  v_fitParameters.at(2).fitRangeSS_lo=20; v_fitParameters.at(2).fitRangeSS_hi=200; 
  v_fitParameters.at(2).parLimitUS_lo[0]=80;  v_fitParameters.at(2).parLimitUS_hi[0]=110;
  v_fitParameters.at(2).parLimitUS_lo[1]=1;  v_fitParameters.at(2).parLimitUS_hi[1]=10;
  v_fitParameters.at(2).parLimitUS_lo[2]=60;  v_fitParameters.at(2).parLimitUS_hi[2]=120;
  v_fitParameters.at(2).parLimitUS_lo[3]=10;  v_fitParameters.at(2).parLimitUS_hi[3]=100;
  v_fitParameters.at(2).parLimitUS_lo[4]=5000; v_fitParameters.at(2).parLimitUS_hi[4]=10000;
  v_fitParameters.at(2).parLimitUS_lo[5]=100;  v_fitParameters.at(2).parLimitUS_hi[5]=2000;
  v_fitParameters.at(2).parLimitSS_lo[0]=80;  v_fitParameters.at(2).parLimitSS_hi[0]=110;
  v_fitParameters.at(2).parLimitSS_lo[1]=5;  v_fitParameters.at(2).parLimitSS_hi[1]=40;
  v_fitParameters.at(2).parLimitSS_lo[2]=30;  v_fitParameters.at(2).parLimitSS_hi[2]=170;
  v_fitParameters.at(2).parLimitSS_lo[3]=20;  v_fitParameters.at(2).parLimitSS_hi[3]=200;
  v_fitParameters.at(2).parLimitSS_lo[4]=2; v_fitParameters.at(2).parLimitSS_hi[4]=15;
  v_fitParameters.at(2).parLimitSS_lo[5]=0;  v_fitParameters.at(2).parLimitSS_hi[5]=5;
  
  v_fitParameters.at(3).pT_lo="80";
  v_fitParameters.at(3).pT_hi="inf";
  v_fitParameters.at(3).eta_lo="0";
  v_fitParameters.at(3).eta_hi="0.8";
  v_fitParameters.at(3).pTbin_lo=10;
  v_fitParameters.at(3).pTbin_hi=12;
  v_fitParameters.at(3).etabin_lo=1;
  v_fitParameters.at(3).etabin_hi=2;
  v_fitParameters.at(3).fitRangeUS_lo=60; v_fitParameters.at(3).fitRangeUS_hi=130; 
  v_fitParameters.at(3).fitRangeSS_lo=20; v_fitParameters.at(3).fitRangeSS_hi=200; 
  v_fitParameters.at(3).parLimitUS_lo[0]=80;  v_fitParameters.at(3).parLimitUS_hi[0]=110;
  v_fitParameters.at(3).parLimitUS_lo[1]=1;  v_fitParameters.at(3).parLimitUS_hi[1]=10;
  v_fitParameters.at(3).parLimitUS_lo[2]=60;  v_fitParameters.at(3).parLimitUS_hi[2]=120;
  v_fitParameters.at(3).parLimitUS_lo[3]=10;  v_fitParameters.at(3).parLimitUS_hi[3]=100;
  v_fitParameters.at(3).parLimitUS_lo[4]=3000; v_fitParameters.at(3).parLimitUS_hi[4]=7000;
  v_fitParameters.at(3).parLimitUS_lo[5]=100;  v_fitParameters.at(3).parLimitUS_hi[5]=2000;
  v_fitParameters.at(3).parLimitSS_lo[0]=80;  v_fitParameters.at(3).parLimitSS_hi[0]=110;
  v_fitParameters.at(3).parLimitSS_lo[1]=1;  v_fitParameters.at(3).parLimitSS_hi[1]=20;
  v_fitParameters.at(3).parLimitSS_lo[2]=30;  v_fitParameters.at(3).parLimitSS_hi[2]=200;
  v_fitParameters.at(3).parLimitSS_lo[3]=30;  v_fitParameters.at(3).parLimitSS_hi[3]=200;
  v_fitParameters.at(3).parLimitSS_lo[4]=2; v_fitParameters.at(3).parLimitSS_hi[4]=15;
  v_fitParameters.at(3).parLimitSS_lo[5]=2;  v_fitParameters.at(3).parLimitSS_hi[5]=5;
  
  v_fitParameters.at(4).pT_lo="30";
  v_fitParameters.at(4).pT_hi="40";
  v_fitParameters.at(4).eta_lo="0.8";
  v_fitParameters.at(4).eta_hi="1.4442";
  v_fitParameters.at(4).pTbin_lo=3;
  v_fitParameters.at(4).pTbin_hi=4;
  v_fitParameters.at(4).etabin_lo=3;
  v_fitParameters.at(4).etabin_hi=3;
  v_fitParameters.at(4).fitRangeUS_lo=60; v_fitParameters.at(4).fitRangeUS_hi=130; 
  v_fitParameters.at(4).fitRangeSS_lo=20; v_fitParameters.at(4).fitRangeSS_hi=200; 
  v_fitParameters.at(4).parLimitUS_lo[0]=80;  v_fitParameters.at(4).parLimitUS_hi[0]=110;
  v_fitParameters.at(4).parLimitUS_lo[1]=1;  v_fitParameters.at(4).parLimitUS_hi[1]=10;
  v_fitParameters.at(4).parLimitUS_lo[2]=60;  v_fitParameters.at(4).parLimitUS_hi[2]=120;
  v_fitParameters.at(4).parLimitUS_lo[3]=10;  v_fitParameters.at(4).parLimitUS_hi[3]=100;
  v_fitParameters.at(4).parLimitUS_lo[4]=4000; v_fitParameters.at(4).parLimitUS_hi[4]=15000;
  v_fitParameters.at(4).parLimitUS_lo[5]=100;  v_fitParameters.at(4).parLimitUS_hi[5]=1000;
  v_fitParameters.at(4).parLimitSS_lo[0]=80;  v_fitParameters.at(4).parLimitSS_hi[0]=110;
  v_fitParameters.at(4).parLimitSS_lo[1]=1;  v_fitParameters.at(4).parLimitSS_hi[1]=15;
  v_fitParameters.at(4).parLimitSS_lo[2]=30;  v_fitParameters.at(4).parLimitSS_hi[2]=200;
  v_fitParameters.at(4).parLimitSS_lo[3]=50;  v_fitParameters.at(4).parLimitSS_hi[3]=200;
  v_fitParameters.at(4).parLimitSS_lo[4]=5; v_fitParameters.at(4).parLimitSS_hi[4]=20;
  v_fitParameters.at(4).parLimitSS_lo[5]=0;  v_fitParameters.at(4).parLimitSS_hi[5]=10;
  
  v_fitParameters.at(5).pT_lo="40";
  v_fitParameters.at(5).pT_hi="60";
  v_fitParameters.at(5).eta_lo="0.8";
  v_fitParameters.at(5).eta_hi="1.4442";
  v_fitParameters.at(5).pTbin_lo=5;
  v_fitParameters.at(5).pTbin_hi=7;
  v_fitParameters.at(5).etabin_lo=3;
  v_fitParameters.at(5).etabin_hi=3;
  v_fitParameters.at(5).fitRangeUS_lo=60; v_fitParameters.at(5).fitRangeUS_hi=130; 
  v_fitParameters.at(5).fitRangeSS_lo=20; v_fitParameters.at(5).fitRangeSS_hi=200; 
  v_fitParameters.at(5).parLimitUS_lo[0]=80;  v_fitParameters.at(5).parLimitUS_hi[0]=110;
  v_fitParameters.at(5).parLimitUS_lo[1]=1;  v_fitParameters.at(5).parLimitUS_hi[1]=10;
  v_fitParameters.at(5).parLimitUS_lo[2]=60;  v_fitParameters.at(5).parLimitUS_hi[2]=120;
  v_fitParameters.at(5).parLimitUS_lo[3]=10;  v_fitParameters.at(5).parLimitUS_hi[3]=100;
  v_fitParameters.at(5).parLimitUS_lo[4]=5000; v_fitParameters.at(5).parLimitUS_hi[4]=10000;
  v_fitParameters.at(5).parLimitUS_lo[5]=100;  v_fitParameters.at(5).parLimitUS_hi[5]=2000;
  v_fitParameters.at(5).parLimitSS_lo[0]=80;  v_fitParameters.at(5).parLimitSS_hi[0]=110;
  v_fitParameters.at(5).parLimitSS_lo[1]=5;  v_fitParameters.at(5).parLimitSS_hi[1]=40;
  v_fitParameters.at(5).parLimitSS_lo[2]=30;  v_fitParameters.at(5).parLimitSS_hi[2]=170;
  v_fitParameters.at(5).parLimitSS_lo[3]=20;  v_fitParameters.at(5).parLimitSS_hi[3]=200;
  v_fitParameters.at(5).parLimitSS_lo[4]=5; v_fitParameters.at(5).parLimitSS_hi[4]=20;
  v_fitParameters.at(5).parLimitSS_lo[5]=2;  v_fitParameters.at(5).parLimitSS_hi[5]=10;
  
  v_fitParameters.at(6).pT_lo="60";
  v_fitParameters.at(6).pT_hi="80";
  v_fitParameters.at(6).eta_lo="0.8";
  v_fitParameters.at(6).eta_hi="1.4442";
  v_fitParameters.at(6).pTbin_lo=8;
  v_fitParameters.at(6).pTbin_hi=9;
  v_fitParameters.at(6).etabin_lo=3;
  v_fitParameters.at(6).etabin_hi=3;
  v_fitParameters.at(6).fitRangeUS_lo=60; v_fitParameters.at(6).fitRangeUS_hi=130; 
  v_fitParameters.at(6).fitRangeSS_lo=20; v_fitParameters.at(6).fitRangeSS_hi=200; 
  v_fitParameters.at(6).parLimitUS_lo[0]=80;  v_fitParameters.at(6).parLimitUS_hi[0]=110;
  v_fitParameters.at(6).parLimitUS_lo[1]=1;  v_fitParameters.at(6).parLimitUS_hi[1]=10;
  v_fitParameters.at(6).parLimitUS_lo[2]=60;  v_fitParameters.at(6).parLimitUS_hi[2]=120;
  v_fitParameters.at(6).parLimitUS_lo[3]=10;  v_fitParameters.at(6).parLimitUS_hi[3]=100;
  v_fitParameters.at(6).parLimitUS_lo[4]=1000; v_fitParameters.at(6).parLimitUS_hi[4]=3000;
  v_fitParameters.at(6).parLimitUS_lo[5]=100;  v_fitParameters.at(6).parLimitUS_hi[5]=2000;
  v_fitParameters.at(6).parLimitSS_lo[0]=80;  v_fitParameters.at(6).parLimitSS_hi[0]=110;
  v_fitParameters.at(6).parLimitSS_lo[1]=1;  v_fitParameters.at(6).parLimitSS_hi[1]=30;
  v_fitParameters.at(6).parLimitSS_lo[2]=20;  v_fitParameters.at(6).parLimitSS_hi[2]=200;
  v_fitParameters.at(6).parLimitSS_lo[3]=20;  v_fitParameters.at(6).parLimitSS_hi[3]=200;
  v_fitParameters.at(6).parLimitSS_lo[4]=2; v_fitParameters.at(6).parLimitSS_hi[4]=20;
  v_fitParameters.at(6).parLimitSS_lo[5]=0;  v_fitParameters.at(6).parLimitSS_hi[5]=5;
  
  v_fitParameters.at(7).pT_lo="80";
  v_fitParameters.at(7).pT_hi="inf";
  v_fitParameters.at(7).eta_lo="0.8";
  v_fitParameters.at(7).eta_hi="1.4442";
  v_fitParameters.at(7).pTbin_lo=10;
  v_fitParameters.at(7).pTbin_hi=12;
  v_fitParameters.at(7).etabin_lo=3;
  v_fitParameters.at(7).etabin_hi=3;
  v_fitParameters.at(7).fitRangeUS_lo=40; v_fitParameters.at(7).fitRangeUS_hi=160; 
  v_fitParameters.at(7).fitRangeSS_lo=20; v_fitParameters.at(7).fitRangeSS_hi=200; 
  v_fitParameters.at(7).parLimitUS_lo[0]=80;  v_fitParameters.at(7).parLimitUS_hi[0]=110;
  v_fitParameters.at(7).parLimitUS_lo[1]=1;  v_fitParameters.at(7).parLimitUS_hi[1]=10;
  v_fitParameters.at(7).parLimitUS_lo[2]=60;  v_fitParameters.at(7).parLimitUS_hi[2]=120;
  v_fitParameters.at(7).parLimitUS_lo[3]=10;  v_fitParameters.at(7).parLimitUS_hi[3]=100;
  v_fitParameters.at(7).parLimitUS_lo[4]=1000; v_fitParameters.at(7).parLimitUS_hi[4]=3000;
  v_fitParameters.at(7).parLimitUS_lo[5]=100;  v_fitParameters.at(7).parLimitUS_hi[5]=2000;
  v_fitParameters.at(7).parLimitSS_lo[0]=80;  v_fitParameters.at(7).parLimitSS_hi[0]=110;
  v_fitParameters.at(7).parLimitSS_lo[1]=1;  v_fitParameters.at(7).parLimitSS_hi[1]=30;
  v_fitParameters.at(7).parLimitSS_lo[2]=20;  v_fitParameters.at(7).parLimitSS_hi[2]=200;
  v_fitParameters.at(7).parLimitSS_lo[3]=20;  v_fitParameters.at(7).parLimitSS_hi[3]=200;
  v_fitParameters.at(7).parLimitSS_lo[4]=2; v_fitParameters.at(7).parLimitSS_hi[4]=15;
  v_fitParameters.at(7).parLimitSS_lo[5]=2;  v_fitParameters.at(7).parLimitSS_hi[5]=5;
  
  v_fitParameters.at(8).pT_lo="30";
  v_fitParameters.at(8).pT_hi="40";
  v_fitParameters.at(8).eta_lo="1.4442";
  v_fitParameters.at(8).eta_hi="2.4";
  v_fitParameters.at(8).pTbin_lo=3;
  v_fitParameters.at(8).pTbin_hi=4;
  v_fitParameters.at(8).etabin_lo=4;
  v_fitParameters.at(8).etabin_hi=5;
  v_fitParameters.at(8).fitRangeUS_lo=60; v_fitParameters.at(8).fitRangeUS_hi=130; 
  v_fitParameters.at(8).fitRangeSS_lo=20; v_fitParameters.at(8).fitRangeSS_hi=200; 
  v_fitParameters.at(8).parLimitUS_lo[0]=80;  v_fitParameters.at(8).parLimitUS_hi[0]=110;
  v_fitParameters.at(8).parLimitUS_lo[1]=1;  v_fitParameters.at(8).parLimitUS_hi[1]=10;
  v_fitParameters.at(8).parLimitUS_lo[2]=60;  v_fitParameters.at(8).parLimitUS_hi[2]=120;
  v_fitParameters.at(8).parLimitUS_lo[3]=10;  v_fitParameters.at(8).parLimitUS_hi[3]=100;
  v_fitParameters.at(8).parLimitUS_lo[4]=2000; v_fitParameters.at(8).parLimitUS_hi[4]=4000;
  v_fitParameters.at(8).parLimitUS_lo[5]=100;  v_fitParameters.at(8).parLimitUS_hi[5]=1000;
  v_fitParameters.at(8).parLimitSS_lo[0]=80;  v_fitParameters.at(8).parLimitSS_hi[0]=110;
  v_fitParameters.at(8).parLimitSS_lo[1]=5;  v_fitParameters.at(8).parLimitSS_hi[1]=40;
  v_fitParameters.at(8).parLimitSS_lo[2]=30;  v_fitParameters.at(8).parLimitSS_hi[2]=110;
  v_fitParameters.at(8).parLimitSS_lo[3]=20;  v_fitParameters.at(8).parLimitSS_hi[3]=100;
  v_fitParameters.at(8).parLimitSS_lo[4]=5; v_fitParameters.at(8).parLimitSS_hi[4]=25;
  v_fitParameters.at(8).parLimitSS_lo[5]=0;  v_fitParameters.at(8).parLimitSS_hi[5]=10;
  
  v_fitParameters.at(9).pT_lo="40";
  v_fitParameters.at(9).pT_hi="60";
  v_fitParameters.at(9).eta_lo="1.4442";
  v_fitParameters.at(9).eta_hi="2.4";
  v_fitParameters.at(9).pTbin_lo=5;
  v_fitParameters.at(9).pTbin_hi=7;
  v_fitParameters.at(9).etabin_lo=4;
  v_fitParameters.at(9).etabin_hi=5;
  v_fitParameters.at(9).fitRangeUS_lo=50; v_fitParameters.at(9).fitRangeUS_hi=130; 
  v_fitParameters.at(9).fitRangeSS_lo=20; v_fitParameters.at(9).fitRangeSS_hi=200; 
  v_fitParameters.at(9).parLimitUS_lo[0]=80;  v_fitParameters.at(9).parLimitUS_hi[0]=110;
  v_fitParameters.at(9).parLimitUS_lo[1]=1;  v_fitParameters.at(9).parLimitUS_hi[1]=10;
  v_fitParameters.at(9).parLimitUS_lo[2]=60;  v_fitParameters.at(9).parLimitUS_hi[2]=120;
  v_fitParameters.at(9).parLimitUS_lo[3]=10;  v_fitParameters.at(9).parLimitUS_hi[3]=100;
  v_fitParameters.at(9).parLimitUS_lo[4]=2000; v_fitParameters.at(9).parLimitUS_hi[4]=5000;
  v_fitParameters.at(9).parLimitUS_lo[5]=100;  v_fitParameters.at(9).parLimitUS_hi[5]=2000;
  v_fitParameters.at(9).parLimitSS_lo[0]=80;  v_fitParameters.at(9).parLimitSS_hi[0]=110;
  v_fitParameters.at(9).parLimitSS_lo[1]=5;  v_fitParameters.at(9).parLimitSS_hi[1]=40;
  v_fitParameters.at(9).parLimitSS_lo[2]=30;  v_fitParameters.at(9).parLimitSS_hi[2]=170;
  v_fitParameters.at(9).parLimitSS_lo[3]=20;  v_fitParameters.at(9).parLimitSS_hi[3]=200;
  v_fitParameters.at(9).parLimitSS_lo[4]=5; v_fitParameters.at(9).parLimitSS_hi[4]=20;
  v_fitParameters.at(9).parLimitSS_lo[5]=2;  v_fitParameters.at(9).parLimitSS_hi[5]=10;
  
  v_fitParameters.at(10).pT_lo="60";
  v_fitParameters.at(10).pT_hi="80";
  v_fitParameters.at(10).eta_lo="1.4442";
  v_fitParameters.at(10).eta_hi="2.4";
  v_fitParameters.at(10).pTbin_lo=8;
  v_fitParameters.at(10).pTbin_hi=9;
  v_fitParameters.at(10).etabin_lo=4;
  v_fitParameters.at(10).etabin_hi=5;
  v_fitParameters.at(10).fitRangeUS_lo=40; v_fitParameters.at(10).fitRangeUS_hi=150; 
  v_fitParameters.at(10).fitRangeSS_lo=20; v_fitParameters.at(10).fitRangeSS_hi=200; 
  v_fitParameters.at(10).parLimitUS_lo[0]=80;  v_fitParameters.at(10).parLimitUS_hi[0]=110;
  v_fitParameters.at(10).parLimitUS_lo[1]=1;  v_fitParameters.at(10).parLimitUS_hi[1]=10;
  v_fitParameters.at(10).parLimitUS_lo[2]=60;  v_fitParameters.at(10).parLimitUS_hi[2]=120;
  v_fitParameters.at(10).parLimitUS_lo[3]=10;  v_fitParameters.at(10).parLimitUS_hi[3]=100;
  v_fitParameters.at(10).parLimitUS_lo[4]=600; v_fitParameters.at(10).parLimitUS_hi[4]=1200;
  v_fitParameters.at(10).parLimitUS_lo[5]=100;  v_fitParameters.at(10).parLimitUS_hi[5]=2000;
  v_fitParameters.at(10).parLimitSS_lo[0]=80;  v_fitParameters.at(10).parLimitSS_hi[0]=120;
  v_fitParameters.at(10).parLimitSS_lo[1]=1;  v_fitParameters.at(10).parLimitSS_hi[1]=25;
  v_fitParameters.at(10).parLimitSS_lo[2]=20;  v_fitParameters.at(10).parLimitSS_hi[2]=200;
  v_fitParameters.at(10).parLimitSS_lo[3]=20;  v_fitParameters.at(10).parLimitSS_hi[3]=200;
  v_fitParameters.at(10).parLimitSS_lo[4]=2; v_fitParameters.at(10).parLimitSS_hi[4]=20;
  v_fitParameters.at(10).parLimitSS_lo[5]=1;  v_fitParameters.at(10).parLimitSS_hi[5]=5;
  
  v_fitParameters.at(11).pT_lo="80";
  v_fitParameters.at(11).pT_hi="inf";
  v_fitParameters.at(11).eta_lo="1.4442";
  v_fitParameters.at(11).eta_hi="2.4";
  v_fitParameters.at(11).pTbin_lo=10;
  v_fitParameters.at(11).pTbin_hi=12;
  v_fitParameters.at(11).etabin_lo=4;
  v_fitParameters.at(11).etabin_hi=5;
  v_fitParameters.at(11).fitRangeUS_lo=20; v_fitParameters.at(11).fitRangeUS_hi=200; 
  v_fitParameters.at(11).fitRangeSS_lo=20; v_fitParameters.at(11).fitRangeSS_hi=200; 
  v_fitParameters.at(11).parLimitUS_lo[0]=80;  v_fitParameters.at(11).parLimitUS_hi[0]=110;
  v_fitParameters.at(11).parLimitUS_lo[1]=10;  v_fitParameters.at(11).parLimitUS_hi[1]=30;
  v_fitParameters.at(11).parLimitUS_lo[2]=90;  v_fitParameters.at(11).parLimitUS_hi[2]=150;
  v_fitParameters.at(11).parLimitUS_lo[3]=40;  v_fitParameters.at(11).parLimitUS_hi[3]=80;
  v_fitParameters.at(11).parLimitUS_lo[4]=100; v_fitParameters.at(11).parLimitUS_hi[4]=500;
  v_fitParameters.at(11).parLimitUS_lo[5]=10;  v_fitParameters.at(11).parLimitUS_hi[5]=200;
  v_fitParameters.at(11).parLimitSS_lo[0]=80;  v_fitParameters.at(11).parLimitSS_hi[0]=100;
  v_fitParameters.at(11).parLimitSS_lo[1]=5;  v_fitParameters.at(11).parLimitSS_hi[1]=40;
  v_fitParameters.at(11).parLimitSS_lo[2]=90;  v_fitParameters.at(11).parLimitSS_hi[2]=150;
  v_fitParameters.at(11).parLimitSS_lo[3]=50;  v_fitParameters.at(11).parLimitSS_hi[3]=100;
  v_fitParameters.at(11).parLimitSS_lo[4]=5; v_fitParameters.at(11).parLimitSS_hi[4]=15;
  v_fitParameters.at(11).parLimitSS_lo[5]=2;  v_fitParameters.at(11).parLimitSS_hi[5]=10;
  
  gStyle->SetOptStat(0);
  
  std::ofstream outfile("qFRDashboard/index.html");
  outfile<<"<html>"<<std::endl;
  outfile<<"<head>"<<std::endl;
  // outfile<<" <link rel=\"StyleSheet\" type=\"text/css\" href=\"PRBDashboardStyle.css\" />"<<std::endl;
  outfile<<"</head>"<<std::endl;
  outfile<<"<body>"<<std::endl;
  outfile<<"<h1 align='center'> VLQ Charge Flip Rate Measurements </h1>"<<std::endl;
  outfile<<"<table border='1'>"<<std::endl;
  for (unsigned int i=0; i<12; ++i)
  {
    fit(h_PromptRate_ee_TT_OS, h_PromptRate_ee_TT_SS, &(v_fitParameters.at(i)), outfile);
  }
  outfile<<"</body>"<<std::endl;
  outfile<<"</html>"<<std::endl;
  outfile.close(); 

}
  

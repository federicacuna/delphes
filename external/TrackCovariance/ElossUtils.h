  std::vector<double> avampv{0.00314338, 0.00164584, 0.00112617, 0.000887539, 0.000679093, 0.000630974, 0.00059563, 0.000573898, 0.000556587,
    0.000547491, 0.000539621, 0.000521902, 0.000527015, 0.000530534, 0.00058391, 0.000647633, 0.000680251, 0.000675453};
  std::vector<double> avbmpv{-0.00114136, -0.0006393, -0.00045854, -0.000374344, -0.000291553, -0.000274131, -0.000260206, -0.000251556, -0.000244366, -0.00024092, -0.000237356,
    -0.000227837, -0.000232123, -0.00022986, -0.000225664, -0.000238779, -0.000259258, -0.000249283};
  std::vector<double> avasgm{0.000205066, 0.000114538, 8.00996e-05, 4.77652e-05, 4.20902e-05, 4.09978e-05, 4.02946e-05, 3.82017e-05, 3.77149e-05, 3.72186e-05, 3.56652e-05,
    3.37644e-05, 3.61425e-05, 3.59534e-05, 3.39436e-05, 3.32454e-05, 3.35233e-05, 3.3515e-05};
  std::vector<double> avbsgm{0.000118624, 6.91577e-05, 5.35544e-05, 3.97809e-05, 4.32185e-05, 4.11169e-05, 3.88939e-05, 3.92737e-05, 3.80544e-05, 3.77401e-05, 3.67337e-05, 3.89669e-05, 3.66332e-05, 3.76185e-05, 4.03574e-05,
    4.9421e-05, 5.19735e-05, 5.32526e-05};
  std::vector<double> bg_sim{0.405123, 0.607685, 0.810246, 1.01281, 1.41793, 1.62049, 1.82305, 2.02562, 2.22818, 2.43074, 2.6333, 3.64611, 4.05123,
    6.07685, 20.2562, 81.0246, 141.793, 202.562};

      std::vector<double> xvaluempv{0.191898,0.21322,0.303842,0.319829,0.364611,
      0.405123,0.426439,0.533049,0.607685,0.746269,0.810246,0.959488,1.01281,1.0661,
      1.17271,1.27932,1.41793,1.62049,1.70576,1.82305,1.89286,1.91898,2.02562,
      2.1322,2.22818,2.43074,2.6333,2.8393,3.24099,3.64611,3.78573,3.84867,
      4.05123,4.73216,5.33049,5.67859,5.73189,6.07685,6.62502,7.57146,8.51789,
      9.31432,9.46432,10.4108,10.661,11.3572,11.4638,12.3036,13.6132,15.1429,
      17.0358,17.9822,18.9286,20.2562,21.322,21.4946,28.393,71.6486,81.0246,
      94.6432,106.61,141.793,189.286,202.562,293.542,352.251,378.573,391.389,
      587.084,662.502,782.779,946.432,978.474,1174.17,1369.86,1761.25,1956.95,2152.64,
      2348.34,2544.03,3131.12,3522.5,3718.2,5870.84,39138.9,78277.9,136986,195695};
      
    std::vector<double> yvaluempv{0.00949555, 0.00666235, 0.00347429, 0.00305658, 0.00243737, 0.00200202, 0.00180494, 0.00123294, 0.00100654, 0.00078075, 0.000667633, 0.000540799, 0.000513195, 0.000484929, 0.000443421,
      0.000412801, 0.00038754, 0.000356843, 0.000344828, 0.000335424, 0.000330722, 0.000327767, 0.000322342, 0.000316176,
      0.000312221, 0.000306571, 0.000302265, 0.000298383, 0.000295586, 0.000294065, 0.000294321, 0.000294783, 0.000294892,
      0.000295182, 0.000298384, 0.000299649, 0.000298202, 0.000300674, 0.000302966, 0.000306942, 0.000311452, 0.000312653,
      0.000315662, 0.000320251, 0.000310405, 0.000322513, 0.000322837, 0.000325249, 0.000331152, 0.00033713, 0.00034315,
      0.000350659, 0.000352483, 0.000358246, 0.00034913, 0.00036068, 0.000370566, 0.000402627, 0.000408854, 0.000412043,
      0.000414483, 0.000420993, 0.000425194, 0.00042617, 0.000430491, 0.000432422, 0.00043299, 0.000432992, 0.000435823,
      0.000435572, 0.00043622, 0.000436652, 0.00043675, 0.000437264, 0.000437809, 0.000438576, 0.000438447,
      0.000438168, 0.000437641, 0.00043821, 0.000437661, 0.000438526, 0.000440099, 0.000438987, 0.00043622, 0.000439298, 0.000439135, 0.000438865};

std::vector<double>xvaluesgm{0.303842,0.319829,0.364611,0.405123,0.426439,0.533049,0.607685,0.810246,
    0.959488,1.01281,1.0661,1.17271,1.27932,1.28968,1.41793,1.43297,1.62049,1.70576,1.82305,1.89286,
    1.91898,2.02562,2.1322,2.22818,2.43074,2.6333,2.8393,3.24099,3.64611,3.78573,3.84867,4.05123,4.73216,
    5.67859,5.73189,6.07685,6.62502,7.57146,8.51789,9.31432,9.46432,10.4108,11.3572,12.3036,15.1429,17.0358,
    81.0246,94.6432,141.793,189.286,202.562,293.542,352.251,378.573,391.389,587.084,662.502,782.779,946.432,978.474,1174.17,
    1369.86,1761.25,1956.95,2152.64,2348.34,2544.03,3131.12,3522.5,3718.2,5870.84,39138.9,78277.9,136986,195695};
std::vector<double>yvaluesgm{0.000540235,0.000446181,0.000385976,0.00032369,0.000293805,0.000220726,0.000183696,0.000133654,
    0.0001213,0.000110032,0.000112333,0.000106138,0.000100757,9.0307e-05,8.75461e-05,8.54347e-05,8.53087e-05,8.91189e-05,
    8.21147e-05,8.13459e-05,8.61202e-05,7.91885e-05,8.39838e-05,7.74754e-05,7.57693e-05,7.49587e-05,7.39794e-05,7.30484e-05,
    7.23989e-05,7.25857e-05,7.27313e-05,7.27757e-05,7.25107e-05,7.32769e-05,7.20314e-05,7.35719e-05,7.38514e-05,7.44724e-05,
    7.54768e-05,7.42085e-05,7.60001e-05,7.66292e-05,7.70694e-05,7.7171e-05,7.82256e-05,7.85329e-05,8.26664e-05,8.3338e-05,
    8.54968e-05,8.61602e-05,8.67676e-05,8.79901e-05,8.88283e-05,8.86416e-05,8.84529e-05,8.95503e-05,9.01781e-05,9.05446e-05,
    9.0461e-05,9.00113e-05,9.03967e-05,9.07454e-05,9.06081e-05,9.09096e-05,9.09226e-05,9.05433e-05,9.12025e-05,9.08368e-05,
    9.11519e-05,9.14264e-05,9.12324e-05,9.05446e-05,9.10745e-05,9.15355e-05,9.11858e-05};





double funcBB(double *xp, double *par)
{

  double x = xp[0] - par[4];
  double func = 0;
  if(x > 0)
  {
    double xsq = x * x;
    double delta = 0.0;
    double xd = TMath::Log10(x);
    double xLog = par[1] * xsq * xsq / (1.0 + par[2] * TMath::Sqrt(1.0 + xsq) + par[2] * par[2] * 0.25);

    if(xLog > 1)
    {
      if(xd >= par[6])
        delta = par[7] * xd - par[8];
      else if(xd < par[6] && xd >= par[5])
        delta = par[7] * xd - par[8] + par[9] * pow((par[6] - xd), par[10]);
      else if(xd < par[5])
      {
        delta = 0.0;
      }
      func = 0.5 * par[0] * (1.0 + xsq) / xsq * (TMath::Log(xLog) - delta) - par[0] + par[3];
    }
  }

  return func;
}
double ExtrapolTD(std::vector<double>xvalue,std::vector<double>yvalue,double bg,unsigned int npt=2){
//    cout<<"extrapolation "<<endl;
    float meand=0.0;
    float meant=0.0;
    float stt=0.0;
    float std=0.0; //varianze
    double tmpresult=0.0;
    npt=(xvalue.size()>npt)?(npt):xvalue.size();
    int start=(xvalue.size()-npt)>0?(xvalue.size()-npt):0;
    for(unsigned int i=start;i<xvalue.size();++i){
        meand+=yvalue[i];
        meant+=xvalue[i];
        stt+=xvalue[i]*xvalue[i];
        std+=xvalue[i]*yvalue[i];
//        cout<<"x "<<xvalue[i]<< " y "<<yvalue[i]<<" std "<<std<<" stt "<<stt<<endl;
    }
//    cout<<" out " <<" std "<<std<<" stt "<<stt<<endl;
    meand/=(float)npt;
    meant/=(float)npt;
    std=std/(float)npt-meand*meant;
    stt=stt/(float)npt-meant*meant;
    tmpresult=(std/stt)*(bg-meant)+meand;
//    cout<<"meand "<< meand<<" meant "<<meant<<" std "<<std<<" stt "<<stt<<" tmpresult "<<tmpresult<<endl;

//    cout<<"extrapolation end"<<endl;
    if(isnan(tmpresult)||isinf(tmpresult)){
        std::cout<<"check ExtrapolTD"<<std::endl;
        std::cout<<"meand "<< meand<<" meant "<<meant<<" std "<<std<<" stt "<<stt<<std::endl;
    }
    return tmpresult;
}

double lincoeff(double bg, std::vector<double> xvalue, std::vector<double> yvalue)
{
    // std::cout<<"interpolation "<<std::endl;

    ROOT::Math::Interpolator *itp;

    double tmpval, tmpbg1, tmpbg2, tmpa1, tmpa2, tmpbg1_1, tmpbg2_1, tmpa1_1, tmpa2_1;
    // if(itp != nullptr)
    // {
        // std::cout<<"interpolation 1"<<std::endl;

        itp = new ROOT::Math::Interpolator(xvalue, yvalue, ROOT::Math::Interpolation::kLINEAR);
    // }
    bool dx = false;
    bool sx = false;

    if(bg > xvalue[xvalue.size()-1])
    {
        dx = true;
        //      std::cout << " out range dx " << bg<<" bg_sim[bg_sim.size()] "<<xvalue[xvalue.size()-1] <<std::endl;
    }
    if(bg < xvalue[0])
    {
        sx = true;
    }

    if(dx)
    {
        tmpbg1 = xvalue[xvalue.size()-2];
        tmpbg2 = xvalue[xvalue.size()-1];
        tmpa1 = yvalue[xvalue.size()-2];
        tmpa2 = yvalue[xvalue.size()-1];
        tmpval = ExtrapolTD(xvalue, yvalue, bg, 10);//(bg - tmpbg1) / (tmpbg2 - tmpbg1) * (tmpa2 - tmpa1) + tmpa1;
        //    std::cout << " out range dx " << bg << " tmpval " << tmpval << std::endl;
    }
    else if(sx)
    {
        tmpbg1_1 = xvalue[0];
        tmpbg2_1 = xvalue[1];
        tmpa1_1 = yvalue[0];
        tmpa2_1 = yvalue[1];
        tmpval = (bg - tmpbg1) / (tmpbg2 - tmpbg1) * (tmpa2 - tmpa1) + tmpa1;
        //    std::cout << " out range sx " << bg <<" tmpval " << tmpval << std::endl;
    }
    else
    {
        tmpval = itp->Eval(bg);
        //    std::cout << " in range " << bg << " tmpval " << tmpval << std::endl;
    }
    return tmpval;
}

Double_t EnergyLoss(Double_t begam)
{
  //
  //energy loss from B&B fit formula
  //
  TF1 *fit_bb = new TF1("fit_bb", funcBB, 0.2, 100, 11);
  fit_bb->SetParameters(4.123e-5, 3.171e5, 0.003083, 1.354e-7, 0.1055, 0.0078, 1.864, 7.035, 4.961, 3.547, 0.6214);
  const double *error = new double[11]{2.345e-8, 1604, 9.4823e-5, 0.0001436, 2.261e-5, 0.001616, 0.009123, 0.01809, 0.01187, 0.002589};
  fit_bb->SetParErrors(error);
  // fit_bb->Draw();
//   std::cout<<fit_bb->Eval(begam)<<std::endl;
  return fit_bb->Eval(begam);
}
Double_t EnergyLossSigma(Double_t begam)
{
  //
  //energy loss sigma from B&B fit formula
  //
  TF1 *fbbsigma = new TF1("fbbsigma", funcBB, 0.1, 100, 11);
  fbbsigma->SetParameters(8.675e-6, 1.934e+5, 0.0653, 1.171e-5, 0.05919, 0.2695, 2.028, 6.637, 4.284, 2.024, 0.5669);
  const double *error = new double[11]{2.755e-08, 7672, 0.003599, 2.069e-7, 0.0009667, 7.834e-06, 0.0001607, 0.01066, 0.04196, 0.03019, 0.02163};
  fbbsigma->SetParErrors(error);
  return fbbsigma->Eval(begam);
}

double coeffa_mpv(float bg)
{
  TF1 *fbb = new TF1("fbb", funcBB, 0.1, 100, 11);
  fbb->SetParameters(4.05e-5, 2.397e11, 2.403, -5.308e-5, 0.02397, 78.38, 1.257, 2.533, 3.34, -3.851e5, 0.6213);
  //    const double *error=new double[11]{2.345e-8,1604,9.4823-5,0.0001436,2.261e-5,0.001616,0.009123,0.01809,0.01187,0.002589};
  //    fbb->SetParErrors(error);
  return fbb->Eval(bg);
}

double coeffb_mpv(float bg)
{
  TF1 *plateau = new TF1("plateau", "[1]-([1]-[0])*TMath::Exp(-[2]*TMath::Log10(x))");
  plateau->SetParameters(-0.0003928, -0.0002381, 2.051);
  return plateau->Eval(bg);
}

Double_t TruncatedMean(std::vector<Double_t> elosses, Double_t truncFrac)
{
  Int_t new_size = Int_t(elosses.size() * (1 - truncFrac));

  // remove outliers and re-compute mean
  elosses.resize(new_size);
  return accumulate(elosses.begin(), elosses.end(), 0.0) / elosses.size();
}
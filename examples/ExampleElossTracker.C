#include "../external/TrackCovariance/SolGeom.h"
#include "../external/TrackCovariance/SolTrack.h"
#include "../external/TrackCovariance/TrkUtil.h"

#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TVectorD.h>
#include <iostream>
#include <sstream>
#include <vector>

void ExampleElossTracker(Int_t Nev = 1000, Double_t ang = 90.)
{
  gStyle->SetOptStat(111111);
  gSystem->Load("libDelphes.so");
  //
  // Constants
  Double_t angRad = ang * TMath::Pi() / 180.; // To radians
  Double_t KpMass = 0.493677; // Charged Kaon mass
  Double_t PiMass = 0.13957018; // Charged Pion mass
  Double_t Bz = 2.0; // Field
  TVector3 x(0.0, 0.0, 0.0); // Track origin
  //
  Double_t pmin = 0.01; // Momentum range
  Double_t pmax = 100.;
  //
  // Setup chamber
  //
  TrkUtil *TU = new TrkUtil(Bz);
  TString Geomdata = "/afs/cern.ch/user/f/fcuna/Delphes_dEdx/dchgeometry.txt";
  TU->setGeomfile(Geomdata);
  //cout << "TU done" << endl;
  Double_t Rmin = 0.35;
  Double_t Rmax = 2.0;
  Double_t Zmin = -2.0;
  Double_t Zmax = 2.0;
  TU->SetDchBoundaries(Rmin, Rmax, Zmin, Zmax);
  SolGeom *fG = new SolGeom();
  fG->ReadDCH(Geomdata);
  // cout<< "some outpoot "<<tk->GetString("DetectorGeometry", "")<<endl;
  // cout << " B " << fG->B() << " Nl " << fG->Nl() << " " << fG->GetZminPos() << " " << fG->GetRmin() << endl;

  char title[50];

  sprintf(title, "Helium 90 - Isobutane 10"); // He-Isobutane

  //
  // Histograms
  //
  TH2F *h_kaon = new TH2F("h_kaon", "Energy loss distribution", 1000, 0.1, 100., 100, 0, 0.005);
  TH2F *h_pion = new TH2F("h_pion", "Energy loss distribution", 1000, 0.1, 100., 100, 0, 0.005);
  TH2F *h_piondedx = new TH2F("h_piondedx", "Energy loss distribution dedx", 1000, 0.1, 100., 100, 0, 0.005);
  TH1F *hNly = new TH1F("hNly", "hNly", 150, 0, 150);
  TH1F *htlen = new TH1F("htlen", "htlen", 200, 0, 200);
  TH1F *hcelllen = new TH1F("hcelllen", "hcelllen", 200, 0, 2);
  TH2F *htlen_nhit = new TH2F("htlen_nhit", "htlen_nhit", 200, 0, 200, 150, 0, 150);
  TH2F *htlen_p = new TH2F("htlen_p", "htlen_p", 200, 0, 200, 100, 0, 100);

  TF1 *fgausP = new TF1("fgausP", "gaus");
  TF1 *fgausK = new TF1("fgausK", "gaus");

  //cout << "histograms booked" << endl;
  //
  // Fill plots
  //
  float d = 0.0;

  // for(Int_t n = 0; n < Nev; n++)
  // {
  //   Double_t R = gRandom->Rndm();
  //   Double_t pval = pmin + R * (pmax - pmin);
  //   Double_t px = pval * TMath::Sin(angRad);
  //   Double_t pz = pval * TMath::Cos(angRad);
  //   Double_t bg = pval / KpMass;
  //   cout<<"bg val "<< bg<<" pval "<< pval<< " kpmass "<<KpMass<<endl;
  //   Double_t bgpi = pval / PiMass;
  //   TVector3 p(px, 0.0, pz);
  //   SolTrack TU2(x, p, fG);
  //   hNly->Fill(TU2.nmHit());
  //   // cout << "nmHit " << TU2.nmHit() << " nHit " << TU2.nHit() << endl;
  //   Double_t Q = 1.0;
  //   TVectorD Par = TrkUtil::XPtoPar(x, p, Q, Bz);
  //   //cout << "Parameters done. pval = "<< pval << endl;
  //   //
  //   hcelllen->Fill(TU->TrkLen(Par)* 1e2 / TU2.nmHit());
  //   htlen_nhit->Fill(TU->TrkLen(Par)* 1e2, TU2.nmHit());
  //   htlen->Fill(TU->TrkLen(Par) * 1e2);
  //   htlen_p->Fill(TU->TrkLen(Par) * 1e2,TMath::Sqrt(px*px+pz*pz));

  //   Double_t Eloss = 0.0;
  //   Double_t dedx = 0.0;

  //   if(TU->IonEnergyLoss(Eloss, KpMass, Par, dedx))
  //   {

  //     // cout << "Kaon extracted. pval = "<< pval << endl;
  //     h_kaon->Fill(pval, Eloss);
  //     // cout << "Kaon filled. Eloss = "<<Eloss << endl;
  //   }
  //   if(TU->IonEnergyLoss(Eloss, PiMass, Par, dedx))
  //   {
  //     // d = (TU->TrkLen(Par) / TU2.nHit())*1e2;
  //     // cout << "ev. "<<n<<" track len pion " << TU->TrkLen(Par) << " n hit  " << TU2.nmHit()<<" d: "<<d<<" Eloss "<<Eloss<<endl;

  //     //cout << "Pion extracted. pval = "<<pval << endl;
  //     h_pion->Fill(pval, Eloss);
  //     h_piondedx->Fill(pval, dedx);

  //     //cout << "Pion filled. Ncl = "<<Ncl << endl;
  //   }
  // }

  // graph K/pi separation
  //
  const Int_t Nint = 501; //501; // 501;
  TH1F **hEtrK = new TH1F *[Nint];
  // TCanvas **cvEtrK = new TCanvas *[Nint];
  // TCanvas **cvEtrPi = new TCanvas *[Nint];
  // TCanvas **cvElossK = new TCanvas *[Nint];
  TH1F **hEtrP = new TH1F *[Nint];
  TH1F **hElossK = new TH1F *[Nint];
  Double_t pmom[Nint];
  // vector<double>pmom;
  Double_t pmn = 0.2;
  Double_t pmx = 100.;
  Double_t stp = (TMath::Log(pmx) - TMath::Log(pmn)) / (Double_t)(Nint - 1);
  for(Int_t i = 0; i < Nint; i++)
  {
    pmom[i] = TMath::Exp(i * stp + TMath::Log(pmn));
    // cout << " pmom "<<pmom[i] << endl;
  }
  Double_t SigDiff[Nint];
  Double_t dedxK, dedxP;
  Double_t ElossK, ElossP;

  vector<double> meanEtrk[Nint], sgmEtrk[Nint], meanEtrPi[Nint], sgmEtrPi[Nint];
  vector<double> ElossKp[Nint], ElossKpErr[Nint], ElossPi[Nint], ElossPiErr[Nint];
   vector<double> meanhistEtrp[Nint],stdhistEtrp[Nint],meanhistEtrk[Nint],stdhistEtrk[Nint];
  for(Int_t i = 0; i < Nint; i++)
  {
    Double_t bgK = pmom[i] / KpMass;
    Double_t bgP = pmom[i] / PiMass;
    // cout << " bgK " << bgK << " bgP " << bgP << endl;
    Double_t px = pmom[i] * TMath::Sin(angRad);
    Double_t pz = pmom[i] * TMath::Cos(angRad);
    // cvEtrK[i] = new TCanvas(Form("cvEtrK%d", i), Form("cvEtrK for kaon %d", i));
    hEtrK[i] = new TH1F(Form("hEtrK%d", i), Form("Energy loss trunc for kaon %d", i), Nev, 0, 1);

    // cvEtrPi[i] = new TCanvas(Form("cvEtrPi%d", i), Form("cvEtrPi for pion%d", i));
    hEtrP[i] = new TH1F(Form("hEtrP%f", bgP), Form("Energy loss trunc for pion %d", i), Nev, 0, 1);

    // cvElossK[i] = new TCanvas(Form("cvElossK%d", i), Form("cvElossK for kaon%d", i));
    hElossK[i] = new TH1F(Form("hElossK%f", bgP), Form("Energy loss for kaon %d", i), Nev, 0,1);

    for(int j = 0; j < Nev; ++j)
    {

      TVector3 p(px, 0.0, pz);
      Double_t Q = 1.0;
      TVectorD Par = TrkUtil::XPtoPar(x, p, Q, Bz);
      Double_t tLen = TU->TrkLen(Par);
      dedxK = 0.0;
      ElossK = 0.0;
      if(TU->IonEnergyLoss(ElossK, KpMass, Par, dedxK))
      {
        // cout<<"dedxk "<<dedxK<<" ElossK "<<ElossK<<endl;
        // ElossKp[i].push_back(dedxK);
        hEtrK[i]->Fill(dedxK);
        // ElossKpErr[i].push_back(TMath::Sqrt(dedxK));
        hElossK[i]->Fill(ElossK);
      }
      ElossP = 0.0;
      dedxP = 0.0;
      if(TU->IonEnergyLoss(ElossP, PiMass, Par, dedxP))
      {
        // ElossPi[i].push_back(dedxP);
        // ElossPiErr[i].push_back(TMath::Sqrt(dedxP));
        hEtrP[i]->Fill(dedxP);
        // cout<<"dedxP "<<dedxP<<" ElossP "<<ElossP<<endl;
      }

      //  if(ElossKp[i] + ElossPi[i] > 0.0) SigDiff[i] = TMath::Abs(ElossPi[i] - ElossKp[i]) / (0.5 * (ElossPiErr[i] + ElossKpErr[i]));
      // cout<<"differences "<<SigDiff[i]<<endl;
    }

    // cvEtrPi[i]->cd();
    hEtrP[i]->Fit(fgausP);
    meanEtrPi[i].push_back(fgausP->GetParameter(1));
    sgmEtrPi[i].push_back(fgausP->GetParameter(2));
    meanhistEtrp[i].push_back(hEtrP[i]->GetMean());
    stdhistEtrp[i].push_back(hEtrP[i]->GetStdDev());

    // cvEtrK[i]->cd();
    hEtrK[i]->Fit(fgausK);
    meanEtrk[i].push_back(fgausK->GetParameter(1));
    sgmEtrk[i].push_back(fgausK->GetParameter(2));
    meanhistEtrk[i].push_back(hEtrK[i]->GetMean());
    stdhistEtrk[i].push_back(hEtrK[i]->GetStdDev());
    // cvElossK[i]->cd();
    // hElossK[i]->Draw();
  }

  for(Int_t i = 0; i < Nint; i++)
  {
    for(int j = 0; j < meanEtrk[i].size(); ++j)
    {

      if(meanEtrk[i].at(j) > 0.0 && meanEtrPi[i].at(j) > 0.0 && sgmEtrPi[i].at(j) > 0.0 && sgmEtrk[i].at(j) > 0.0)
      {
        // SigDiff[i] = TMath::Abs(meanEtrk[i].at(j) - meanEtrPi[i].at(j)) / (0.5 * (sgmEtrPi[i].at(j) + sgmEtrk[i].at(j))); ///(TMath::Power(112,-0.43)*TMath::Power(1.45,-0.32));
        SigDiff[i] = TMath::Abs(meanhistEtrk[i].at(j) - meanhistEtrp[i].at(j)) / (0.5 * (stdhistEtrp[i].at(j) + stdhistEtrk[i].at(j))); ///(TMath::Power(112,-0.43)*TMath::Power(1.45,-0.32));
        // cout << "sig diff " << SigDiff[i] << " at " << pmom[i] << " " << meanhistEtrk[i].at(j) << " " << meanhistEtrp[i].at(j) << " " << stdhistEtrp[i].at(j) << " " << stdhistEtrk[i].at(j) << endl;
      }
    }
  }

  //
  //plot
  //

  TCanvas *cPidEtr = new TCanvas("cPidEtr", "cPidEtr", 1000, 1000);
  cPidEtr->cd();
  cPidEtr->SetLogx();
  TGraph *gr = new TGraph(Nint, pmom, SigDiff);
  gr->SetMinimum(0.0);
  gr->SetMaximum(15.0);
  gr->GetXaxis()->SetLimits(0., 100.);
  gr->SetTitle("K/#pi separation in nr. of #sigma");
  gr->GetXaxis()->SetTitle("Momentum (GeV)");
  gr->SetLineColor(kBlue);
  // gr->Draw("APC");

  // TCanvas *cvnly = new TCanvas("cvnly", "cvnly");
  // hNly->Draw();
  // TCanvas *cvtlen = new TCanvas("cvtlen", "cvtlen");
  // cvtlen->cd();
  // htlen->Draw();
  TCanvas *cnv = new TCanvas("cnv", title, 50, 50, 800, 500);
  cnv->cd();
  gStyle->SetOptStat(0);
  h_kaon->SetMarkerColor(kRed);
  h_kaon->SetLineColor(kRed);
  h_kaon->SetMarkerStyle(22);
  h_kaon->GetXaxis()->SetTitle("Momentum (GeV)");
  h_kaon->Draw();
  h_pion->SetMarkerColor(kBlack);
  h_pion->SetLineColor(kBlack);
  h_pion->SetMarkerStyle(29);
  h_pion->Draw("SAME");
  TLegend *lg = new TLegend(0.1, 0.9, 0.3, 0.70);
  TString LgTitle = "Particle type:";
  lg->SetHeader(LgTitle);
  lg->AddEntry(h_pion, "#pi", "L");
  lg->AddEntry(h_kaon, "#color[2]{K}", "L");
  lg->Draw();
  TCanvas *cvtlen_nhit = new TCanvas("cvtlen_nhit", "cvtlen_nhit");
  htlen_nhit->Draw("colz");
  TCanvas *cvtlen_p = new TCanvas("cvtlen_p", "cvtlen_p");
  htlen_p->Draw("colz");
  TFile *MyFile = new TFile(Form("test.root"), "recreate");
  MyFile->cd();
  gr->Write();
  cnv->Write();
  for(int i = 0; i < Nint; ++i)
  {
    hEtrP[i]->Write();
    hEtrK[i]->Write();
    hElossK[i]->Write();
  }
  htlen->Write();
  hNly->Write();
  hcelllen->Write();
  cvtlen_nhit->Write();
  cvtlen_p->Write();
  MyFile->Close();
}
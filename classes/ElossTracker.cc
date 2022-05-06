/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2020  Universite catholique de Louvain (UCLouvain), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \class ElossTracker
 *
 *  Energy loss in drift chambers
 *
 *  \authors F. Bedeschi - INFN Pisa
*            P. Demin - UCLouvain, Louvain-la-Neuve
 *           M. Selvaggi - CERN
 *           F. Cuna - INFN
 *           G.F.Tassielli - INFN
 *
 *
 */

#include "modules/ElossTracker.h"
#include "TrackCovariance/SolGeom.h"
#include "TrackCovariance/SolTrack.h"
#include "classes/DelphesClasses.h"

#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TVectorD.h"




using namespace std;

//------------------------------------------------------------------------------

ElossTracker::ElossTracker() :
  fTrackUtil(0)
{
  fTrackUtil = new TrkUtil();
}

//------------------------------------------------------------------------------

ElossTracker::~ElossTracker()
{
  if(fTrackUtil) delete fTrackUtil;
}

//------------------------------------------------------------------------------

void ElossTracker::Init()
{

  // geometric acceptance
  fRmin = GetDouble("Rmin", 0.);
  fRmax = GetDouble("Rmax", 0.);
  fZmin = GetDouble("Zmin", 0.);
  fZmax = GetDouble("Zmax", 0.);

  // magnetic field
  fBz = GetDouble("Bz", 0.);

  fGasOption = GetInt("GasOption", 0);

  // initialize drift chamber geometry and gas mix
  fTrackUtil->SetBfield(fBz);
  fTrackUtil->SetDchBoundaries(fRmin, fRmax, fZmin, fZmax);
  fTrackUtil->SetGasMix(fGasOption);

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "TrackMerger/tracks"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void ElossTracker::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void ElossTracker::Process()
{
  Candidate *candidate, *mother, *particle;
  Double_t mass, trackLength, Eloss;
  Double_t dedx;
  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {

    // converting to meters
    particle = static_cast<Candidate *>(candidate->GetCandidates()->At(0));

    // converting to meters
    const TLorentzVector &candidatePosition = particle->Position * 1e-03;
    const TLorentzVector &candidateMomentum = particle->Momentum;

    TVectorD Par = TrkUtil::XPtoPar(candidatePosition.Vect(), candidateMomentum.Vect(), candidate->Charge, fBz);
    mass = candidateMomentum.M();

    trackLength = fTrackUtil->TrkLen(Par);
   
    mother = candidate;
    candidate = static_cast<Candidate *>(candidate->Clone());

    Eloss = 0.;
    dedx=0.;
     
    if(fTrackUtil->IonEnergyLoss(Eloss, mass, Par,dedx))
    {
      candidate->EnergyLoss = Eloss;
    
      candidate->dEdx = dedx;
    }

    candidate->AddCandidate(mother);

    fOutputArray->Add(candidate);
  }
}




//------------------------------------------------------------------------------
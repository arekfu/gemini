#include "CNucleus.h"
#include "CAngle.h"
#include <iostream>
#include <sstream>
#include "TFile.h"
#include "TTree.h"

// this is an example of using GEMINI CNucleus class to give the
//statistical decay of a compound nucleus

int main(int argc, char *argv[])
{
  if(argc!=6) {
    std::cerr << "Syntax: testDecay Z A EStar J nEvents" << std::endl;
    return 1;
  }
  std::stringstream ss;

  int iZCN; // proton number of compound nucleus
  ss.str(argv[1]);
  ss >> iZCN;
  ss.clear();

  int iACN; // mass number of compound nucleus
  ss.str(argv[2]);
  ss >> iACN;
  ss.clear();

  CNucleus CN(iZCN,iACN); //constructor

  float fEx; //excitation energy of compound nucleus
  ss.str(argv[3]);
  ss >> fEx;
  ss.clear();

  float fJ; // spin of compound nucleus
  ss.str(argv[4]);
  ss >> fJ;
  ss.clear();

  CN.setCompoundNucleus(fEx,fJ); //specify the excitation energy and spin

  CN.setVelocityCartesian(); // set initial CN velocity to zero
  CAngle spin(CNucleus::pi/2,(float)0.); 
  CN.setSpinAxis(spin); //set the direction of the CN spin vector

  int nEvents; // number of events
  ss.str(argv[5]);
  ss >> nEvents;
  ss.clear();

  const int maxParticles = 300;
  Short_t nParticles, A[maxParticles], Z[maxParticles];
  Float_t EKin[maxParticles], px[maxParticles], py[maxParticles], pz[maxParticles], theta[maxParticles], phi[maxParticles];
  TFile *file = new TFile("testDecay.root","recreate");
  TTree *t = new TTree("t","t");
  t->Branch("nParticles", &nParticles, "nParticles/S");
  t->Branch("A", A, "A[nParticles]/S");
  t->Branch("Z", Z, "Z[nParticles]/S");
  t->Branch("EKin", EKin, "EKin[nParticles]/F");
  t->Branch("px", px, "px[nParticles]/F");
  t->Branch("py", py, "py[nParticles]/F");
  t->Branch("pz", pz, "pz[nParticles]/F");
  t->Branch("theta", theta, "theta[nParticles]/F");
  t->Branch("phi", phi, "phi[nParticles]/F");

  for (int i=0;i<nEvents;i++) {
    if(i%1000==0)
      cout << "event = " << i << endl;
    //CN.setWeightIMF();// turn on enhanced IMF emission
    CN.decay(); //decay the compound nucleus

    if (CN.abortEvent)
    {
      CN.reset();
      continue;
    }

    nParticles = CN.getNumberOfProducts();

    for(int j=0; j<nParticles; ++j) {
      CNucleus * product = CN.getProducts(j);

      A[j] = product->iA;
      Z[j] = product->iZ;
      EKin[j] = product->getKE();
      CAngle angle = product->getAngleDegrees();
      theta[j] = angle.theta;
      phi[j] = angle.phi;
      px[j] = product->getMomentumVector()[0];
      py[j] = product->getMomentumVector()[1];
      pz[j] = product->getMomentumVector()[2];
    }

    t->Fill();
    CN.reset();
  }
  t->Write();
  file->Close();

}

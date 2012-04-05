#include "CRunThick.h"

/**
 * this constructor is everything at the moment.
/param iZcn - CN proton number
/param iAcn is the mass number
/param fEx is the excitation energy in MeV
/param l0 is the critical spin , specifying the max spin in fusion (unit hbar)
/param d0 is the diffuseness of the spin distribution (units hbar)
/param lmax is the maximum spin value considered. (<l0 for just er info)
/param plb is pi-lambda-bar squared in mb
/param numTot is number of Monte Carlo simulations
/param title0 is name of output root file (without ".root" extension)

 */
CRunThick::CRunThick(int iZcn, int iAcn, float fEx_min,float fEx_max,float l0_min, 
	   float l0_max, float d0, int lmax, float plb,int nBins,
	   int numTot,string title0,float vcm/*=0.*/, float thetaDetMin/*=0.*/,
           float thetaDetMax/*=360*/)
{

  cout << title0 << endl;
  string title = title0 + ".root";
  bool residue;
  bool residueDet;


  float ExArray[nBins];
  for (int i=0;i<nBins;i++) ExArray[i] = fEx_min + (fEx_max-fEx_min)/
    ((float) nBins)*((float)i+0.5);




  float prob[nBins][lmax+1];
  for (int i=0;i<nBins;i++)
    {
      float l0 = l0_min + (l0_max-l0_min)/((float)nBins)*((float)i+0.5);

     float sum = 0.;
     for (int l=0;l<=lmax;l++)
       {
         prob[i][l] =  (float)(2*l+1);
         if (d0 > 0.) prob[i][l] /= (1.+exp(((float)l-l0)/d0));
         else if ( l > l0) prob[i][l] = 0.;
         sum += prob[i][l];
       }
     for (int l=0;l<=lmax;l++)
       {
         prob[i][l] /= sum;
         if (l > 0) prob[i][l] += prob[i][l-1];
       }
    }


  float sum = 0.;
     for (int l=0;l<=lmax;l++)
       {
         float fact =  (float)(2*l+1);
         if (d0 > 0.) fact /= (1.+exp(((float)l-(l0_min+l0_max)/2.)/d0));
         else if ( l > (l0_min+l0_max)/2.) fact = 0.;
         sum += fact;
       }


  TFile *f = new TFile(title.c_str(),"RECREATE");
  CNucleus CN(iZcn,iAcn);

  //CNucleus::setSolution(1);
  //CNucleus::setFissionScaleFactor(7.38);
  //CNucleus::setAddToFisBarrier(-1.);
  //CNucleus::setNoIMF();
  //  CNucleus::setAddToFisBarrier(4.);
  //CNucleus::setLestone();
  //CLevelDensity::setAfAn(1.036);
  //CLevelDensity::setAimfAn(1.05);
  //CNucleus::setTimeTransient(1.);
  //CTlBarDist::setBarWidth(1.);
  //CTlBarDist::setBarWidth(0.);
  //CYrast::forceSierk();


  CN.setVelocityCartesian((float)0.,(float)0.,(float)0.);
  CAngle spin((float)0.,(float)0.);
  CN.setSpinAxis(spin);

  CN.printParameters();
  
  float asy[20]={0.};
  float asyMultPre[20] = {0.};
  float asyMultPost[20] = {0.};
  float asyMultTot[20] = {0.};
  float Nres = 0.;
  float NresDet = 0.;
  float sumAres = 0.;
  float sumAresDet = 0.;
  float Nfiss = 0.;
  float NfissLost = 0.;
  float LfissLost = 0.;
  float NpreSad = 0.;
  float NpreScis = 0.;
  float Npost = 0.;
  float Nalpha = 0.;
  float Nproton = 0.;  
  float Nneutron = 0.;
  float NLi6 = 0.;
  float NLi7 = 0.;
  float NBe7 = 0.;
  float Mfis = 0;
  float M2fis = 0.;
  float M0fis = 0.;
  int numberA = 0;
  int averageA = 0;

  TH1F histEgamma("Egamma","",100,0,50);
  TH1F histER("histER","",90,0,90);
  TH1F histERxn("histERxn","",90,0,90);
  TH1F histFis("histFis","",90,0,90);
  TH1F histFus("histFus","",90,0,90);
  TH1F histA("histA","",230,0,230);
  TH2F histAA("histAA","",200,0,200,200,0,200);
  TH1F histAFis("histAFis","",230,0,230);
  TH1F histAFisPrimary("histAFisPrimary","",230,0,230);
  TH1F histAFisPrimaryVel("histAFisPrimaryVel","",230,0,230);
  TH1F histZ("histZ","",92,0,92);
  TH1F histZ_fis("histZ_fis","",92,0,92);
  TH1F histZ_nofis("histZ_nofis","",92,0,92);
  TH1F histN("histN","",132,0,132);
  TH1F angle("angle","",180,0,180);
  TH2F histZN("histZN","",151,0,151,151,0,151);
  TH1F keFF("keFF","",150,0,150);
  TH1F kePreSad("kePreSad","",100,0,30);
  TH1F kePreSS("keSS","",100,0,30);
  TH1F kePreSc("kePreSc","",100,0,30);
  TH1F kePost("kePost","",100,0,30);
  TH1F keEvap("keEvap","",100,0,30);
  TH1F velFF("velFF","",100,0,4.);
  TH1F keAlpha("keAlpha","",50,0,50);
  TH1F keProton("keProton","",50,0,50);
  TH1F keNeutron("keNeutron","",50,0,50);
  TH1F keLi6("keLi6","",60,0,60);
  TH1F keLi7("keLi7","",60,0,60);
  TH1F keBe7("keBe7","",60,0,60);
  TH1F histFis2("histFis2","",100,0,7000);
  TH2F histAL("histAL","",250,0,250,100,0,100);
  TH2F histxnEx("histxnEx","",50,0,50,50,0,20);
  TH2F histxnExA("histxnExA","",16,200,216,50,0,20);

  bool f14=1;
  bool f12=1;
  bool f34=1;
  for (int i=0;i<numTot;i++)
    {
      float weight = 1.;
      //cout <<"event= " <<  i << endl;
      if (i > numTot*.25 && f14)
	{
	  cout << " 25%" << endl;
          f14 = 0;
	} 
      if (i > numTot*.5 && f12)
	{
	  cout << " 50%" << endl;
          f12 = 0;
	} 
      if (i > numTot*.75 && f34)
	{
	  cout << " 75%" << endl;
          f34 = 0;
	} 
      int mbin = nBins+1;
      for (;;)
	{
	  mbin= (int)floor(CN.ran.Rndm()*(float)nBins);
	  if (mbin < nBins) break;
	}

      float fEx = ExArray[mbin];


      float ran = CN.ran.Rndm();
      int l = 0;
      for (;;)
	{
	  if (ran < prob[mbin][l]) break;
          l++;
	}
      //l = 2; //rjc
      //fEx = 77.83;

      /*
      if (i==60) //rjc
	{
	  cout << "here" << endl;
	}
      */

      CN.setCompoundNucleus(fEx,(float)l);      

	{
	  //if (i%100==0)cout << "l= "<< l  << " i= " << i << endl;
	  CN.setWeightIMF();
         CN.decay();

        if (CN.abortEvent)
          {
	    CN.reset();
            continue;
          }

	histEgamma.Fill(CN.getSumGammaEnergy());
         int iStable = CN.getNumberOfProducts();

         CNucleus *productER = CN.getProducts(iStable-1);
	 weight *= productER->getWeightFactor();

	 if(productER->iZ == iZcn)
	   {
	     averageA += productER->iA;
             numberA += 1;
	   }

         int iZres = productER->iZ;
         float resEx = productER->fEx;
         float resJ = productER->fJ;
         int iAres = productER->iA;
         int multTot = 0;
         int iZ, iA;

         CNucleus * product = CN.getProducts(0);
         histFus.Fill(l,weight);

         if (CN.isResidue()) 
	   {
             float * vv;
             vv = productER->getVelocityVector();
             vv[2] += vcm;
             float vvv = sqrt(pow(vv[0],2)+pow(vv[1],2)+pow(vv[2],2));
             float AngleDeg = acos(vv[2]/vvv)*180./3.14159;
             if (AngleDeg > thetaDetMin && AngleDeg < thetaDetMax) 
	       residueDet = 1;
             else residueDet = 0;
	     residue = 1;
             histER.Fill(l,weight);
	     if (iZres == iZcn)
	         {
                    histERxn.Fill(l,weight);
		    histxnEx.Fill(resJ,resEx,weight);
		    histxnExA.Fill(iAres,resEx,weight);
	          }

             Nres += weight;
             sumAres += weight*(float)productER->iA;
             if (residueDet) 
	       {
               NresDet += weight;
               sumAresDet += weight*(float)productER->iA;
	       }
	   }
	 else 
	   {
	     residue = 0;
             residueDet = 0;
	     Nfiss += weight;
             histFis.Fill(l,weight);
             histFis2.Fill(l*l,weight);
	   }


         int iAmax = 0;
         int iAnext = 0;
         float vmax = 0.;
         float vnext = 0.;
         float emax = 0.;
         float enext = 0.;
         for (int j=0;j<iStable;j++)
	   {
             iZ = product->iZ;
             iA = product->iA;
	     //if (CN.SymmetricisFission())cout << iZ << " " << iA << endl;  //rjc
             if (iA > iAmax) 
	       {
                 iAnext = iAmax;
                 vnext = vmax;
                 enext = emax;
                 iAmax = iA;
		 emax = product->getKE();
		 vmax = product->getVelocity();
	       }
	     else if (iA > iAnext)
	       {
                 iAnext = iA;
                 enext = product->getKE();
                 vnext = product->getVelocity();
	       }
             //cout << iZ << " " << iA << endl;
             if (product->getTime() < 0.) 
	       {
		 cout << "negative time" << endl;
                 cout << iZ << " " << iA << " " << 
		   product->getParent()->iZ << " " << 
		   product->getParent()->iA << " " 
		      << product->getParent()->fEx << " " 
		      << product->getParent()->fJ << endl;
	       }

             histZ.Fill(iZ,weight);
             if (CN.isSymmetricFission())histZ_fis.Fill(iZ,weight);
             else histZ_nofis.Fill(iZ,weight);
             histA.Fill(iA,weight);

             histAL.Fill(iA,l,weight);
             histN.Fill(iA-iZ,weight);
             histZN.Fill(iA-iZ,iZ,weight);
             if (iZ == 0 && iA == 1) 
	       {
		 if (residueDet) //iARes >= Ares)
		   {
   		   keNeutron.Fill(product->getKE(),weight);
                   Nneutron += weight;
		   }
                 multTot++;
                 if (CN.isSymmetricFission())
		   {
                    if (product->origin == 0)
		      { 
		      kePreSad.Fill(product->getKE(),weight);
		      NpreSad += weight;
		      }
                    if (product->origin == 1) 
                      kePreSS.Fill(product->getKE(),weight);
                    if (product->origin <= 1) 
		      {
			NpreScis += weight;
                      kePreSc.Fill(product->getKE(),weight);
		      }
		    if (product->origin > 1) 
		      {
			Npost += weight;
                      kePost.Fill(product->getKE(),weight);
		      }
		   }
		 else keEvap.Fill(product->getKE(),weight);
	       }
	     else if (iZ == 1 && iA == 1 && residueDet) //iARes >= Ares)
	       {
		keProton.Fill(product->getKE(),weight);
                Nproton += weight;
	       }
	     else if (iZ == 2 && iA == 4)
	       {
		 if(residueDet) //iARes >=Ares )
		   {

		     /*
                    //rjc
                     if (product->getKE() < 15.)
		       {
                       cout << " i = " << i << endl;
		       cout << product->getParent()->iZ << " " <<
			 product->getParent()->iA << " " <<
			 product->getParent()->fEx << " " <<
			 product->getParent()->fJ << " " <<
			 product->getKE() << " " <<
			 product->getParent()->daughterHeavy->fJ  << " " <<
                         product->getParent()->daughterHeavy->fEx << endl;
		       }
		     */
		     //cout << "alpha " << product->getKE() << endl; //rjc
	             keAlpha.Fill(product->getKE(),weight);
                     Nalpha += weight;
		   }
		 }
	     else if (iZ == 3 && iA == 6)
	       {
		 if(residueDet) //iARes >=Ares )
		   {

	             keLi6.Fill(product->getKE(),weight);
                     NLi6 += weight;
		   }
		 }
	     else if (iZ == 3 && iA == 7)
	       {
		 if(residueDet) //iARes >=Ares )
		   {

	             keLi7.Fill(product->getKE(),weight);
                     NLi7 += weight;
		   }
		 }
	     else if (iZ == 4 && iA == 7)
	       {
		 if(residueDet) //iARes >=Ares )
		   {

	             keBe7.Fill(product->getKE(),weight);
                     NBe7 += weight;
		   }
		 }
             if (iZ > 1 && iZ < 5)
	       {
                 angle.Fill(product->getThetaDegrees(),weight);
	       }
             if (iZ > 5 && CN.isSymmetricFission()) 
	       {
                 keFF.Fill(product->getKE(),weight);
                 velFF.Fill(product->getVelocity(),weight);
                 Mfis += (float)product->iA;
                 M2fis += pow((float)product->iA,2);
                 M0fis += 1.;
                 histAFis.Fill(product->iA,weight);
	       }
	     product=CN.getProducts();
	    }
         if (CN.isSymmetricFission())
	   {
             if ((float)iAmax > 0.77*(float)CN.iA)
	       {
                NfissLost += weight;
                LfissLost += weight*CN.fJ;
	       }                
             float A2 = emax/(emax+enext)*(float)CN.iA;
             float A1 = (float)CN.iA - A2;

             histAFisPrimary.Fill(A1,weight);
             histAFisPrimary.Fill(A2,weight);
             
             A2 = vmax/(vmax+vnext)*(float)CN.iA;
             A1 = (float)CN.iA - A2;

             histAFisPrimaryVel.Fill(A1,weight);
             histAFisPrimaryVel.Fill(A2,weight);
             //cout << iAmax << " " << iAnext << " " << A2 << " " << A1 << endl;
	     histAA.Fill((float)iAnext,A2,weight);
	     histAA.Fill((float)iAmax,A1,weight);
	   }

         float Amax = (float)iAmax/(float)(iAmax+iAnext)*162.;
         float Anext = (float)iAnext/(float)(iAmax+iAnext)*162.;
	 int iasy = (int)(Amax/10);
         asy[iasy] += weight;
         asyMultPre[iasy] += weight*(float)CN.getMultPre();
         asyMultPost[iasy] += weight*(float)CN.getMultPost();
         asyMultTot[iasy] += weight*(float)multTot;
	 iasy = (int)(Anext/10);
         asy[iasy] += weight;
         asyMultPre[iasy] += weight*(float)CN.getMultPre();
         asyMultPost[iasy] += weight*(float)CN.getMultPost();
         asyMultTot[iasy] += weight*(float)multTot;
   

        CN.reset();
	}
      }

  title = title0+"M.dat";
  ofstream ofFile(title.c_str());
  for (int i=0;i<20;i++)
    {
      if (asy[i] == 0) continue;
      ofFile << i*10 + 5 << " " << asyMultPre[i]/asy[i] << " " <<
	asyMultPost[i]/asy[i] << " " << asyMultTot[i]/asy[i] << endl;
    }

  histA.Scale(plb/(float)numTot*sum);
  histZ.Scale(plb/(float)numTot*sum);
  histZ_fis.Scale(plb/(float)numTot*sum);
  histZ_nofis.Scale(plb/(float)numTot*sum);
  histN.Scale(plb/(float)numTot*sum);
  histAFis.Scale(plb/(float)numTot*sum);
  histAFisPrimary.Scale(plb/(float)numTot*sum);
  histAFisPrimaryVel.Scale(plb/(float)numTot*sum);
  histZN.Scale(plb/(float)numTot*sum);

  histER.Write();
  histERxn.Write();
  histxnEx.Write();
  histxnExA.Write();
  histFis.Write();
  histFis2.Write();
  histFus.Write();
  angle.Write();
  histZ.Write();
  histZ_fis.Write();
  histZ_nofis.Write();
  histA.Write();
  histAA.Write();
  histAFisPrimary.Write();
  histAFisPrimaryVel.Write();
  histAL.Write();
  histN.Write();
  histZN.Write();
  keFF.Write();
  kePreSad.Write();
  kePreSad.Write();
  kePreSc.Write();
  kePost.Write();
  keEvap.Write();
  histEgamma.Write();
  keAlpha.Scale(1./NresDet);
  keProton.Scale(1./NresDet);
  keNeutron.Scale(1./NresDet);
  keLi6.Scale(1./NresDet);
  keLi7.Scale(1./NresDet);
  keBe7.Scale(1./NresDet);
  keAlpha.Write();
  keProton.Write();
  keNeutron.Write();
  keLi6.Write();
  keLi7.Write();
  keBe7.Write();
  velFF.Write();
  f->Write();
  cout << "NresDet= " << NresDet << " Nneut= " << Nneutron << " NProt= " <<
    Nproton << " Nalpha= " << Nalpha << " NLi6= " << NLi6 << " NLi7= " << NLi7
       << " NBe7= " << NBe7 << endl;
  cout << "Li6 mult = " << NLi6/NresDet << endl;
  cout << "Li7 mult = " << NLi7/NresDet << endl;
  cout << "Be7 mult = " << NBe7/NresDet << endl;
  cout << "neutron mult= " << Nneutron/NresDet << endl;
  cout << "proton mult= " << Nproton/NresDet << endl;
  cout << "alpha mult= " << Nalpha/NresDet << endl;

  cout << " mean ER A = " << sumAres/Nres << endl;
  cout << " for det res = " << sumAresDet/NresDet << endl;

  float xER = Nres/(float)numTot*sum*plb;
  float xFiss2 = Nfiss/(float)numTot*sum*plb;
  if (NfissLost > 0) LfissLost /= NfissLost;
  float xFissLost = NfissLost/(float)numTot*sum*plb;
  cout << "sigmaER = " << xER << " mb " << endl;

  float xFus = 0.;
  for (int l=0;l<200;l++)
    {
      float xx =  (float)(2*l+1);
      if (d0 > 0.) xx /=(1.+exp(((float)l-(l0_max+l0_min)/2.)/d0));
      else if (l > (l0_min+l0_max)/2.) break;
      xFus += xx;
    }
  xFus *= plb;
  float xFis = xFus - xER;

  cout << "fusion xsec= " << xFus << " mb" << endl;
  cout << "fission xsec= " << xFis << " mb " << xFiss2 <<  " " 
       << xFissLost << " " << LfissLost<<   endl;
  cout << "preSaddle neut mult = " << NpreSad/Nfiss << endl;
  cout << "preScis neut mult = " << NpreScis/Nfiss << endl;
  cout << "post neut mult = " << Npost/Nfiss << endl;

  float Mav = Mfis/M0fis;
  cout << "mean fission mass = " << Mav << endl;
  float sigma2 = M2fis/(M0fis-1) - M0fis/(M0fis-1)*pow(Mav,2);
  //float sigma = sqrt(sigma2);
  cout << "sigma2M= " << sigma2 << endl;

  if (numberA > 0) cout << "average x for xn products is " << (float)iAcn-(float)averageA/(float)numberA << endl;
}

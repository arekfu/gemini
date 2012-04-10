#include "CNucleus.h"
bool const Isig  = 1;    //in weisskopf, use parametrized Inverse Section
                         // otherwise calculated them from transmission coeff.
int  CNucleus::iHF = 2; //1=Hauser Feshback,0= Weisskopf,2=switches from 0 to 1
bool const fissionMassScission = 0; //1=scission-point fission mass dist
                                    //0=saddle-point fission mass dist
float const scaleImf = 1.;
bool const Isaddle = 0;
bool const Iscission = 1;
//float const WignerScaled = 0.; //scaling factor for the Wigner Energy
float const WignerAdd = 7.; //adding factor for the Wigner Energy
//the following line needed in the ROOT version
//ClassImp(CNucleus)

bool          CNucleus::noIMF = 0;
bool          CNucleus::BohrWheeler = 1;
CRandom       CNucleus::ran;
CYrast        CNucleus::yrast; 
CScission     CNucleus::scission;
CLevelDensity CNucleus::levelDensity;
CAngleDist    CNucleus::angleDist;
float         CNucleus::de = 1.0;

int const CNucleus::nStore=5000;
int const CNucleus::nSub = 800;
int const CNucleus::nSubTl = 80;

int const CNucleus::Nproducts =180;
int const CNucleus::Nstable=100;
CNucleus* CNucleus::allProducts[Nproducts]={0};
CNucleus* CNucleus::stableProducts[Nstable]={0};
int CNucleus::iProducts=0;
int CNucleus::iStable=0;
float const CNucleus::r0=1.16;
float const CNucleus::sep=2.;
float const CNucleus::pi=acos(-1.);
short unsigned CNucleus::Zshell = 2;
short unsigned const CNucleus:: lMaxQuantum = 50; 
CEvap CNucleus::evap;
float const CNucleus::gammaInhibition[3]={0.,.025,9.};
float const CNucleus::wue[3]={0.,6.8e-8,4.9e-14};//gives weisskopf units in MeV
int const CNucleus::nGamma = 4000;
float const CNucleus::viscosity_scission = 1.;
float const CNucleus::viscosity_saddle = 1.5;
float CNucleus::timeTransient = 0.;
//float CNucleus::fissionScaleFactor = 2.46;
float CNucleus::fissionScaleFactor = 1.00;
float CNucleus::sumGammaEnergy = 0.;
float CNucleus::barAdd = 0.;
float const CNucleus::kRotate = 41.563;
bool const CNucleus::noSymmetry = 1; // no symmetric fission calculated in
                                 // asyFissionWidth if there is a fission peak
int CNucleus::iPoint = -1;                           
float CNucleus::threshold = .001;
//*******************************************
/**
 * constructor specifies the isotope.
\param iZ0 is the proton number
\param iA0 is the mass number
 */
CNucleus::CNucleus(int iZ0, int iA0) : CNuclide(iZ0,iA0)
{
  bStable = 0;
  saddleToSciss = 0;
  timeSinceSaddle = 0.;
  velocity[0] = 0.;
  velocity[1] = 0.;
  velocity[2] = 0.;
  spin = CAngle((float)0.,(float)0.); 
  daughterLight = NULL;
  daughterHeavy = NULL;
  parent = NULL;
  abortEvent = 0;
}
//*******************************************************
/**
 * destructor
 */
CNucleus::~CNucleus()
{
  //destructor
   if (daughterLight != NULL)
     {
       delete daughterLight;
       daughterLight = NULL;

     }
  if (daughterHeavy != NULL)
     {
        delete daughterHeavy;
        daughterHeavy = NULL;
     }
} 

//**************************************************
/**
 * Constructor specifying more parameters
\param iZ0 is proton number
\param iA0 is mass number
\param fEx0 is excitation energy in MeV
\param fJ0 is spin in hbar
*/
CNucleus::CNucleus(int iZ0, int iA0, float fEx0, float fJ0) 
    : CNuclide(iZ0,iA0)
{
  //alternative constructor
  bStable = 0;
  saddleToSciss = 0;
  timeSinceSaddle = 0.;
  velocity[0] = 0.;
  velocity[1] = 0.;
  velocity[2] = 0.;
  spin = CAngle((float)0.,(float)0.); 
  daughterLight = NULL;
  daughterHeavy = NULL;
  parent = NULL;
  abortEvent = 0;
  setCompoundNucleus(fEx0,fJ0);
}
//**************************************************
/**
 * prints out the parameters of the nucleus
 */
void CNucleus::print()
{
  //prints out paramters of nucleus
  cout << strName << endl;
  cout << "Z= " << iZ << " A= " << iA << endl;
  cout << "mass excess = " << fExpMass << endl;
  cout << strChemName << endl;
  cout << "Ex= " << fEx << " J= " << fJ << endl;
  cout << "velocity= " << velocity[0] << " " << velocity[1] << " " 
       << velocity[2] << " cm/ns" << endl;
  cout << "KE= " << getKE() << " MeV" << endl;
  cout << "Ex = " << fEx << endl;
  cout << "spin axis, theta = " << spin.theta*180./pi << " phi= "
       << spin.phi*180/pi << " deg" << endl;
  if (origin == 0) cout << " preSaddle emission" << endl; 
  else if (origin == 1) cout << "saddle-to-scission emission" << endl;
  else if (origin == 2) cout << "from light  fission fragment" << endl;
  else if (origin == 3) cout << "from heavy  fission fragment" << endl;


  cout << "emission time = " << timeSinceStart << " zs" << endl; 

}
//******************************************************
/**
 * produces a single binary decay of the nucleus from statistical-model widths
 */

void CNucleus::binaryDecay()
{
  if (bStable) return; //decay not possible



  if (notStatistical)
    {
     daughterLight = new CNucleus(evap.decay[notStatisticalMode].Z1,
                               evap.decay[notStatisticalMode].A1);
     daughterHeavy = new CNucleus(iZ-evap.decay[notStatisticalMode].Z1,
                                iA-evap.decay[notStatisticalMode].A1);
     daughterLight->origin = origin;
     daughterHeavy->origin = origin;

     daughterLight->parent = this;
     daughterHeavy->parent = this;


     float decayTime =  
       ran.expDecayTime(evap.decay[notStatisticalMode].gamma);
     daughterLight->timeSinceStart = timeSinceStart + decayTime;
     daughterHeavy->timeSinceStart = timeSinceStart + decayTime;



     daughterLight->excite(0.,evap.decay[notStatisticalMode].S1);
     daughterHeavy->excite(0.,evap.decay[notStatisticalMode].S2);
     daughterLight->bStable = 1;
     daughterLight->iWeight = 0;
     daughterHeavy->iWeight = 0;
     daughterLight->runningWeight = runningWeight;
     daughterHeavy->runningWeight = runningWeight;
     daughterLight->fact = fact;
     daughterHeavy->fact = fact;

      
     //make sure 8Be fragments decay
     if (iZ-evap.decay[notStatisticalMode].Z1 == 4 && 
         iA-evap.decay[notStatisticalMode].A1 == 8 ) 
       {
         daughterHeavy->notStatistical = 1;
         daughterHeavy->notStatisticalMode = 20;
         daughterHeavy->bStable = 0;        
       }
     else daughterHeavy->bStable = 1;
     EvapLPlusS = evap.decay[notStatisticalMode].lPlusS1;
     EvapL = evap.decay[notStatisticalMode].L;
     EvapS1 = evap.decay[notStatisticalMode].S1;
     EvapS2 = evap.decay[notStatisticalMode].S2;
     EvapEk = evap.decay[notStatisticalMode].Ek;
     EvapA1 = evap.decay[notStatisticalMode].A1;
     EvapA2 = iA - EvapA1;   
     

     angleEvap(); // find the emission angles of the two fragments
     allProducts[iProducts] = daughterLight;
     iProducts++;
     if (iProducts == Nproducts)
       {
         cout << "increase Nproducts" << endl;
	 abort();
       }

     allProducts[iProducts] = daughterHeavy;
     iProducts++;
     if (iProducts == Nproducts)
       {
         cout << "increase Nproducts" << endl;
	 abort();
       }

     return;          
    }

  //"statistical" decays into alpha plus ligher partner
  if (fEx > 0. && (iZ == 3 || (iZ == 2&& iA == 5)))
    {
      daughterLight = new CNucleus(iZ-2,iA-4);
      daughterHeavy = new CNucleus(2,4);

      double Ek = fExpMass + fEx - daughterLight->fExpMass - 
	daughterHeavy->fExpMass;

      if (Ek <= 0.)
	{
	  delete daughterLight;
          delete daughterHeavy;
          daughterLight = NULL;

          daughterHeavy = new CNucleus(iZ,iA);
          daughterHeavy->origin = origin;
          daughterHeavy->parent = this;
          daughterHeavy->timeSinceStart = timeSinceStart;
          daughterHeavy->excite(0.,0.);
          daughterHeavy->iWeight = iWeight;
          daughterHeavy->runningWeight = runningWeight;
          daughterHeavy->fact = fact;
          daughterHeavy->bStable = 1;
	  sumGammaEnergy += fEx;
          allProducts[iProducts] = daughterHeavy;
          iProducts++;
          if (iProducts == Nproducts)
             {
               cout << "increase Nproducts" << endl;
	       abort();
              }

          return;
	}
      else 
	{

         daughterLight->origin = origin;
         daughterHeavy->origin = origin;

         daughterLight->parent = this;
         daughterHeavy->parent = this;



         daughterLight->timeSinceStart = timeSinceStart;
         daughterHeavy->timeSinceStart = timeSinceStart;



         daughterLight->excite(0.,0.);
         daughterHeavy->excite(0.,0.);
         daughterLight->bStable = 1;
         daughterHeavy->bStable = 1;
         daughterLight->iWeight = 0;
         daughterHeavy->iWeight = 0;
         daughterLight->runningWeight = runningWeight;
         daughterHeavy->runningWeight = runningWeight;
         daughterLight->fact = fact;
         daughterHeavy->fact = fact;


         //randomise angle
         float theta = acos(1.-2.*ran.Rndm());
         float phi = 2.*pi*ran.Rndm();
	 float vrel = sqrt(2.*Ek*(float)iA/(((float)iA-4.)*4.))*0.9794;
         float v1 = vrel*4./(float)iA;
         daughterLight->velocity[0] = v1*sin(theta)*cos(phi);
         daughterLight->velocity[1] = v1*sin(theta)*sin(phi);
         daughterLight->velocity[2] = v1*cos(theta);

	 float v2 = vrel - v1;
	 theta = pi - theta;
         phi += pi;

         daughterHeavy->velocity[0] = v2*sin(theta)*cos(phi);
         daughterHeavy->velocity[1] = v2*sin(theta)*sin(phi);
         daughterHeavy->velocity[2] = v2*cos(theta);

         for (int i=0;i<3;i++)
          {
           daughterLight->velocity[i] += velocity[i];
           daughterHeavy->velocity[i] += velocity[i];
          }

          allProducts[iProducts] = daughterLight;
          iProducts++;
          if (iProducts == Nproducts)
            {
              cout << "increase Nproducts" << endl;
	      abort();
            }

          allProducts[iProducts] = daughterHeavy;
          iProducts++;
          if (iProducts == Nproducts)
            {
              cout << "increase Nproducts" << endl;
	      abort();
             }
	 return;
	}

    }



  // start change

  float widthAsyFission = 0.;
  needSymmetricFission = false;
  // don;t call asymmetric fission if all asymmetries are handled 
  // by the evaporation folmalism
  if (iZ/2 > evap.maxZ)
    {
     if (noIMF)   needSymmetricFission = true;
     else widthAsyFission = asyFissionWidthZA();
    }
  float widthSymFission = 0.;
  if (pow((float)iZ,2)/(float)iA < 20.) needSymmetricFission = false;

  //end change

  float widthEvaporation = evaporationWidth();
  
  float widthGamma = gammaWidth();
  float width = widthAsyFission + widthEvaporation + widthGamma;
  float decayTime = timeSinceStart;
  

  
  width = widthEvaporation + widthGamma; 
  decayTime += ran.expDecayTime(width);
  if (decayTime < timeTransient)
    {
      widthAsyFission = 0.;
    }
  else     {
     width += widthAsyFission;
     decayTime  = timeTransient + ran.expDecayTime(width);
    }
  

  //transient delay for symmetry fission only
  if(timeSinceStart < timeTransient && needSymmetricFission)
    {
      decayTime += ran.expDecayTime(width);
      if (decayTime < timeTransient) 
	{
         needSymmetricFission = 0;
	}
      else decayTime = timeTransient;
    }

  if (needSymmetricFission) 
    {
      if (BohrWheeler)widthSymFission = BohrWheelerWidth();
      else widthSymFission = LestoneFissionWidth();
     width += widthSymFission;
     decayTime += ran.expDecayTime(width);
    }
  //end transient delay


  if (isnan(width))
    {
      cout << " a non valid decay with width (nan) was obtained for " <<
	" Z= " << iZ << " A= " << iA << " Ex= " << fEx << " j= " << fJ << endl;
      cout << "this event will be aborted " << endl;
      cout << "please report this to Bob Charity " << endl;
      abortEvent = 1;
      return; 
   }

  if (width == 0)
    {
      bStable = 1;
      return;
    }
  if (widthGamma/width > 0.99) //terminate particle evaporation.
    {
      bStable = 1;
      sumGammaEnergy += fEx;
      return;
    }


  float xran = ran.Rndm();

  if (fEx/(float)iA > 6.) widthSymFission = 0.; 
                                        //  I switch off fission at excitation
                                        //  energies above 6 MeV.  Above 
                                        //  this we get conceptable problems
                                        //  with most, or all, of the mass lost
                                        //  during the saddle-to-scission 
                                        //  transition
 
  int iChan = chooseChannel(widthEvaporation,widthAsyFission,widthSymFission,widthGamma,xran);


  if (iChan == 0) //light particle evaporation
    {

     // choose light [particle channelchannel
     float xran = ran.Rndm();
     int i = 0;
     for (;;)
       {
         float prob = evap.prob[i]/widthEvaporation;
         if (prob > xran) break;
         if ( i == evap.nLight-1) break;
         i++;
       }

      lightP = evap.lightP[i];

     EvapZ1 = lightP->iZ;
     EvapA1 = lightP->iA;
     EvapZ2 = iZ - EvapZ1;
     EvapA2 = iA - EvapA1;
     EvapMode = i;
     EvapEx1 = lightP->fEx;
     EvapS1 = lightP->fJ;

     getSpin((bool)0);

     // find the sum of L + S1
     float lPlusSMin1 = fabs((float)EvapL-EvapS1);
     float lPlusSMin2 = fabs(fJ-EvapS2);
     float lPlusSMin = max(lPlusSMin1,lPlusSMin2);
     float lPlusSMax1 = (float)EvapL+EvapS1;
     float lPlusSMax2 = fJ + EvapS2;
     float lPlusSMax = min(lPlusSMax1,lPlusSMax2);

     //the following two lines makes sure random!=1.000
      float random = 1.5;
      while(floor(random) > 0.5)random = ran.Rndm();
      EvapLPlusS = lPlusSMin + floor(random*(lPlusSMax - lPlusSMin + 1));


     daughterLight = new CNucleus(EvapZ1,EvapA1);
     daughterHeavy = new CNucleus(EvapZ2,EvapA2);
     

     daughterLight->origin = origin;
     daughterHeavy->origin = origin;

     daughterLight->iWeight = 0;
     daughterHeavy->iWeight = iWeight;
     daughterLight->runningWeight = runningWeight;
     daughterHeavy->runningWeight = runningWeight;
     daughterLight->fact = fact;
     daughterHeavy->fact = fact;


     daughterLight->parent = this;
     daughterHeavy->parent = this;

     daughterLight->timeSinceStart = timeSinceStart + decayTime;
     daughterHeavy->timeSinceStart = daughterLight->timeSinceStart;

     daughterLight->excite(EvapEx1,EvapS1);
     daughterHeavy->excite(EvapEx2,EvapS2);


     if (evap.decay[EvapMode].Ek > 0)
       {
         daughterLight->notStatistical = 1;
         daughterLight->notStatisticalMode = EvapMode;
         daughterLight->bStable = 0;
       }
     else daughterLight->bStable = 1;

     angleEvap(); // find the emission angles of the two fragments
     allProducts[iProducts] = daughterLight;
     iProducts++;
     if (iProducts == Nproducts)
       {
         cout << "increase Nproducts" << endl;
	 abort();
       }

     allProducts[iProducts] = daughterHeavy;
     iProducts++;
     if (iProducts == Nproducts)
       {
         cout << "increase Nproducts" << endl;
	 abort();
       }


    }
  else if (iChan == 1) //complex fragment or asymmetric fission decay
    {

      getCompoundNucleus()->bResidue = false;
      getCompoundNucleus()->bAsymmetricFission = true;
     daughterLight = new CNucleus(fissionZ,fissionA);
     daughterHeavy = new CNucleus(iZ-fissionZ,iA-fissionA);

     daughterLight->origin = 2;//origin;
     daughterHeavy->origin = 3;//origin;

     daughterLight->iWeight = 0;
     daughterHeavy->iWeight = 0;
     daughterLight->runningWeight = runningWeight;
     daughterHeavy->runningWeight = runningWeight;
     daughterLight->fact = fact;
     daughterHeavy->fact = fact;

     daughterLight->parent = this;
     daughterHeavy->parent = this;

     daughterLight->timeSinceStart = timeSinceStart + decayTime;
     daughterHeavy->timeSinceStart = daughterLight->timeSinceStart;

     // Coulomb part of kinetic energy
     scission.init(iZ,iA,fJ,iChan);
     scission.getFissionKineticEnergy(fissionZ,fissionA);
     float Ek = scission.ekTot;
     levelDensity.getLogLevelDensitySpherical(iA,fissionU,(float)0.,(float)0.,
				       fJ,fMInertia,2);
     //add in fluctuations to Coulomb barrier
     float sigma = Ek*sqrt(levelDensity.getTemp())*.1;
     //float sigma = 0.;
     // fluctuations are in Coulomb energy are gaussian, but make sure the
     //Coulomb energy is always positive
     // and thermal excitation energy is positive



     float Qvalue = fExpMass 
       - daughterLight->fExpMass - daughterHeavy->fExpMass
       -  scission.Erotate1 - scission.Erotate2;  //Qvalue apart from EK

     for (;;)
       {
       Ek = ran.Gaus(scission.ekTot,sigma);
       fissionU  = fEx + Qvalue - Ek;
       if (Ek > 0. && fissionU > 0.) break;
       }
     //Qvalue -= Ek;
    
     // add the thermal part of kinetic energy
     int const nE=200;
     float sumE[nE];
     int ie = 0;
     float probStart = levelDensity.
       getLogLevelDensitySpherical(iA,fissionU,(float)0.,(float)0.,(float)0.,
           fMInertia);
     for (;;)
       {
	 float de = (float)ie*0.5 + 0.25;
	 if (de >= fissionU) break;
         float prob = exp(levelDensity.
	 getLogLevelDensitySpherical(iA,fissionU-de,0.,0.,0.,fMInertia)
         - probStart);
         sumE[ie] = prob;
         if (ie > 0) sumE[ie] += sumE[ie-1];
         ie ++;
         if (prob < .01) break;
         if (ie >= nE) break;

       }
     float de = 0.;
     if (ie > 0)
       {
        xran = ran.Rndm(); 
        int iie = 0;
        for (;;)
          {
	    if (xran < sumE[iie]/sumE[ie-1]) break;
            iie++;
          }
         de = (float)iie*0.5 + 0.25;
         Ek += de;
         fissionU -= de;
       }
     scission.ekTot = Ek;


     asyFissionDivide(); // find spin and excitation energy of two fragments


     allProducts[iProducts] = daughterLight;
     iProducts++;
     if (iProducts == Nproducts)
       {
         cout << "increase Nproducts" << endl;
	 abort();
       }

     allProducts[iProducts] = daughterHeavy;
     iProducts++;
     if (iProducts == Nproducts)
       {
         cout << "increase Nproducts" << endl;
	 abort();
       }


    }
  else if (iChan == 2) //fission decay
    {
      getCompoundNucleus()->bResidue = false;
      getCompoundNucleus()->bSymmetricFission = true;

      // start saddle to scission transition
      float Esaddle = yrast.getSymmetricSaddleEnergy(iZ,iA,fJ);
      Esaddle -= fPairing + fShell;
      scission.init(iZ,iA,fJ,2);
      float Escission = scission.getScissionEnergy() - fExpMass;


      daughterHeavy = new CNucleus(iZ,iA);

      //if mass division is determined from energy at saddle point
      if (fissionMassScission == 0) 
	{
	  //asyFissionWidthBW(); //try this instead of massAsymmetry() 
	                         //if you want 
	                         //mass distribution from Sierk's 
	                         //asymmetric barriers

	 massAsymmetry(Isaddle);
         bool const notSym = 0;
         daughterHeavy->fissionZ = fissionZ;
         daughterHeavy->fissionA = fissionA;
         daughterHeavy->fissioningA = iA;
         daughterHeavy->fissioningZ = iZ; 
         daughterHeavy->sigma2 = sigma2;


         daughterHeavy->iWeight = iWeight;
         daughterHeavy->runningWeight = runningWeight;
         daughterHeavy->fact = fact;



         daughterHeavy->exciteScission(fEx,fJ,notSym);

         if (Escission >= Esaddle) timeScission = 0.;
         else timeScission = (Esaddle-Escission)*viscosity_saddle;

	}
      //else the mass division is determined at scission
      else 
	{
         daughterHeavy->exciteScission(fEx,fJ);
         if (Escission >= Esaddle) timeScission = 0.;
         else timeScission = (Esaddle-Escission)*viscosity_scission;
	}

      daughterHeavy->timeSinceSaddle = 0.;
      daughterHeavy->timeSinceStart = timeSinceStart + decayTime;
      daughterHeavy->origin = 1;
      daughterHeavy->parent = this;
      daughterHeavy->setVelocityCartesian(velocity[0],velocity[1],velocity[2]);
      daughterHeavy->setSpinAxis(spin);
      daughterHeavy->timeScission = timeScission;
      daughterLight = NULL;      


     allProducts[iProducts] = daughterHeavy;
     iProducts++;
     if (iProducts == Nproducts)
       {
         cout << "increase Nproducts" << endl;
	 abort();
       }


    }
  else //gamma decay
   {

     sumGammaEnergy += fEx - GammaEx;
     daughterLight = NULL;
     daughterHeavy = new CNucleus(iZ,iA);
     daughterHeavy->origin = origin;
     daughterHeavy->parent = this;
     daughterHeavy->timeSinceStart = timeSinceStart + decayTime;
     daughterHeavy->excite(GammaEx,GammaJ);
     daughterHeavy->iWeight = iWeight;
     daughterHeavy->runningWeight = runningWeight;
     daughterHeavy->fact = fact;
     
     angleGamma();
     allProducts[iProducts] = daughterHeavy;
     iProducts++;
     if (iProducts == Nproducts)
       {
         cout << "increase Nproducts" << endl;
	 abort();
       }


   }

}

//************************************************************
/**
 * fission mass division - uses the Rusanov systemtics to
 * determine the fission-fragment mass distribution.
 * \param saddleOrScission is true then the third Rusanov systematics
 * are used, i.e., variance is determined from the temperature at the 
 * scission point
 * if false, the second Rusanov systematics is used, i.e, the 
 * variance is determined from the temperature at the saddle point
 */
void CNucleus::massAsymmetry(bool saddleOrScission)
{
  // variance from systemacs
 
  if (saddleOrScission == 1)//asymmetry determined at scission point
     sigma2 = scission.sigmaFissionSystematicsScission(iZ,iA,fJ,fU0);
  else //asymmetry determined at saddle point
    {
     sigma2
      = scission.sigmaFissionSystematicsSaddle(iZ,iA,fJ,fEx-symSaddlePoint);
     fU0 = fEx - symSaddlePoint;
    }
  
  //at this point fU0 is the thermal excitation energy above the 
  //configuration for symmetric division. 

 
  int const nStore = 1000;
  SStore store[nStore];
  int iStore = 0;
  int iZ1,iZ2,iA1,iA2;
  iA1 = iA/2;
  float Amax = 0.;
  //float gammaAOld = 1e32;
  double gammaTot = 0.;
  for (;;)
    {

      iA2 = iA - iA1;

      int iZ1Start= (int)((float)iA1/(float)iA*(float)iZ);
      iZ1 = max(3,iZ1Start-8);
      float Zmax = 0.;
      float gammaA = 0.;
      for (;;)
        {
          iZ2 = iZ - iZ1;
          if (iZ2 <= (float)iZ/17.) break;
          if (iA1 == iA2 && iZ1 > iZ2) break;  //no double counting
          //see if fragments are listed in mass table
          if (mass->chart.getIndex(iZ1,iA1) < 0 ||   
              mass->chart.getIndex(iZ2,iA2) < 0)     
	    {
	      iZ1++;
              continue;
	    }

          float scissionE = scission.getScissionEnergy(iZ1,iA1);
          // fU0 has the asymmetric saddle or scission energy removed

          float U = fU0 - scissionE + scission.Esymmetric;
          if (U <= 0.) 
	    {
              if (iZ1 > iZ1Start) break;
	      iZ1++;
              continue;
	    }

          float logLevelDensitySciss = 
	    levelDensity.getLogLevelDensityScission(iA,U,10.);
 

          if (logLevelDensitySciss == 0) 
       	    {
              if (iZ > iZ1Start) break;
	      iZ1++;
              continue;
	    }
          float gamma= exp(logLevelDensitySciss-logLevelDensity);
          if (iA1 == iA2 && iZ1==iZ2 ) gamma /= 2.;
          gammaA += gamma;

          store[iStore].gamma = (double)gamma;
          store[iStore].iZ = iZ1;
          store[iStore].iA = iA1;
          iStore++;
          if (iStore >= nStore) 
	    {
	      cout << " increase nStore in saddleToScission " << endl;
              abort();
	    }
          Zmax = max(gamma,Zmax);
          if (gamma < Zmax*.01) break;
          iZ1++;
        }

      // rescale widths so the sum for each A has a gausiian distribution
      //with the widths from systematics 
      if (gammaA > 0.)
         {
	   float prob = exp(-pow((float)(iA1-iA2)/2.,2)/2./sigma2);
           if (iA1 == iA2) prob /= 2.;

           double scale = prob/gammaA;
           int jj = iStore-1;
           for (;;)
	     {
               if (jj < 0) break;
	       if (store[jj].iA != iA1) break;
               store[jj].gamma *= scale;
               jj--;
             }

          gammaTot += prob;
          Amax = max(Amax,prob);
          if (prob < Amax*0.001) break;
        }
      else break;
      iA1--;
      if (iA1 < (float)iA/17.) break;
    }
  if (iStore == 0) 
    {
      //this shouldn't be possible, but this statement is here just in 
      //case or else the program would get stuck in the following loop
      fissionZ = iZ/2;
      fissionA = iA/2;
      cout << "iStore == 0 in saddle to scission" << endl;
      cout << " A = " << iA << " Z = " << iZ << " J= " << fJ <<
        " Ex= " << fEx << " U= " << fU0 << endl; 

    }
  else
    {
      for(int i=1;i<iStore;i++) store[i].gamma += store[i-1].gamma;


    double xran = ran.Rndm();
    int i = 0;
    for (;;)
      {
        double prob = store[i].gamma/store[iStore-1].gamma;
        if (prob > xran) break;
        if (i == iStore-1) break;
        i++; 
      }
    fissionZ = store[i].iZ;
    fissionA = store[i].iA;
    }



  iZ1 = fissionZ;
  iZ2 = iZ - iZ1;
  iA1 = fissionA;
  iA2 = iA - iA1;


}

//************************************************************
/**
 * Treats the saddle to scission evaporations
 */
void CNucleus::saddleToScission()
{


  EvapEx2 = 0.;
  float widthEvaporation = evaporationWidthSS();
  float timeDecay = ran.expDecayTime(widthEvaporation);
  float newTime;
  if (widthEvaporation == 0) newTime = 100000.;
  else newTime = timeSinceSaddle + timeDecay;


  if (newTime < timeScission && EvapEx2 > 15.)
    {
      daughterLight = new CNucleus(EvapZ1,EvapA1);
      daughterLight->origin = 1;
      daughterLight->iWeight = 0;
      daughterLight->runningWeight = runningWeight;
      daughterLight->fact = fact;
 
      daughterLight->timeSinceStart = timeSinceStart + timeDecay;
      daughterLight->timeSinceSaddle = newTime;
      daughterHeavy = new CNucleus(EvapZ2,EvapA2);
      daughterHeavy->timeSinceSaddle = newTime; 
      daughterHeavy->origin = 1;
      daughterHeavy->iWeight = 0;
      daughterHeavy->runningWeight = runningWeight;
      daughterHeavy->fact = fact;

      daughterHeavy->timeSinceStart = daughterLight->timeSinceStart;
      daughterHeavy->timeScission = timeScission;
      daughterLight->excite(EvapEx1);

      daughterHeavy->sigma2 = sigma2;
      daughterHeavy->fissionZ = fissionZ;
      daughterHeavy->fissionA = fissionA;
      daughterHeavy->fissioningZ = fissioningZ;
      daughterHeavy->fissioningA = fissioningA;
      daughterHeavy->exciteScission(EvapEx2,fJ,fissionMassScission);

      daughterLight->parent = this;
      daughterHeavy->parent = this;

      if (evap.decay[EvapMode].Ek > 0)
        {
          daughterLight->notStatistical = 1;
          daughterLight->notStatisticalMode = EvapMode;
          daughterLight->bStable = 0;
        }
      angleIsotropic(); // find the emission angles of the two fragments
      CAngle angle(0.,0.);
      daughterLight->setSpinAxis(angle); 
      daughterHeavy->setSpinAxis(spin);
      allProducts[iProducts] = daughterLight;
      iProducts++;
     if (iProducts == Nproducts)
       {
         cout << "increase Nproducts" << endl;
	 abort();
       }

      allProducts[iProducts] = daughterHeavy;
      iProducts++;
     if (iProducts == Nproducts)
       {
         cout << "increase Nproducts" << endl;
	 abort();
       }

    }
  else 
    {

      if(fissionMassScission)massAsymmetry(Iscission);
      else
	{
          float Z = (float)fissionZ/(float)fissioningZ*(float)iZ;
          float A = (float)fissionA/(float)fissioningA*(float)iA;   
          fissionZ = (int)Z;
          fissionA = (int)A;
          if (ran.Rndm() < Z - (float)fissionZ) fissionZ++;
          if (ran.Rndm() < A - (float)fissionA) fissionA++;

	}
      
      // go back to the separation that gives us the correct
      // fission kinetic energy
      scission.sep = scission.sep0;
      float scissionE = scission.getScissionEnergy(fissionZ,fissionA);
      // for later us calculate fission kinetic energy and rotational
      // energy of each fragment
      scission.getFissionKineticEnergy(fissionZ,fissionA);


      scissionE -= fExpMass;
      fissionU = fEx - scissionE;
      fissionU = max((float)0.,fissionU);

    

      daughterLight = new CNucleus(fissionZ,fissionA);
      daughterLight->origin = 2;
      daughterLight->parent = this;
      daughterLight->iWeight = 0;
      daughterLight->runningWeight = runningWeight;
      daughterLight->fact = fact;


      daughterLight->timeSinceStart = timeSinceStart 
                    - timeSinceSaddle + timeScission;

      daughterHeavy = new CNucleus(iZ-fissionZ,iA-fissionA);
      daughterHeavy->origin = 3;
      daughterHeavy->parent = this;
      daughterHeavy->iWeight = 0;
      daughterHeavy->runningWeight = runningWeight;
      daughterHeavy->fact = fact;

      daughterHeavy->timeSinceStart = daughterLight->timeSinceStart;

      
      asyFissionDivide();

      allProducts[iProducts] = daughterLight;
      iProducts++;
     if (iProducts == Nproducts)
       {
         cout << "increase Nproducts" << endl;
	 abort();
       }

      allProducts[iProducts] = daughterHeavy;
      iProducts++;
     if (iProducts == Nproducts)
       {
         cout << "increase Nproducts" << endl;
	 abort();
       }


    }


}




//************************************************************
  /**
   * forces decay of 8Be
   */
void CNucleus::force8Be()
{
  daughterLight = new CNucleus(2,4);
  daughterHeavy = new CNucleus(2,4);
  daughterLight->origin = origin;
  daughterHeavy->origin = origin;
  daughterLight->parent = this;
  daughterHeavy->parent = this;
  daughterHeavy->notStatistical = false;
  daughterLight->notStatistical = false;
  daughterLight->excite(0.,0.);
  daughterHeavy->excite(0.,0.);
  daughterLight->bStable = 1;      
  daughterHeavy->bStable = 1;      

  float width;
  if (fEx < 1.)width = 5e-6; // use ground state experimental width
  else width = 1.53; //use first excited state experimental width
  float decayTime = timeSinceStart+ran.expDecayTime(width);
  daughterLight->timeSinceStart = decayTime;
  daughterHeavy->timeSinceStart = decayTime;

  EvapLPlusS = fJ;
  EvapL = (int)fJ;
  EvapS1 = 0.;
  EvapS2 = 0.;
  EvapEk = fEx + 0.0918; 
  EvapA1 = 4;
  EvapA2 = 4;
     
 
  angleEvap(); // find the emission angles of the two fragments
  allProducts[iProducts] = daughterLight;
  iProducts++;
  if (iProducts == Nproducts)
    {
      cout << "increase Nproducts" << endl;
      abort();
    }

  allProducts[iProducts] = daughterHeavy;
  iProducts++;
  if (iProducts == Nproducts)
    {
      cout << "increase Nproducts" << endl;
      abort();
    }

  return;
}          
//************************************************************
  /**
   * forces decay of 5Li
   */
void CNucleus::force5Li()
{
  daughterLight = new CNucleus(1,1);
  daughterHeavy = new CNucleus(2,4);
  daughterLight->origin = origin;
  daughterHeavy->origin = origin;
  daughterLight->parent = this;
  daughterHeavy->parent = this;
  daughterHeavy->notStatistical = false;
  daughterLight->notStatistical = false;
  daughterLight->excite(0.,0.5);
  daughterHeavy->excite(0.,0.);
  daughterLight->bStable = 1;      
  daughterHeavy->bStable = 1;

  float const width = 1.23;
  float decayTime = timeSinceStart+ran.expDecayTime(width);
  daughterLight->timeSinceStart = decayTime;
  daughterHeavy->timeSinceStart = decayTime;

      
  EvapLPlusS = fJ;
  EvapL = (int)(fJ-0.5);
  EvapS1 = 0.;
  EvapS2 = 0.;
  EvapEk = fEx + 1.69; 
  EvapA1 = 1;
  EvapA2 = 4;
     

  angleEvap(); // find the emission angles of the two fragments
  allProducts[iProducts] = daughterLight;
  iProducts++;
  if (iProducts == Nproducts)
    {
      cout << "increase Nproducts" << endl;
      abort();
    }

  allProducts[iProducts] = daughterHeavy;
  iProducts++;
  if (iProducts == Nproducts)
    {
      cout << "increase Nproducts" << endl;
      abort();
    }
  return;
}          
//************************************************************
  /**
   * forces decay of 5He
   */
void CNucleus::force5He()
{
  daughterLight = new CNucleus(0,1);
  daughterHeavy = new CNucleus(2,4);
  daughterLight->origin = origin;
  daughterHeavy->origin = origin;
  daughterLight->parent = this;
  daughterHeavy->parent = this;
  daughterHeavy->notStatistical = false;
  daughterLight->notStatistical = false;
  daughterLight->excite(0.,0.5);
  daughterHeavy->excite(0.,0.);
  daughterLight->bStable = 1;      
  daughterHeavy->bStable = 1;      

  float const width = .648;
  float decayTime = timeSinceStart+ran.expDecayTime(width);
  daughterLight->timeSinceStart = decayTime;
  daughterHeavy->timeSinceStart = decayTime;


  EvapLPlusS = fJ;
  EvapL = (int)(fJ-0.5);
  EvapS1 = 0.;
  EvapS2 = 0.;
  EvapEk = fEx + 0.798; 
  EvapA1 = 1;
  EvapA2 = 4;
     

  angleEvap(); // find the emission angles of the two fragments
  allProducts[iProducts] = daughterLight;
  iProducts++;
  if (iProducts == Nproducts)
    {
      cout << "increase Nproducts" << endl;
      abort();
    }

  allProducts[iProducts] = daughterHeavy;
  iProducts++;
  if (iProducts == Nproducts)
    {
      cout << "increase Nproducts" << endl;
      abort();
    }
  return;
}          
//************************************************************
  /**
   * forces decay of 9B
   */
void CNucleus::force9B()
{
  daughterLight = new CNucleus(1,1);
  daughterHeavy = new CNucleus(4,8);
  daughterLight->origin = origin;
  daughterHeavy->origin = origin;
  daughterLight->excite(0.,0.5);
  daughterHeavy->excite(0.,0.);
  daughterLight->bStable = 1;      
  daughterHeavy->bStable = 1;  
  daughterLight->parent = this;
  daughterHeavy->parent = this;    
  daughterHeavy->notStatistical = false;
  daughterLight->notStatistical = false;

  float const width = .54E-3;
  float decayTime = timeSinceStart+ran.expDecayTime(width);
  daughterLight->timeSinceStart = decayTime;
  daughterHeavy->timeSinceStart = decayTime;

  EvapLPlusS = fJ;
  EvapL = (int)(fJ-0.5);
  EvapS1 = 0.;
  EvapS2 = 0.;
  EvapEk = fEx + 0.1851; 
  EvapA1 = 1;
  EvapA2 = 8;
     

  angleEvap(); // find the emission angles of the two fragments
  allProducts[iProducts] = daughterLight;
  iProducts++;
  if (iProducts == Nproducts)
    {
      cout << "increase Nproducts" << endl;
      abort();
    }

  allProducts[iProducts] = daughterHeavy;
  iProducts++;
  if (iProducts == Nproducts)
    {
      cout << "increase Nproducts" << endl;
      abort();
    }
  return;
}          
//**********************************************************
/**
 * recursive function does multiple binary decays until excitation energy 
 * is exhausted.
 * 
 * After executation, the pointer array allProducts- points to each of the
 *intermediate and final products produced. The array stableProducts points
 * to just the final stable products. This array can be accessed to get 
 * these fragments 
 */

void CNucleus:: recursiveDecay()
{

  if (!bStable)
    {
     if (saddleToSciss) saddleToScission();
     else binaryDecay();
  
     if (abortEvent) return;
     }


  if (bStable) {

   //first force the decay of fragments that are unstable in their ground state
    if (iZ==4 && iA == 8) force8Be();
    else if (iZ==5 && iA == 9) force9B();
    else if (iZ == 3 && iA == 5) force5Li();
    else if (iZ == 2 && iA == 5) force5He();
    else
      {


        stableProducts[iStable]=this;
        iStable++;
       if (iStable == Nstable)
         {
           cout << "increase Nstable" << endl;
	   abort();
         }
  
        return;
      }
  }

  //daughterLight->print();
  if (daughterLight != NULL) 
    {
      daughterLight->recursiveDecay();
      if (daughterLight->abortEvent) 
         {
           abortEvent = 1;
           return;
         }
    }

  //daughterHeavy->print();
  //cout << daughterHeavy->iZ << " " << daughterHeavy->iA << " " << daughterHeavy->spin.theta << endl;
  //if (daughterHeavy->iA > 180) cout <<  daughterHeavy->spin.theta << endl;

  daughterHeavy->recursiveDecay();

  if (daughterHeavy->abortEvent)
    {
      abortEvent = 1;
      return;
    }

}
//**********************************************************
 /**
 * set the pointers to the static arrays of products to NULL.
 * 
 */
void CNucleus::resetGlobal()
{

  for (int i = iProducts-1;i>=0;i--)
    {
      allProducts[i]->daughterLight = NULL;
      allProducts[i]->daughterHeavy = NULL;
      delete allProducts[i];
      allProducts[i] = NULL;
    }
  for (int i=0;i<iStable;i++) stableProducts[i] = NULL;
  iProducts = 0;
  iStable = 0;
  iPoint = -1;
}
//*********************************************************
/**
 * This reset function should be used before starting another statistical
 * decay.
 */
void CNucleus::reset()
{
  //local reset
  resetGlobal();
  daughterLight = NULL;
  daughterHeavy = NULL;
  bStable = 0;

}
//***************************************************
/**
 * Prints out the information on all the stable decay products
 * produced in the statistical decay
 */

void CNucleus::printStableProducts()
{

  for (int i=0;i<iStable;i++)
    {
     stableProducts[i]->print();
    }
}
//*******************************************************************
/**
 * Prints out information on all products (stable and intermediates)
 * formed in the statistical decay
 */

void CNucleus::printAllProducts()
{
// prints out the information on qll of the  decay products
  for (int i=0;i<iProducts;i++)
    {
     allProducts[i]->print();
    }
}
//****************************************************************
/**
 * sets the excitation energy and spin of the compound nucleus.
 \param fEx0 is the excitation energy in MeV
 \param dJ0 is the spin in units of hbar
 */ 
void CNucleus::excite(float fEx0, double dJ0)
{
  excite(fEx0,(float)dJ0);
}
//****************************************************************
/**
 * sets the excitation energy and spin of the compound nucleus.
 \param dEx0 is the excitation energy in MeV
 \param dJ0 is the spin in units of hbar
 */ 
void CNucleus::excite(double dEx0, double dJ0)
{
  excite((float)dEx0,(float)dJ0);
}
//****************************************************************
/**
 * sets the excitation energy and spin of the compound nucleus.
 \param dEx0 is the excitation energy in MeV
 \param fJ0 is the spin in units of hbar
 */ 
void CNucleus::excite(double dEx0, float fJ0)
{
  excite((float)dEx0,fJ0);
}
//*****************************************************************
/**
 * sets the excitation energy and spin of the compound nucleus.
 \param fEx0 is the excitation energy in MeV
 \param fJ0 is the spin in units of hbar
 */ 
void CNucleus::excite(float fEx0, float fJ0)
{
  //initialized the excitation of the nucleus
  //calculates the level density
  

  notStatistical = 0;
  fEx = fEx0;
  //make sure we have either interg spin for even mass or half interger
  // spin for odd mass
  fJ = floor(fJ0);
  if (iA%2 == 1) fJ += 0.5;

  if (fEx == 0.) 
    {
      bStable = 1;
      return;
    }
  EdefScission = 0.;
  fPairing = mass->getPairing(iZ,iA); 
  fShell = mass->getShellCorrection(iZ,iA);  
  //fShell = fExpMass - mass.getLiquidDropMass(iZ,iA) - fPairing;
  if (iZ <= Zshell)
     {
      fPairing = 0.;
      fShell = 0.;
     }

  //note1
  //old code with standard SM treatment of yrast line
   Erot=yrast.getYrast(iZ,iA,fJ);
  //new code
  //spheroid.init(iZ,iA);
  //deformation = spheroid.getDeformationProlate(fJ);
  //Erot = spheroid.getRotatePlusDefEnergy(fJ);
  //end change

  fMInertia =  0.4*pow(r0,2)*pow((float)iA,(float)(5./3.));



  fU0 = fEx - Erot;
  if (fU0 < 0.) 
    {
      bStable = 1;
      return;
    }

  //decide if to use Hauser-Feshback calculation
  HF = 0;
  if ( (Erot > fU0/4. && fJ > 10) || iHF == 1) HF = 1;

  logLevelDensity = levelDensity.getLogLevelDensitySpherical(
		   	     iA,fU0,fPairing,fShell,fJ,fMInertia);

  if (fU0 < 0. || fU0 > 2000.)
    {
      cout << "Ex= " << fEx << "Erot= " << Erot << " iZ = " << iZ
	   << " iA= " << iA << " excite" << endl;
    }

  temp = levelDensity.getTemp();

  if (logLevelDensity == 0.)
    {
      bStable = 1;
      return;
    }
   
}
//*****************************************************************
/**
 * sets the excitation energy and of the compound nucleus
 * for a Weisskopf calculation (Spin not considered).
 \param fEx0 is the excitation energy in MeV
 */ 
void CNucleus::excite(float fEx0)
{
  //initialized the excitation of the nucleus
  //calculates the level density


  HF = 0;
  notStatistical = 0;
  fEx = fEx0;
  fJ = 0.;

  fU0 = fEx; 
  if (fU0 < 0.) 
    {
      bStable = 1;
      return;
    }

  EdefScission = 0.;
  fPairing = mass->getPairing(iZ,iA);   
  fShell = mass->getShellCorrection(iZ,iA);  
  //fShell = fExpMass - mass.getLiquidDropMass(iZ,iA) - fPairing;
  if (iZ <= Zshell)
     {
      fPairing = 0.;
      fShell = 0.;
     }
  Erot = 0.;



  logLevelDensity = 
      levelDensity.getLogLevelDensitySpherical(iA,fU0,fPairing,fShell);

  if (fU0 < 0. || fU0 > 2000.)
    {
      cout << "Ex= " << fEx << "Erot= " << Erot << " iZ = " << iZ
	   << " iA= " << iA << " excite " << endl;
    }

  temp = levelDensity.getTemp();

  if (logLevelDensity == 0.)
    {
      bStable = 1;
      return;
    }
   
}
//**********************************************************************
/**
 *  Initializes the excitation of the nucleus at its scission point
 * and calculates the level density
 */
void CNucleus::exciteScission(float fEx0,float fJ0,bool sym/*=1*/)
{

  saddleToSciss = 1;
  notStatistical = 0;
  HF = 0;
  fEx = fEx0;
  //make sure we have either interg spin for even mass or half interger
  // spin for odd mass
  fJ = floor(fJ0);
  if (iA%2 == 1) fJ += 0.5;

  if (sym)
    {
      scission.init(iZ,iA,fJ,2);
     EdefScission = scission.getScissionEnergy();
    }
  else
    {
     float Z1 = (float)fissionZ/(float)fissioningZ*(float)iZ;
     float A1 = (float)fissionA/(float)fissioningA*(float)iA;
 
     scission.init(iZ,iA,fJ,2,Z1,A1);

     EdefScission = scission.getScissionEnergy();   
    }

  EdefScission -= fExpMass;
  fPairing = 0.;
  fShell = 0.;
  Erot = 0.;
  if (fEx -EdefScission <= 0.) 
    {
       bStable = 1;
       return;
    }
  //thermal excitation energy above the scission-point configuration
  fU0 = fEx - Erot - EdefScission; 
  if (fU0 < 0.) 
    {
      bStable = 1;
      return;
    }

  logLevelDensity = levelDensity.getLogLevelDensityScission(iA,fU0);

  if (fU0 < 0 || fU0 > 2000)
    {
      cout << "fU0 = " << fU0 << " fEx= " << fEx << " Erot= " << Erot <<
	"EdefScission= " << EdefScission << " iZ = " << iZ << "iA= " << 
        iA << " fJ= " << fJ << " exciteScission " << endl;
    }

  temp = levelDensity.getTemp();

  if (logLevelDensity == 0.)
    {
      bStable = 1;
      return;
    }
   
}
//**********************************************************************
/**
 *  Returns the transition-state decay width.
 \param saddlePoint is the saddle-point energy in MeV
 \param iAfAn indicates that the saddle-point ld parameter is increased by afan
*/
/*
float CNucleus::getWidthZA(float saddlePoint,short iAfAn)
{

  if (saddlePoint > fEx) return 0.;
  float U = fEx - saddlePoint;
  float saddleLD = levelDensity.
           getLogLevelDensitySpherical(iA,U,(float)0.,(float)0.,
				       fJ,fMInertia,iAfAn);
  if (saddleLD == 0.) return 0.;
  if (saddleLD - logLevelDensity < -65.) return 0.;
  float gamma = exp(saddleLD-logLevelDensity)
                    *levelDensity.getTemp()/2./pi;
  return gamma;
}
*/
float CNucleus::getWidthZA(float saddlePoint,short iAfAn)
{

  if (saddlePoint > fEx) return 0.;
  float U = fEx - saddlePoint;
  float const deltaEk = 0.2;
  float gamma = 0.;
  float gammaMax = 0.;
  for(;;)
    { 
      float saddleLD = levelDensity.
           getLogLevelDensitySpherical(iA,U,(float)0.,(float)0.,
				       fJ,fMInertia,iAfAn);
     if (saddleLD == 0.) break;
     if (saddleLD - logLevelDensity < -65.) break;
     float gamma0 = exp(saddleLD-logLevelDensity)/2./pi;
     gammaMax = max(gammaMax,gamma0);
     gamma += gamma0*deltaEk;
     if (gamma0 < gammaMax/200.) break;
     
     U -= deltaEk;
     if (U <= 0.) break;
    }
  return gamma;
}

//*************************************************************************
/**
 * Calculates the asymmetric decay width in MeV from the gammaZ formalism.
 *
 * if noSymmetry = 1 then it only includes channels outside of the fission peak
 */

float CNucleus::asyFissionWidth()
{
  short iAfAn = 2;
  int const nStore =1000;
  SStore store[1000];
  needSymmetricFission = 0;
  int iStore = 0;
  float gammaTot = 0;
  float saddlePointOld = -10000.;
  yrast.prepareAsyBarrier(iZ,iA,fJ);
  //float Wigner0 = yrast.WignerEnergy(iZ,iA);
  scission.init(iZ,iA,fJ,1);

  for (int iZ1=evap.maxZ+1;iZ1<=iZ/2;iZ1++)
    {
      int iZ2 = iZ - iZ1;
      float A1 = (float)iZ1/(float)iZ*float(iA);
      float saddlePoint = yrast.getSaddlePointEnergy(A1);
      saddlePoint = (saddlePoint-Erot)*scaleImf + Erot;
      //add on a Wigner energy
      //saddlePoint += WignerScaled*Wigner0 + WignerAdd;
      saddlePoint +=  WignerAdd;
      if (iZ > Zshell) saddlePoint -= fPairing + fShell;

	if (noSymmetry && saddlePoint < saddlePointOld && iZ1 > 8)
	{
	  if (iZ/2 - iZ1 > 5 && iA > 120) 
	    {
              needSymmetricFission = 1;
              break;
	    }
	}
      else saddlePointOld = saddlePoint;

      float gammaZ0 = getWidthZA(saddlePoint,iAfAn);
      if (iZ1 == iZ - iZ1) gammaZ0 /= 2.;
      if (gammaZ0 <= 0.) continue;
      int iaMin1 = mass->chart.getAmin(iZ1);  
      int iaMin2 = iA - mass->chart.getAmax(iZ-iZ1);   
      int iaMin = max(iaMin1,iaMin2);
      int iaMax1 = mass->chart.getAmax(iZ1);   
      int iaMax2 = iA - mass->chart.getAmin(iZ-iZ1);  
      int iaMax = min(iaMax1,iaMax2);
      if (iZ1 == iZ - iZ1) iaMax = min(iA/2,iaMax);
      float gammaZ = 0;
      int iStart = iStore;
      for (int iA1 = iaMin;iA1<=iaMax;iA1++)
	{
          int iA2 = iA - iA1;
	  //make sure I can conserve energy latter
          float mass1 = mass->getCalMass(iZ1,iA1);
          float mass2 = mass->getCalMass(iZ2,iA2);
          scission.getFissionKineticEnergy(iZ1,iA1);  
	  float Qvalue = fExpMass - mass1 - mass2 - scission.ekTot 
          - scission.Erotate1 - scission.Erotate2;
          if( fEx + Qvalue <= 0.) continue;  

          float saddlePoint = yrast.getSaddlePointEnergy(iZ1,iA1);
	  saddlePoint = (saddlePoint-Erot)*scaleImf + Erot; 
	  //add on a congruence energy


          //float Wigner1 = yrast.WignerEnergy(iZ1,iA1);
          //float Wigner2 = yrast.WignerEnergy(iZ2,iA2);
   
          //saddlePoint += WignerScaled*(Wigner1+Wigner2-Wigner0) + WignerAdd;
          saddlePoint +=  WignerAdd;
          if (iZ > Zshell) saddlePoint -= fPairing + fShell; 


          float gamma = getWidthZA(saddlePoint,iAfAn);
          if (iZ1 == iZ - iZ1 && iA1 == iA - iA1) gamma /= 2.;
          if (gamma <= 0.) continue;
          store[iStore].gamma = gamma;
          store[iStore].iZ = iZ1;
          store[iStore].iA = iA1;

          iStore++;
	  if (iStore == nStore)
	    {
              cout << "increase nStore in asyFissionWidth" << endl;
	      abort();
	    }
          gammaZ += gamma;
	}

      // normalize gamma
     
      if (gammaZ == 0.) iStore=iStart;
      else 
	{
	  for (int j=iStart;j<iStore;j++) store[j].gamma *= gammaZ0/gammaZ;
	  gammaTot += gammaZ0;

	}
     
    }

  if (iStore == 0) return 0.;

  for (int j=0;j<iStore;j++)
    {
      store[j].gamma /= gammaTot;
      if (j > 0) store[j].gamma += store[j-1].gamma;
    }


  //call random nummber
  float x = ran.Rndm();
  int j=0;
  for (;;)
    {
      if (x < store[j].gamma) break;
      if (j == iStore - 1) break;
      j++;
      
    }
  fissionZ = store[j].iZ;
  fissionA = store[j].iA;
  float saddlePoint = yrast.getSaddlePointEnergy(fissionZ,fissionA);
  saddlePoint = (saddlePoint-Erot)*scaleImf + Erot; 

  /*
  //add on a congruence energy
  int iA2 = iA - fissionA;
  int iZ2 = iZ - fissionZ;
  float Wigner1 = yrast.WignerEnergy(fissionZ,fissionA);
  float Wigner2 = yrast.WignerEnergy(iZ2,iA2); 
  saddlePoint += WignerScaled*(Wigner1+Wigner2-Wigner0) + WignerAdd;
  */
  saddlePoint +=  WignerAdd;
  if (iZ > Zshell) saddlePoint -= fPairing + fShell;
  fissionU = fEx - saddlePoint;

  return gammaTot;
}
//********************************************************************
/**
 * Calculates the complex fragments decay widths in MeV where the total fission
 * width all channels is normalised to the Bohr-Wheeler result.
 */
float CNucleus::asyFissionWidthBW()
{


  short  iAfAn = 1;
  int const nnstore=4000;
  SStore store[nnstore];
  needSymmetricFission = 0;
  int iStore = 0;
  int iFission = -1;
  float gammaTot = 0;
  float saddlePointOld = -10000.;
  yrast.prepareAsyBarrier(iZ,iA,fJ);
  float A1 = (float)iA/2.;


  //the following is the differece between the symmetric saddle points 
  //Sierk symmetric calculation and my interpolated and extrapolated
  //asymmetric barriers using Sierk asymmetry calculations as input
  // at J=0, they are garenteed to be identical
  // I fill just shift down the asymmetric barrirs by delta now

  float delta = yrast.getSaddlePointEnergy(A1) - symSaddlePoint;

  for (int iZ1=evap.maxZ+1;iZ1<=iZ/2;iZ1++)
    {

      float A1 = (float)iZ1/(float)iZ*float(iA);
      float saddlePoint = yrast.getSaddlePointEnergy(A1) - delta;
      if (iZ > Zshell) saddlePoint -= fPairing + fShell; 

      float A2 = (float)(iZ1-1)/(float)(iZ)*float(iA);
      float A3 = (float)(iZ1+1)/(float)(iZ)*float(iA);
      float A4 = (float)(iZ1-2)/(float)(iZ)*float(iA);
      float A5 = (float)(iZ1+2)/(float)(iZ)*float(iA);
      float saddlePoint2 = yrast.getSaddlePointEnergy(A2) - delta;
      float saddlePoint3 = yrast.getSaddlePointEnergy(A3) - delta;
      float saddlePoint4 = yrast.getSaddlePointEnergy(A4) - delta;
      float saddlePoint5 = yrast.getSaddlePointEnergy(A5) - delta;

      saddlePoint = (saddlePoint+saddlePoint2+saddlePoint3
                    +saddlePoint4+saddlePoint5)/5.;


      if (saddlePoint < saddlePointOld && iZ1 > 8 && iFission == -1)
	{
	  if (iZ/2 - iZ1 > 5) 
	    {
              iFission = iStore;
               
	    }
	}
      else saddlePointOld = saddlePoint;

      


      float gammaZ0 = getWidthZA(saddlePoint,iAfAn);
      if (iZ1 == iZ - iZ1) gammaZ0 /= 2.;
      if (gammaZ0 <= 0.) continue;
      int iaMin1 = mass->chart.getAmin(iZ1);  
      int iaMin2 = iA - mass->chart.getAmax(iZ-iZ1);  
      int iaMin = max(iaMin1,iaMin2);
      int iaMax1 = mass->chart.getAmax(iZ1);  
      int iaMax2 = iA - mass->chart.getAmin(iZ-iZ1);  
      int iaMax = min(iaMax1,iaMax2);
      if (iZ1 == iZ - iZ1) iaMax = min(iA/2,iaMax);
      float gammaZ = 0;
      int iStart = iStore;
      for (int iA1 = iaMin;iA1<=iaMax;iA1++)
	{

          float saddlePoint = yrast.getSaddlePointEnergy(iZ1,iA1);
          if (iZ > Zshell) saddlePoint -= fPairing + fShell; 

          float gamma = getWidthZA(saddlePoint,iAfAn);
          if (iZ1 == iZ - iZ1 && iA1 == iA - iA1) gamma /= 2.;
          if (gamma <= 0.) continue;
          store[iStore].gamma = gamma;
          store[iStore].iZ = iZ1;
          store[iStore].iA = iA1;

          iStore++;
          if (iStore == nnstore) abort();
          gammaZ += gamma;

	}

      // normalize gamma

      if (gammaZ == 0.) iStore=iStart;
      else 
	{
	  for (int j=iStart;j<iStore;j++) store[j].gamma *= gammaZ0/gammaZ;
	  gammaTot += gammaZ0;

	}
    }


  if (iStore == 0) 
    {
     fissionZ = iZ/2;
     fissionA = iA/2;
     fissionU = 0.;
     return 0.;
    }
  


  // to normalise fission peak to Bohr-wheeler width
  if (iFission > 0)
    {
      float gammaFission = 0.;
      for (int i=iFission;i<iStore;i++)
	{
	  gammaFission += store[i].gamma;
	}
      float gammaBW = BohrWheelerWidth();
      float ratio = gammaBW/gammaFission;
      gammaTot = gammaBW;
      for (int i=iFission;i<iStore;i++)
	{
          store[i].gamma *= ratio;
	}

      //set other imf's to zero
      for (int i=0;i<iFission;i++)
	{
          store[i].gamma = 0.;
	}
      
    }

  for (int j=0;j<iStore;j++)
    {
      store[j].gamma /= gammaTot;
      if (j > 0) store[j].gamma += store[j-1].gamma;
    }


  //call random nummber
  float x = ran.Rndm();
  int j=0;
  for (;;)
    {
      if (x < store[j].gamma) break;
      if (j == iStore - 1) break;
      j++;
    }
  fissionZ = store[j].iZ;
  fissionA = store[j].iA;
  float saddlePoint = yrast.getSaddlePointEnergy(fissionZ,fissionA);
  if (iZ > Zshell) saddlePoint -= fPairing + fShell;
  fissionU = fEx - saddlePoint;


  return gammaTot;
}

//*************************************************************************
/**
 * Uses the GammaZA formalism to get complex fragment decay widths as in the 
 * original GEMINI.
 *
 * The asymmetric fission width are in MeV.
 */
float CNucleus::asyFissionWidthZA()
{
  short iAfAn = 2;
  SStore store[nStore];
  needSymmetricFission = 0;
  int iStore = 0;
  float gammaTot = 0;
  float saddlePointOld = -10000.;
  yrast.prepareAsyBarrier(iZ,iA,fJ);
  //yrast.printAsyBarrier();
  //float Wigner0 = yrast.WignerEnergy(iZ,iA);
  scission.init(iZ,iA,fJ,1);

  for (int iZ1=evap.maxZ+1;iZ1<=iZ/2;iZ1++)
    {
      int iZ2 = iZ - iZ1;
      float A1 = (float)iZ1/(float)iZ*float(iA);
      float saddlePoint = yrast.getSaddlePointEnergy(A1);
      saddlePoint = (saddlePoint-Erot)*scaleImf + Erot; 
      //add on a congruence energy

      //saddlePoint += WignerScaled*Wigner0+WignerAdd;
      saddlePoint += WignerAdd;

      if (iZ > Zshell) saddlePoint -= fPairing + fShell; 


      
      if (noSymmetry && saddlePoint < saddlePointOld && iZ1 > 8
	  && iZ/2 - iZ1 > 5)
	{
	  needSymmetricFission = 1.;
	  break;
	}
      else saddlePointOld = saddlePoint;
      

      int iaMin1 = mass->chart.getAmin(iZ1);  
      int iaMin2 = iA - mass->chart.getAmax(iZ-iZ1);  
      int iaMin = max(iaMin1,iaMin2);
      int iaMax1 = mass->chart.getAmax(iZ1);  
      int iaMax2 = iA - mass->chart.getAmin(iZ-iZ1);   
      int iaMax = min(iaMax1,iaMax2);
      if (iZ1 == iZ - iZ1) iaMax = min(iA/2,iaMax);

      for (int iA1 = iaMin;iA1<=iaMax;iA1++)
	{

          int iA2 = iA - iA1;
	  //make sure I can conserve energy latter
          float mass1 = mass->getCalMass(iZ1,iA1);
          float mass2 = mass->getCalMass(iZ2,iA2);
          scission.getFissionKineticEnergy(iZ1,iA1);  
	  float Qvalue = fExpMass - mass1 - mass2 - scission.ekTot 
          - scission.Erotate1 - scission.Erotate2;
          if( fEx + Qvalue <= 0.) continue;  

          float saddlePoint = yrast.getSaddlePointEnergy(iZ1,iA1);
          saddlePoint = (saddlePoint-Erot)*scaleImf + Erot; 
          //add on a congruence energy
	  /*
          float Wigner1 = yrast.WignerEnergy(iZ1,iA1);
          float Wigner2 = yrast.WignerEnergy(iZ2,iA2);
          saddlePoint += WignerScaled*(Wigner1+Wigner2-Wigner0)+WignerAdd;
	  */
          saddlePoint += WignerAdd;
          if (iZ > Zshell) saddlePoint -= fPairing + fShell; 



          float gamma = getWidthZA(saddlePoint,iAfAn);
 
          if (iZ1 == iZ - iZ1 && iA1 == iA - iA1) gamma /= 2.;
          if (gamma <= 0.) continue;

	  /*
          // the following paragraph of code corrects the Moretto 
          //decay width such that a Lestone like logic is used. 
          //This was tried at 
          //some point, but for symmetric fission we obtained better agreement
          //for fusion and splattion reactions with the standard 
          //Bohr-Wheeler value. It therefore was inconsistent to used 
          //the Lestone logical for asymmetric fission
          //
	  // start Lestone correction
          //float R1 = pow((float)iA1,(float)(1./3.))*r0;
          //float R2 = pow((float)(iA-iA1),(float)(1./3.))*r0;

          //float momInertia1 = 0.4*(float)iA*pow(R1,2);
          //float momInertia2 = 0.4*(float)(iA-iA1)*pow(R2,2);
          //float Areduced = (float)(iA1*(iA-iA1))/float(iA);
          //float momInertiaPerp = momInertia1 + momInertia2 + 
	  //  Areduced*(R1+R2+sep);
          //float momInertiaPara = momInertia1 + momInertia2;
	  //float momInertiaEff = 1./(1./momInertiaPara - 1./momInertiaPerp);
          //float corr = LestoneCorrection(fEx-saddlePoint,momInertiaEff,iAfAn);
          //gamma *= corr;
	  //end lestone correction
	  */          

          gammaTot += gamma;
          store[iStore].gamma = gamma;
          store[iStore].iZ = iZ1;
          store[iStore].iA = iA1;

          iStore++;
	  if (iStore == nStore)
	    {
              cout << "increase nStore in asyFissionWidthZA" << endl;
	      abort();
	    }
	}
    }

  if (iStore == 0) return 0.;

  for (int j=0;j<iStore;j++)
    {
      store[j].gamma /= gammaTot;
      if (j > 0) store[j].gamma += store[j-1].gamma;
    }


  //call random nummber
  float x = ran.Rndm();
  int j=0;
  for (;;)
    {
      if (x < store[j].gamma) break;
      if (j == iStore -1) break;
      j++;
    }
  fissionZ = store[j].iZ;
  fissionA = store[j].iA;
  float saddlePoint = yrast.getSaddlePointEnergy(fissionZ,fissionA);
  saddlePoint = (saddlePoint-Erot)*scaleImf + Erot; 
  //add on a congruence energy
  /*
  int iA2 = iA - fissionA;
  int iZ2 = iZ - fissionZ;
  float Wigner1 = yrast.WignerEnergy(fissionZ,fissionA);
  float Wigner2 = yrast.WignerEnergy(iZ2,iA2);
  saddlePoint += WignerScaled*(Wigner1+Wigner2-Wigner0) + WignerAdd;
  */

   
  saddlePoint +=  WignerAdd;
  if (iZ > Zshell) saddlePoint -= fPairing + fShell;
  fissionU = fEx - saddlePoint;

  return gammaTot;
}
//********************************************************************
float CNucleus::LestoneFissionWidth()
/**
  * Gives the fission decay width from Lestone in units of MeV.
  * Lestone gives an extention of the BohrWheeler width with the inclusion 
  * of one of the angular momentum degrees of freedom (tilting mode).
  * See PRC 59 (1999) 1540. 
  */
{
  float gamma = BohrWheelerWidth();
  if (gamma <= 0.) return gamma;

  short iAfAn = 1;

  float momInertiaPerp = yrast.getMomentOfInertiaSierk(fJ);
  float momInertiaPara = yrast.momInertiaMin;
  float momInertiaEff = 1.0/(1.0/momInertiaPara-1.0/momInertiaPerp);

  float U = fEx - symSaddlePoint; 

  float corr = LestoneCorrection(U,momInertiaEff,iAfAn);
  return gamma*corr;
}


//*************************************************************************
float CNucleus::BohrWheelerWidth()
/**
  * Gives the Bohr-Wheeler decay width for fission in units of MeV
  */
{

  symSaddlePoint = yrast.getSymmetricSaddleEnergy(iZ,iA,fJ) + barAdd;
  if (iZ > Zshell) symSaddlePoint -= fPairing + fShell; 

  /*float fJJ = fJ;
  if (fJ > yrast.Jmax) fJJ = yrast.Jmax;*/
  short iAfAn = 1;
  float gamma = getWidthZA(symSaddlePoint,iAfAn)*fissionScaleFactor;
  return gamma;
}
//************************************************
/**
 *  Responsible for subdividing spin and excitation energy for
 * asymmetry fission.
 *
 * its uses a two-sphere approximation to get collective modes
 */

void CNucleus::asyFissionDivide()
{




  //first a statement to handle a bad case
 if (fissionU <= 0.)
   {
     daughterLight->excite(0.,0.);
     daughterHeavy->excite(0.,0.);

     float phi= 2.*pi*ran.Rndm();
     float theta = acos(1.0-2.0*ran.Rndm());
     CAngle symmetryCM(theta,phi);
     split(symmetryCM);
     return;  
   }



  float A1 = (float)fissionA;
  float A2 = (float)(iA-fissionA);
  float Ared = A1*A2/(A1+A2);
  float U1 = fissionU*A1/(A1+A2);
  float U2 = fissionU - U1;
  float r1 = scission.r0*pow(A1,(float)(1./3.));
  float r2 = scission.r0*pow(A2,(float)(1./3.));
  float MInertia1 = 0.4*A1*pow(r1,2);
  float MInertia2 = 0.4*A2*pow(r2,2);

  float MInertiaOrbit = Ared*pow(r1+r2+scission.sep,2);
  float MInertia12 = MInertia1 + MInertia2;
  float MInertiaSaddle = MInertia12 + MInertiaOrbit;


  float aden1 = levelDensity.getLittleA(fissionA,U1);
  //float entropy1 = levelDensity.getEntropy();
  float aden2 = levelDensity.getLittleA(iA-fissionA,U2);
  //float entropy2 = levelDensity.getEntropy();
  float aden = aden1 + aden2;
  //float entropy0 = entropy1 + entropy2;

  float temp = sqrt(fissionU/aden);
  float entropy0 = 2.*sqrt(aden*fissionU);

  float fact10 =  MInertiaOrbit*(MInertia2-MInertia1)
    /(2.0*MInertia1*MInertia2);
  float fact1 = sqrt(pow(fact10,2)+1);
  float theta = atan(fact10+fact1);
  float fact3  = MInertia12/(MInertia1*MInertia2);
  float c = cos(theta);
  float s = sin(theta);
  float cc = pow(c,2);
  float ss = pow(s,2);
  float fact4 =  2.0*s*c/MInertiaOrbit;
  float aBending = (cc/MInertia1+ss/MInertia2
		    +1.0/MInertiaOrbit-fact4)*kRotate;
 float aWriggling = (ss/MInertia1+cc/MInertia2 
		    +1.0/MInertiaOrbit+fact4)*kRotate;
 float aTilting = MInertiaOrbit/MInertia12/MInertiaSaddle*kRotate;
 float aTwisting = fact3*kRotate;


 float K,JTwisting,JBendingX,JBendingZ,JWrigglingX,JWrigglingZ;
 float deltaE;



 int tries = 0;
 for (;;)
   {
     //choose K tilting
     K = selectJ(aTilting,aden,entropy0,fJ);

     //choose twisting
     float sig = 5.*sqrt(temp/aTwisting);
     JTwisting = selectJ(aTwisting,aden,entropy0,sig);


     //choose bending for X direction
     sig = 5.0*sqrt(temp/aBending);
     JBendingX = selectJ(aBending,aden,entropy0,sig);

     //choose Jbending for Z direction
     JBendingZ = selectJ(aBending,aden,entropy0,sig);


     //choose Jwriggling for x direction
     sig = 5.0*sqrt(temp/aWriggling);
     JWrigglingX=selectJ(aWriggling,aden,entropy0,sig);

     //choose Jwriggling for Z direction
     JWrigglingZ=selectJ(aWriggling,aden,entropy0,sig);


    deltaE = 0.5*(aTilting*pow(K,2) + aTwisting*pow(JTwisting,2)
    + aBending*(pow(JBendingX,2)+pow(JBendingZ,2))
	      + aWriggling*(pow(JWrigglingX,2)+pow(JWrigglingZ,2)));

    if (deltaE <= fissionU) break;

    tries++;
    if (tries == 15) 
      {
        // not able to find a solution for some reason.
	//so just turn off the depolarization modes.
        //cout << "tries too big in asyFissionDivide" << endl;
        //abort();
        K = 0.;
	JTwisting = 0.;
        JBendingX = 0.;
        JBendingZ = 0.;
        JWrigglingX = 0.;
        JWrigglingZ = 0.;
        deltaE = 0.;
        break;
     }
   }


 //calculate components of r1, depolarizing mode of fragment 1
 float R1x = JBendingX*c + JWrigglingX*s;
 float R1y = JTwisting;
 float R1z = JBendingZ*c + JWrigglingZ*s;

 //       calculate components of r2
 float R2x = JWrigglingX*c - JBendingX*s;
 float R2y = -JTwisting;
 float R2z = JWrigglingZ*c - JBendingZ*s;


 float omega_t = ran.Rndm()*2.*pi;
 float J0xz = sqrt(pow(fJ,2)-pow(K,2));
 float J0x = J0xz*sin(omega_t);
 float J0z = J0xz*cos(omega_t);


 //calculate components of s1
 float ratio1 = MInertia1/MInertiaSaddle;
 float S1x = J0x*ratio1 + R1x;
 float S1y = MInertia1/MInertia12*K + R1y;
 float S1z = J0z*ratio1 + R1z;
 float S1 = sqrt(pow(S1x,2) + pow(S1y,2) + pow(S1z,2));

  //calculate components of s2
 float ratio2 = MInertia2/MInertiaSaddle;
 float S2x = J0x*ratio2 + R2x;
 float S2y = MInertia2/MInertia12*K + R2y;
 float S2z = J0z*ratio2 + R2z;
 float S2 = sqrt(pow(S2x,2) + pow(S2y,2) + pow(S2z,2));

 //note that fissionU was calculated as the thermal energy above
 //a rigidly rotation configuration. The rotational energy of the 
 //individual frgament's in that configuration has already been 
 //subtracted. As we now about to assign new rotational energies
 //to these fragment, this component must be corrected for

 //old rotational energies
 float Erotate10 = scission.Erotate1;
 float Erotate20 = scission.Erotate2;


 float Erotate1 = pow(S1,2)*kRotate/2./MInertia1;
 float Erotate2 = pow(S2,2)*kRotate/2./MInertia2; 
 fissionU += Erotate10 + Erotate20 - Erotate1 - Erotate2;
 if (fissionU < 0.) fissionU = 0.;

 temp = sqrt(fissionU/aden);
 float sig = sqrt(2.0*pow(temp,3)*aden1*aden2/aden);
 U1 = fissionU*A1/(A1+A2);
 U2 = fissionU - U1;
 entropy0 = 2.0*sqrt(aden*fissionU);
 float demax = min((float)5.0*sig,U1);
 float demin = min((float)5.0*sig,U2);

 tries = 0;
 for (;;) 
   {
     deltaE = (demin+demax)*ran.Rndm() - demax;
    float y = exp(2.0*sqrt(aden1*(U1+deltaE)) 
		  +2.0*sqrt(aden2*(U2-deltaE))-entropy0);
    float yy = ran.Rndm();
    // cout << deltaE << " " << y << " " << yy << endl;
    if (yy <= y) break;
    tries++;
    if (tries == 10)
      {
        //cout << "2nd tries too big in asyFissionDivide" << endl;
        deltaE = 0.;
        break;
      }
   }

 float Ex1 = U1 + deltaE + Erotate1;
 float Ex2 = U2 - deltaE + Erotate2;


 // the following should be zero if energy is conserved
 //cout << fEx + fExpMass - daughterLight->fExpMass - daughterHeavy->fExpMass
 //- scission.ekTot -Ex1 - Ex2 << endl;

 float fJ1,fJ2;
 if (fissionA%2 == 0) 
   {
     fJ1 = round(S1);
   }
 else
   {
     fJ1 = round(S1+.5) - .5;
   }
 daughterLight->excite(Ex1,fJ1);

 if ((iA-fissionA)%2 == 0) 
   {
     fJ2 = round(S2);
   }
 else
   {
     fJ2 = round(S2+.5) - .5;
   }
 daughterHeavy->excite(Ex2,fJ2);


 // this part of the subroutines calculates the final
 // centre-of-mass velocity vectors of the complex fragments and
 // their spin vectors.  vectors are specified as:
 //                x_y_z
 // x is polar angular coordinate for vector (theta or phi)
 // y  identifies the vector
 // options  :  V1 -  velocity of fragment 1
 //             S1 -  spin of fragment 1
 //             V2 -  velocity of fragment 2
 //             S2 -  spin of fragment 2
 //            symmetry - specifies the symmetry axis of
 //            saddle-scission point configuration
 //             J0 - spin vector of fragment 0 (decaying system)
 //             V0 - velocity vector of fragment 0
 // Z is frame in which polar coordinates are specified
 //             r - rotation frame of saddle point configuration
 //             see Schmitt and Pacheco Nucl Phys. A379 (1982) 313
 //             J0 - frame with Vector J0 as Z axis
 //             cm - centre of mass frame of reaction.



 // in the frame of the di-nuclear complex the total angular
 // momentum vector processes. The theta and phi angles at
 // scission are. see schmitt and pacheco for coordinate system.
 CAngle J0r;
 if (fJ == 0.0) 
  {
    J0r.phi= 2.*pi*ran.Rndm();
    J0r.theta = -acos(1.0-2.0*ran.Rndm());
  }
 else
  {
    J0r.theta = -acos(J0z/fJ);
    J0r.phi = atan2(K,J0x);
  }

 // the symmetry axis occurs at x=0,z=0 in the above coordinate
 // system i.e. theta_symmetry_r=90., phi_symmetry_r =90. degrees, find
 // angles of the symmetry axis relative to the j0 vector
 // (.i.e. rotate everything so that J0 is parallel to z axis)
 CAngle symmetryJ0;
 CAngle perp(pi/2.,pi/2.); //rjc
 symmetryJ0 = CAngle::transform(perp,J0r);


 // randomize the phi coordinate.
 float dphi = 2.*pi*ran.Rndm();
 symmetryJ0.phi += dphi;
 if (symmetryJ0.phi > 2.*pi) symmetryJ0.phi -= 2.*pi;


 // find symmetry axis in lab frame
 CAngle symmetryCM;
 symmetryCM = CAngle::transform(symmetryJ0,spin);

 // in the rotating frame of dinuclear system, the spin vector of
 // fragment 1 is s1X,s1Y,s1Z, relative to di nuclear coordinate system
 // mentioned above
 CAngle S1r;
 if (S1 == 0.)
   {
     S1r.theta = acos(1.0-2.0*ran.Rndm());
     S1r.phi = 2.*pi*ran.Rndm();
   }
 else
   {
     S1r.theta = acos(S1z/S1);
     S1r.phi = atan2(S1y,S1x);
   }

 // rotate to frame of J0 vector (.i.e. z axis parallel to J0)
 CAngle S1J0;
 S1J0 = CAngle::transform(S1r,J0r);

 S1J0.phi += dphi;
if (S1J0.phi > 2.*pi) S1J0.phi -= 2.*pi;
  // rotate to lab frame 
 daughterLight->spin = CAngle::transform(S1J0,spin);     

  //  do the same for fragment 2
 CAngle S2r;
 if (S2 == 0.0)
  {
    S2r.theta = acos(1.0-2.0*ran.Rndm());
    S2r.phi = 2.*pi*ran.Rndm();
  }
 else
   {
     S2r.theta = acos(S2z/S2);
     S2r.phi = atan2(S2y,S2x);
   }


 // rotate to frame of J0 vector (.i.e. z axis parallel to J0)
 CAngle S2J0;
 S2J0 = CAngle::transform(S2r,J0r);
 S2J0.phi += dphi;
 if (S2J0.phi > 2.*pi) S2J0.phi -= 2.*pi;

 // rotate to lab frame
 daughterHeavy->spin = CAngle::transform(S2J0,spin);      

  // calculate velocities in frame of emitting system
  // first find the kinetic energy of separation
  //Ek = 1.44*(float)fissionZ*(float)(iZ-fissionZ)/(r1+r2+sep)
  //generlization of Viola see Hinde Nucl Phys A472, 318 (1987).
 //float Ek = yrast.viola((float)fissionZ,A1,(float)(iZ-fissionZ),A2);
 //float Ek = 0.755*(float)fissionZ*(float)(iZ-fissionZ)/
 // (pow(A1,(float)(1./3.))+pow(A2,(float)(1./3.))) 
 //  + 7.3; 


 // find orbital angular momentum
 //float ratio_l = MInertiaOrbit/MInertiaSaddle;
 //float lx = J0x*ratio_l - R1x - R2x;
 //float lz = J0z*ratio_l - R1z - R2z;
 //float l = sqrt(pow(lx,2) + pow(lz,2));
 //ek = ek + 20.8993*l**2.0/M_inertia_orbit

 split(symmetryCM);

}
//************************************************
/**
 *  Randomly selects the spin associated with a fission normal modes such as
 * wriggling, twisting, etc.
 */
float CNucleus::selectJ(float ac, float aden, float entropy0,float Jmax)
{

  float Jmax1 = sqrt(2.0*fissionU/ac);
  float Jmax2 = min(Jmax,Jmax1);

  float J;
  int tries = 0;
  for (;;)
  {
    J = Jmax2*(1.0-2.0*ran.Rndm());
    float diff = max( fissionU-0.5*ac*J*J, 0.0 );
    float y1 = exp(2.0*sqrt(aden*diff)-entropy0);
    float y2 = ran.Rndm();
    if (y2 < y1) break;
    tries++;
    if (tries == 15)
      {
        //after 15 tries, return zero
	//cout << "tries too big in selectJ" << endl;
        //abort();
        return 0.;
      }

  }
  
  return J;
}

//****************************************
/**
 *  Sets the angle of the compound nucleus spin axis. 
 *
  \param spin0 is the (theta,phi) angles in radians
*/
void CNucleus::setSpinAxis(CAngle spin0)
{

  spin = spin0;
}
//****************************************
/**
 *  Sets the angle of the compound nucleus spin axis. 
 *
  \param spin0 is the (theta,phi) angles in degrees
*/
void CNucleus::setSpinAxisDegrees(CAngle spin0)
{
  spin = spin0;
  spin.theta *= CAngle::pi/180.;
  spin.phi *= CAngle::pi/180.;
}
//***************************************
/**
 * Sets the velocity of the fragment in polar coordinates
 *
 /param vel is velocity in units of cm/ns
 /param theta is theta angle in radians
 /param phi is phi angle in radians
*/
void CNucleus::setVelocityPolar(float vel/*=0.*/, 
                                float theta/*=0.*/, float phi/*=0.*/)
{

  //default values are zero
  velocity[0] = vel*sin(theta)*cos(phi);
  velocity[1] = vel*sin(theta)*sin(phi);
  velocity[2] = vel*cos(theta);
} 
//***************************************
/**
 * Sets the velocity of the fragment in cartesian coordinates.
 *
 * with no input parameter, all components are set to zero.
 \param vx is x component of velocity in cm/ns
 \param  vy is y component of velocity in cm/ns
 \param vz is z component of velocity in cm/ns
 */
void CNucleus::setVelocityCartesian(float vx/*=0.*/, float vy/*=0.*/, 
                                    float vz/*=0.*/)
{
  //default values are zero
  velocity[0] = vx;
  velocity[1] = vy;
  velocity[2] = vz;
}
//*****************************************
/**
 *  Calculates the total decay widths in MeV for light particle evaporation.
 * using the Hauser-Feshbach formulism 
 */
float CNucleus::evaporationWidth()
{

  float width = 0.;
  EcostMin = 1000.;
  for (int i=0;i<evap.nLight;i++) 
    {
      lightP = evap.lightP[i];


      if (HF) evap.prob[i] = hauserFeshbach(i);
      else  evap.prob[i] = weiskopf((bool)0);

      evap.prob[i] *= lightP->suppress;//user-given suppression factor
      width += evap.prob[i];
      if (i > 0) evap.prob[i] += evap.prob[i-1];
    }

  if (width <= 0.) 
    {
      return 0.;
    }



  return width;
}
//*****************************************
/**
 * Calculated the Hauser-Feshbach decay width in MeV for a given channel
\param iChannel is the channel of the evaporated particle
*/
float CNucleus::hauserFeshbach(int iChannel)
{
  //returns the Hauser-Feshbach evaporation width for a given channel
  //width is in units of MeV
  if (exp(-lightP->fEx/temp) < 0.01  && iChannel > 1) return 0.;

  lightP->residue.iZ = iZ - lightP->iZ;
  if (lightP->iZ >= lightP->residue.iZ) return 0.;
  lightP->residue.iA = iA - lightP->iA;
  if (lightP->iA >= lightP->residue.iA) return 0.;
  if (lightP->residue.iZ >= lightP->residue.iA) return 0.;
  //products must lie on the chart of nuclides used by gemini
  if (lightP->residue.iA < mass->chart.getAmin(lightP->residue.iZ)) return 0.;
  if (lightP->residue.iA > mass->chart.getAmax(lightP->residue.iZ)) return 0.;

  //note1, new treament of yrast line 
  //spheroid.init(lightP->residue.iZ,lightP->residue.iA);
  //spheroid.setDeformationProlate(deformation);



  //odd A nuclei have 1/2 integer spin
  lightP->odd =0;
  if (lightP->residue.iA%2 == 1) lightP->odd = 1; 

  if (lightP->residue.iZ > Zshell)
  lightP->fPair = mass->getPairing(lightP->residue.iZ,lightP->residue.iA);
  else lightP->fPair = 0.;

  lightP->residue.fExpMass = mass->getCalMass(lightP->residue.iZ,lightP->residue.iA);  
  lightP->fShell = 0.;

  if (lightP->residue.iZ > Zshell) lightP->fShell = 
             mass->getShellCorrection(lightP->residue.iZ,lightP->residue.iA); 

  lightP->separationEnergy = lightP->residue.fExpMass + lightP->fExpMass - fExpMass;


  rResidue = r0*pow((float)lightP->residue.iA,(float)(1./3.));
  lightP->fMInertia = 0.4*(float)lightP->residue.iA*pow(rResidue,2);
  lightP->fMInertiaOrbit = (float)(lightP->residue.iA*lightP->iA)/
    (float)(lightP->residue.iA+lightP->iA)*pow(rResidue+lightP->rLight+2.,2);

  // decide if decay mode is too costly` (very rare) if so abort calculation
  Ecoul = 1.44*(float)(lightP->iZ*lightP->residue.iZ)/(rResidue+lightP->rLight+2.);
  float Ecost = lightP->separationEnergy + Ecoul;
  if (Ecost < EcostMin) EcostMin = Ecost;
  if (exp(-(Ecost-EcostMin)/temp) < threshold && lightP->iA > 1) return 0.;

  S2Start = roundf(lightP->fMInertia/
           (lightP->fMInertia+lightP->fMInertiaOrbit)*fJ);
  if (lightP->odd) S2Start += 0.5;


  //prepare transmission coeff
  lightP->tlArray->prepare(lightP->residue.iZ);

  lightP->width = S2Loop(-1.);
  
  return lightP->width;
}
//***************************************************************
/**
 * Calculates \f$\sum_{\ell=\ell_{min}}^{\ell_{max}} T_{\ell}(\varepsilon)\f$
\param ek is \f$\varepsilon\f$, the kinetic energy in MeV
\param temp is the temperature on the residue in MeV
*/ 
float CNucleus::getSumTl(float ek,float temp)
{
  //function used by hauserFeshbach
  float sumTl = 0.;
  float tlWeightMax = 0.;

  storeSub = new SStoreSub[nSubTl];
  int iSub = 0;
  int iL = lMin;
  float scale = 1.;
  if (lightP->iA == 1 && lightP->iZ == 1)
    {
      //mass at the attractor line
      float Aeal = 2.045*(float)lightP->residue.iZ 
         + 3.57e-3*pow((float)lightP->residue.iZ,2);
      float r_eal = 1.2*pow(Aeal,float(1./3.)) + 1.167;
      float r     = 1.2*pow((float)lightP->residue.iA,(float)(1./3.)) + 1.167;
      scale = r/r_eal;  //rjc
    }
  else if (lightP->iA == 4)
    {
      //mass at the attractor line
      float Aeal = 2.045*(float)lightP->residue.iZ 
         + 3.57e-3*pow((float)lightP->residue.iZ,2);
     float deltaA = (float)lightP->residue.iA - Aeal;
     scale = 1.+deltaA*(-5.412976E-04-9.176213E-02/(float)lightP->residue.iZ
			+ 3.916151E-06*(float)lightP->residue.iZ);
     scale = 1./scale; 
    }   

  for (;;)
    {
      float tl = lightP->tlArray->getTl(iL,ek*scale,temp); 
      //cout << iL << " " << ek << " " << lightP->residue.iZ << " " << tl << endl;
      float maxLplusS = min(iL+lightP->fJ,lPlusSMax);
      float minLplusS = max(fabs(iL-lightP->fJ),lPlusSMin);
      float tlWeight = tl*(maxLplusS - minLplusS + 1.);
      storeSub[iSub].weight = tlWeight;
      if (iSub > 0) storeSub[iSub].weight += storeSub[iSub-1].weight; 
      storeSub[iSub].L = iL;
      iSub++;


      if (iSub > nSubTl)
	{
          cout << " increase nSubTl" << endl;
	  abort();
	}
      sumTl += tlWeight;
      tlWeightMax = max(tlWeight,tlWeightMax);
      if (tlWeight < 0.01*tlWeightMax) break;
      iL++;
      if (iL > lMax) break;
    }


  float xran = ran.Rndm();
  int i = 0;
  for (;;)
    {
      float prob = storeSub[i].weight/sumTl;
      if (prob > xran) break;
      if ( i == iSub-1) break;
      i++;
    }
  EvapL = storeSub[i].L;
  

  delete [] storeSub;
  return sumTl;
}
//**************************************************************************
/**
 * Selects the angle and velocity of an evaporated product
 */
void CNucleus::angleEvap()
{
  //chooses the angles of evaporated particles
  // find orientation of residue spin vector with respect to s0
  // vector (i.e. z axis parallel to s0)
  CAngle residueCM;
 if (EvapLPlusS > 0.0 && fJ > 0.0) 
    residueCM.theta = acos((pow(fJ,2) + pow(EvapLPlusS,2) 
			    - pow(EvapS2,2))/(2.0*EvapLPlusS*fJ));
 else if (EvapLPlusS == 0) residueCM.theta = 0;
 else residueCM.theta = acos(1.-2.*ran.Rndm());

 residueCM.phi = 2.*pi*ran.Rndm();
 // rotate to lab frame 
 daughterHeavy->spin = CAngle::transform(residueCM,spin);
     
 // find direction of emitted fragment
 // first the orientation of the vector l_plus_s with respect to s0
 // (i.e. z axis parallel to s0)

 CAngle lPlusSCM;
 if (EvapS2 > 0.0 && fJ > 0.0)
    lPlusSCM.theta = acos((pow(fJ,2) + pow(EvapS2,2) - pow(EvapLPlusS,2))
			/(2.0*EvapS2*fJ));
 if (EvapS2 == 0) lPlusSCM.theta = 0.; 
 else lPlusSCM.theta = acos(1.0-2.0*ran.Rndm());
 lPlusSCM.phi = pi + residueCM.phi;
 if (lPlusSCM.phi > 2.*pi) lPlusSCM.phi -= 2.*pi;

 //rotate to lab frame
 CAngle lPlusSLab = CAngle::transform(lPlusSCM,spin);

 //angle of L vector - orbital angular momentum
 CAngle LLLab;

if (EvapS1 > 0.0)
  {
    //   find orientation of L vector with respect to l_plus_s vector
    //   (i.e. z axis parallel to l_plus_s)
    CAngle LL;
    if (EvapLPlusS > 0.0)
      LL.theta = acos((pow(EvapLPlusS,2) + pow((float)EvapL,2) - pow(EvapS1,2))
		      /(2.0*EvapLPlusS*(float)EvapL));
    else LL.theta = acos(1.-2.*ran.Rndm()); 
    LL.phi = 2.*pi*ran.Rndm();

    // rotate to lab frame 
    LLLab = CAngle::transform(LL,lPlusSLab);

    // orientation of s1 vector  with respect to l_plus_s
    //          (i.e. z axis parallel to l_plus_s)
    CAngle SS;
    if (EvapLPlusS > 0.0 && EvapL > 0) 
          SS.theta = acos((pow(EvapLPlusS,2) + pow((float)EvapL,2) - 
               pow((float)EvapS1,2))/(2.0*EvapLPlusS*(float)EvapL));
    else if (EvapL == 0) SS.theta = 0.;   
    else SS.theta = acos(1.-2.*ran.Rndm());
    SS.phi =  LL.phi+pi;
    if (SS.phi > 2.*pi) SS.phi -= 2.*pi;
    //   rotate to lab frame
    daughterLight->spin = CAngle::transform(SS,lPlusSLab);
    }
 else // spin zero particle
   {
     LLLab.theta = lPlusSLab.theta;
     LLLab.phi   = lPlusSLab.phi;
     CAngle S1Lab(0.,0.);
     daughterLight->spin = S1Lab;
   }

// emission vector is perpendicular to l vector for classical mechanics
  //   a good approximation can be obtained if theta
  //   is chosen from the legendre distribution sqaured  for m=l
 CAngle emitL;
 emitL.theta = angleDist.getTheta(EvapL);
 emitL.phi = 2.*pi*ran.Rndm();
 
 // rotate into the lab frame 
 CAngle emitLab= CAngle::transform(emitL,spin);

 if (EvapEk < 0.0) 
   { 

     cout << "EvapEk < 0" << " " << EvapEk << " " << iZ << " " << iA 
	  << " " << fEx << " " << fJ << " " << daughterLight->iZ << " " 
	  << daughterLight->iA << " " << daughterLight->fEx << " " <<
       daughterLight->fJ << " " << daughterHeavy->iZ << " " << 
       daughterHeavy->iA << " " << daughterHeavy->fEx << " " <<
       daughterHeavy->fJ << endl;
     CNucleus * parent;
     parent = getParent();
     cout << parent->iZ << " " << parent->iA << " " << parent->fEx << " " << parent->fJ << endl;
     parent = parent->getParent();
     cout << parent->iZ << " " << parent->iA << " " << parent->fEx << " " << parent->fJ << endl;

    abort();
   }

 // get velocities
  float vrel = sqrt(2.0*EvapEk*(float)(EvapA1+EvapA2)/
		    (float)(EvapA1*EvapA2))*0.9794;
  float v1 = (float)EvapA2/(float)(EvapA2+EvapA1)*vrel;

 daughterLight->velocity[0] = v1*sin(emitLab.theta)*cos(emitLab.phi);
 daughterLight->velocity[1] = v1*sin(emitLab.theta)*sin(emitLab.phi);
 daughterLight->velocity[2] = v1*cos(emitLab.theta);
 float v2 = vrel - v1;
 emitLab.theta = pi - emitLab.theta;
 emitLab.phi += pi;
 if (emitLab.phi > 2.*pi) emitLab.phi -= 2.*pi;
 daughterHeavy->velocity[0] = v2*sin(emitLab.theta)*cos(emitLab.phi);
 daughterHeavy->velocity[1] = v2*sin(emitLab.theta)*sin(emitLab.phi);
 daughterHeavy->velocity[2] = v2*cos(emitLab.theta);

 for (int i=0;i<3;i++)
   {
     daughterLight->velocity[i] += velocity[i];
     daughterHeavy->velocity[i] += velocity[i];
   }


}
//*********************************************************
/**
 * selects the velocity of evaporated particles, the emission angles 
 * is assumed isotropic
 */
void CNucleus::angleIsotropic()
{
  //chooses angles and velocities of evaporated particles from isotropic 
  //distribution - used with Weiskopf formalism 
  float theta = acos(1.-2.*ran.Rndm());
  float phi = 2.*pi*ran.Rndm();


  float vrel = sqrt(2.0*EvapEk*(float)(EvapA1+EvapA2)/
		    (float)(EvapA1*EvapA2))*0.9794;
  float v1 = (float)EvapA2/(float)(EvapA2+EvapA1)*vrel;

 daughterLight->velocity[0] = v1*sin(theta)*cos(phi);
 daughterLight->velocity[1] = v1*sin(theta)*sin(phi);
 daughterLight->velocity[2] = v1*cos(theta);
 float v2 = vrel - v1;
 theta = pi - theta;
 phi += pi;
 if (phi > 2.*pi) phi -= 2.*pi;
 daughterHeavy->velocity[0] = v2*sin(theta)*cos(phi);
 daughterHeavy->velocity[1] = v2*sin(theta)*sin(phi);
 daughterHeavy->velocity[2] = v2*cos(theta);

 for (int i=0;i<3;i++)
   {
     daughterLight->velocity[i] += velocity[i];
     daughterHeavy->velocity[i] += velocity[i];
   }


}
//*********************************************************
/**
 * used to check momentum conservation
 *
 * prints out the center-of-mass velocity of all decay products
 */
void CNucleus::vCMofAllProducts()
{
  float momTot[3] = {0.,0.,0.};
  for (int i=0;i<iStable;i++)
    {
      for (int j=0;j<3;j++)
        momTot[j] += stableProducts[i]->velocity[j]*(float)stableProducts[i]->iA;
    }
  float vcm[3];
  for (int j=0;j<3;j++) vcm[j] = momTot[j]/(float)iA;

  cout << "Vcm= " << vcm[0] << " " << vcm[1] << " " << vcm[2] << endl;
}
//************************************************************************
/**
 * Returns the total gamma-ray decay width in MeV
 *
 *  Contributions from E1's and E2's only
 */
float CNucleus::gammaWidth()
{

  float widthTot = 0;
  float width[3];
  float Ex[3];
  float S0[3];
  for (int iMode =1;iMode<=2;iMode++)
    {
      //if (iMode == 1) width[iMode]=gammaWidthE1GDR();
      //else width[iMode] = gammaWidthMultipole(iMode); 
      width[iMode] = gammaWidthMultipole(iMode); 
      Ex[iMode] = GammaEx;
      S0[iMode] = GammaJ;
      widthTot += width[iMode];
    }
  if (widthTot <= 0.) return 0.;
  int iMode;
  if (ran.Rndm() < width[1]/widthTot)iMode = 1 ;
  else iMode = 2;
  GammaEx = Ex[iMode];
  GammaJ = S0[iMode];
  GammaL = iMode;
  return widthTot;
}
//**************************************************************
/**
 * Returns the gamma-ray decay width in MeV for a specified multipole.
 * The width is from Blatt and Weisskopf, "Theoretical Nuclear Physics"
 * (Wiley, New York, 1958) Page=649 scaled by the factors gammaInhibition[]
 * Values of the latter are taken from Phys. Rev. C39, 516 (1989).
\param iMode 1 is E1, 2=E2
*/

float CNucleus::gammaWidthMultipole(int iMode)
{

  //iMode =1 is E1 and Imode=2 is E2
  //width given in units of MeV
  storeEvap = new SStoreEvap [nGamma];

  float deltaE = 1;
  if (fEx-Erot <= 15.0) deltaE = 0.2;

  float  exponent = (float)(2.*iMode)/3.;
  float constant = pow((float)(iA),exponent)*wue[iMode]*pi*2.0;
  constant *= gammaInhibition[iMode];
  float S0 = fJ;
  float S0Min = fabs(S0-(float)iMode);
  float S0Max = S0 + (float)iMode;
  S0 = S0Min;
  float sumGammaMax = 0.0;
  float width = 0.;
  int jj = 0;
  for (;;) //loop over S0 nucleus spin
    {
      float Erotation = yrast.getYrast(iZ,iA,S0);
      float U2Max = fEx - Erotation;
      if (U2Max < 0.01) break;

      float e = 0.1;
      float sumGamma = 0.0;
      float gammaMax = 0.0;
      for (;;) //loop over excitation energy
        {
          float U2 = U2Max - e;
          if (U2 < 0.01) break; 
          float logLevelDensityDaughter = 
	    levelDensity.getLogLevelDensitySpherical
	    (iA,U2,fPairing,fShell,S0,fMInertia);
          if (logLevelDensityDaughter == 0) break;

            float gamma = pow(e,2*iMode+1)*constant*deltaE*
            exp(logLevelDensityDaughter-logLevelDensity)/2./pi;


          storeEvap[jj].gamma = gamma;
          storeEvap[jj].spin = (int)S0;
          storeEvap[jj].energy = e;

          gammaMax = max(gamma,gammaMax);
          sumGamma += gamma;
          if (jj > 0) storeEvap[jj].gamma += storeEvap[jj-1].gamma;
          jj++;
          if (jj > nGamma) 
             {
               cout << " increase nGamma" << endl;
               abort();
             }
          if (gamma < 0.01*gammaMax) break;
          e+=deltaE;
        }

    sumGammaMax = max(sumGamma,sumGammaMax);
    if (sumGamma < 0.01*sumGammaMax) break;
    width += sumGamma;
    S0++;
    if (S0 > S0Max) break;  
   }

  if (jj == 0  || width <= 0.) 
    {
      delete [] storeEvap;
      return 0.;
    }
  int i = 0;
  float xran = ran.Rndm();
  for (;;)
    {
      float prob = storeEvap[i].gamma/width;
      if (xran < prob) break;
      i++;
      if (i == jj) break;
    }

  //choose gamma energy and excitation energy of residual
  float e; // gamma ray energy
  GammaEx = -1.;
  int iTry = 0;
  for (;;)
     {
     e = storeEvap[i].energy  + deltaE*(0.5*ran.Rndm());
     GammaEx = fEx - e;
     if (e > 0. && GammaEx > 0.) break;
     iTry++;
     if (iTry > 10)
       {
	 e = storeEvap[i].energy;
         GammaEx = fEx - e;
         break;
       }
     }
  GammaJ = storeEvap[i].spin;
  if (iA%2 == 1) GammaJ += 0.5;
  delete [] storeEvap;
  return width;
}
//*****************************************************************
/**
 * Returns the gamma-ray decay width in MeV for statistical E1
 * taking into acount GDR strength function
 * see see D.R. Chakrabarty et al. Phys Rev. C36 (1987) 1886
*/
float CNucleus::gammaWidthE1GDR()
{
  //gives the statistical E1 gamma width  in units of MeV
  storeEvap = new SStoreEvap [nGamma];

  float deltaE = 1;
  if (fEx-Erot <= 15.0) deltaE = 0.2;

  //giant dipole strength function
  //see D.R. Chakrabarty et al. Phys Rev. C36 (1987) 1886
  //for details
  float EGDR = 77.0/pow((float)iA,(float)(1./3.));
  float gammaGDR = 4.8+0.0026*pow(fEx-Erot,(float)1.6);
  float Fe1 = 2.09e-5*(float)iZ*(float)(iA-iZ)/(float)(iA)*gammaGDR;

  float S0 = fJ;
  float S0Min = fabs(S0-1.);
  float S0Max = S0 + 1.;
  S0 = S0Min;
  float sumGammaMax = 0.0;
  float width = 0.;
  int jj = 0;
  for (;;) //loop over S0 nucleus spin
    {
      float Erotation = yrast.getYrast(iZ,iA,S0);
      float U2Max = fEx - Erotation;
      if (U2Max < 0.01) break;

      float e = 0.1;
      float sumGamma = 0.0;
      float gammaMax = 0.0;
      for (;;) //loop over excitation energy
        {
          float U2 = U2Max - e;
          if (U2 < 0.01) break; 
          float logLevelDensityDaughter = 
	    levelDensity.getLogLevelDensitySpherical
	    (iA,U2,fPairing,fShell,S0,fMInertia);
          if (logLevelDensityDaughter == 0) break;

 
        float gamma = Fe1*pow(e,4)/
            (pow(pow(e,2)-pow(EGDR,2),2)+pow(gammaGDR*e,2))
            *deltaE*exp(logLevelDensityDaughter-logLevelDensity)/2./pi;


          storeEvap[jj].gamma = gamma;
          storeEvap[jj].spin = (int)S0;
          storeEvap[jj].energy = e;

          gammaMax = max(gamma,gammaMax);
          sumGamma += gamma;
          if (jj > 0) storeEvap[jj].gamma += storeEvap[jj-1].gamma;
          jj++;
          if (jj > nGamma) 
             {
               cout << " increase nGamma" << endl;
               abort();
             }
          if (gamma < 0.01*gammaMax) break;
          e+=deltaE;
        }

    sumGammaMax = max(sumGamma,sumGammaMax);
    if (sumGamma < 0.01*sumGammaMax) break;
    width += sumGamma;
    S0++;
    if (S0 > S0Max) break;  
   }

  if (jj == 0  || width <= 0.) 
    {
      delete [] storeEvap;
      return 0.;
    }

  int i = 0;
  float xran = ran.Rndm();
  for (;;)
    {
      float prob = storeEvap[i].gamma/width;
      if (xran < prob) break;
      i++;
      if (i == jj) break;
    }
  float e = storeEvap[i].energy;
  e += deltaE*(0.5*ran.Rndm());
  if (e < 0.) e=0.;
  GammaEx = fEx - e;
  GammaJ = storeEvap[i].spin;
  if (iA%2 == 1) GammaJ += 0.5;
  delete [] storeEvap;
  return width;
}
//*****************************************************************
/**
 * Selects the angle and sin axis of a daughter product after gamma
 * emission
 */
void CNucleus::angleGamma()
{
  for (int i=0;i<3;i++) daughterHeavy->velocity[i] = velocity[i];

  // find orientation of residue spin vector with respect to s0
  // vector (i.e. z axis parallel to s0)
  CAngle residueCM;
 if (GammaL > 0 && fJ > 0.0) 
   residueCM.theta = acos((pow(fJ,2) + pow((float)GammaL,2) 
			   - pow(GammaJ,2))/(2.0*(float)GammaL*fJ));
 else if (GammaL == 0) residueCM.theta = 0;
 else residueCM.theta = acos(1.-2.*ran.Rndm());

 residueCM.phi = 2.*pi*ran.Rndm();
 // rotate to lab frame 
 daughterHeavy->spin = CAngle::transform(residueCM,spin);
}
//*******************************************************************
/**
 * Return the theta angle of the fragments in radians
 */

float CNucleus::getTheta()
{
  float V = 0.;
  for (int i=0;i<3;i++) V += pow(velocity[i],2);
  return acos(velocity[2]/sqrt(V));
}
//*******************************************************************

/**
 * Returns the theta and phi angle of the fragments in radians
 */

CAngle CNucleus::getAngle()
{
  
  float V = 0.;
  for (int i=0;i<3;i++) V += pow(velocity[i],2);

  float theta;
  float phi;

  if (V > 0.)
    {
      theta = acos(velocity[2]/sqrt(V));
      phi = atan2(velocity[1],velocity[0]);
    }
  else
    {
      theta = 0.0;
      phi = 0.0;
    }

  return CAngle(theta,phi);
}
//*******************************************************************

/**
 * Return the theta angle of the fragments in degrees
 */
float CNucleus::getThetaDegrees()
{
  //get the theta angle of a fragments in degrees
  return getTheta()*180./pi;
}
//************************************************
  /**
   * returns the theta and phi angles in degrees
   */
CAngle CNucleus::getAngleDegrees()
{
  CAngle angle = getAngle();
  angle.theta *= 180./pi;
  angle.phi *= 180./pi;
  return angle;
}
//*****************************************
/**
 * Calculated the evaporation decay width of a channel in MeV using
 * the Weisskopf formulism
 \param saddle bool true=saddleToScission decay false=normal evaporation
*/
float CNucleus::weiskopf( bool saddle)
{
  //calculates the evaporation decay width with the Weiskopf formalism
  lightP->width = 0.;
  lightP->iStore = 0;
  if (exp(-lightP->fEx/temp) < 0.01) return 0.;

  lightP->residue.iZ = iZ - lightP->iZ;
  if (lightP->iZ >= lightP->residue.iZ) return 0.;
  lightP->residue.iA = iA - lightP->iA;
  if (lightP->iA >= lightP->residue.iA) return 0.;
  if (lightP->residue.iZ >= lightP->residue.iA) return 0.;
  //products must lie on the chart of nuclides used by gemini
  if (lightP->residue.iA < mass->chart.getAmin(lightP->residue.iZ)) return 0.;
  if (lightP->residue.iA > mass->chart.getAmax(lightP->residue.iZ)) return 0.;
  rResidue = r0*pow((float)lightP->residue.iA,(float)(1./3.));
  lightP->fMInertia = 0.4*(float)lightP->residue.iA*pow(rResidue,2);
  lightP->fMInertiaOrbit = (float)(lightP->residue.iA*lightP->iA)/
    (float)(lightP->residue.iA+lightP->iA)*pow(rResidue+lightP->rLight+2.,2);

  if (saddle)
    {
      if (fissionMassScission)scission.init(lightP->residue.iZ,lightP->residue.iA,fJ,2);
      else 
	{
         float Z1 = (float)fissionZ/(float)fissioningZ*(float)lightP->residue.iZ;
         float A1 = (float)fissionA/(float)fissioningA*(float)lightP->residue.iA;
         scission.init(lightP->residue.iZ,lightP->residue.iA,fJ,2,Z1,A1);
	}
    }

  //odd A nuclei have 1/2 integer spin
  lightP->odd =0;
  if (lightP->residue.iA%2 == 1) lightP->odd = 1; 

  if (lightP->residue.iZ > Zshell)
  lightP->fPair = mass->getPairing(lightP->residue.iZ,lightP->residue.iA);   
  else lightP->fPair = 0.;

  lightP->residue.fExpMass = mass->getCalMass(lightP->residue.iZ,lightP->residue.iA); 
  lightP->fShell = 0.;
  if (lightP->residue.iZ > Zshell) lightP->fShell = mass->getShellCorrection(lightP->residue.iZ,lightP->residue.iA);   

  lightP->separationEnergy = lightP->residue.fExpMass + lightP->fExpMass - fExpMass;

  float Edef = 0.;
  if (saddle)
    {
       Edef = scission.getScissionEnergy();
       Edef -= lightP->residue.fExpMass;
    }
  else Edef = yrast.getYrast(iZ,iA,fJ);




  float maxGamma = 0.;
  float constant = (2.*lightP->fJ+1.)*de/(2.*pi);
  if (!saddle) constant *= (2*fJ+1.);

  //prepare transmission coeff


  float Ecoul;
  if (Isig)
    {
     lightP->sigBarDist->prepare((float)lightP->residue.iZ,(float)lightP->residue.iA);
     Ecoul = lightP->sigBarDist->getBarrier();
    }
  else 
    {
     lightP->sigBarDist->prepare((float)lightP->residue.iZ,(float)lightP->residue.iA);
     Ecoul = lightP->sigBarDist->getBarrier();
     lightP->tlArray->prepare(lightP->residue.iZ);
    }
  float fEk = Ecoul;
  bool up = 1;
  for (;;)
    {
      float U = fEx - lightP->separationEnergy - fEk - Edef;
      if (U <= 0.) break;

      float logLevelDensity2 = 0.;
      if (saddle) logLevelDensity2 =
	levelDensity.getLogLevelDensityScission(iA,U);
      else logLevelDensity2 = levelDensity.getLogLevelDensitySpherical
	(lightP->residue.iA,U,lightP->fPair,lightP->fShell,-fJ,lightP->fMInertia);

      if (logLevelDensity2 == 0.) break;
      float temp = levelDensity.getTemp();
      float sigmaInverse;
      if (Isig) sigmaInverse = lightP->sigBarDist->getInverseXsec(fEk,temp);
      else sigmaInverse = 
                lightP->tlArray->getInverseXsec(fEk,temp);
      float gamma = sigmaInverse*exp(logLevelDensity2-logLevelDensity)
                    *constant;
      lightP->width += gamma;
      lightP->storeEvap[lightP->iStore].gamma = gamma;
      if (lightP->iStore > 0) lightP->storeEvap[lightP->iStore].gamma += lightP->storeEvap[lightP->iStore-1].gamma;
      lightP->storeEvap[lightP->iStore].energy = fEk;
      lightP->iStore++;
      if (lightP->iStore == lightP->nStore)
	{
	  cout << "increase nStore in nucleus.h "<< endl;
          abort();
	}
      maxGamma = max(maxGamma,gamma);

      if (gamma < maxGamma*.01) 
	{
	  if (up)
	    {
	      up = 0;
              fEk = Ecoul - de;
	      if (fEk <= 0.) break;
	      else continue; 
	    }
	  else break;
          
	}
      if (up) fEk += de;
      else 
	{
         fEk -= de;
	 if (fEk <- 0.) break;
	}
    }





  //  float del = 0.;
  //if (lightP->residue.iA%2 == 1)del =  0.5;



  return lightP->width;
}
//*****************************************
/**
 * Calculates the decay width for evaporation at the scission point using
 * the Weiskopf formalism 
 *
 * The width is given in units of MeV
 */
float CNucleus::evaporationWidthSS()
{
  float width = 0.;

  EcostMin = 1000.;
  bool isaddle = 1;

  for (int i=0;i<evap.nLight;i++) 
    {
      lightP = evap.lightP[i];
      evap.prob[i] = weiskopf(isaddle);
      width += evap.prob[i];
      if (i > 0) evap.prob[i] += evap.prob[i-1];
    }
  if (width <= 0.) return 0.;


  //if decay possible choose channel
  float xran = ran.Rndm();
  int i = 0;
  for (;;)
    {
      float prob = evap.prob[i]/width;
      if (prob > xran) break;
      if ( i == evap.nLight-1) break;
      i++;
    }
  lightP = evap.lightP[i];

  EvapZ1 = lightP->iZ;
  EvapA1 = lightP->iA;
  EvapZ2 = iZ - EvapZ1;
  EvapA2 = iA - EvapA1;
  EvapMode = i;
  EvapEx1 = lightP->fEx;
  EvapS1 = lightP->fJ;

  getSpin(bool(1));


  return width;
}
//*****************************************
/**
 * Initializes the compound nucleus excitation and spin.
 \param fEx0 is the compound nucleus excitation energy in MeV
 \param fJ0 is the compound nucleus spin in hbar
*/
void CNucleus::setCompoundNucleus(float fEx0, float fJ0)
{
  //initialize as a compound nucleus ready for decay
  excite(fEx0,fJ0);
  origin = 0;
  timeSinceStart = 0.;
  sumGammaEnergy = 0.;
  iPoint = -1;
  runningWeight = 1.;
  iWeight = 0;
  bResidue  = true;
  bSymmetricFission = false;
  bAsymmetricFission = false;
}
//*****************************************
/**
 * Initializes the compound nucleus excitation and spin.
 \param dEx0 is the compound nucleus excitation energy in MeV
 \param fJ0 is the compound nucleus spin in hbar
*/
void CNucleus::setCompoundNucleus(double dEx0, float fJ0)
{
  setCompoundNucleus((float)dEx0,fJ0);
}
//********************************************
/**
 * Initializes the compound nucleus excitation and spin.
 \param fEx0 is the compound nucleus excitation energy in MeV
 \param dJ0 is the compound nucleus spin in hbar
*/
void CNucleus::setCompoundNucleus(float fEx0, double dJ0)
{
  setCompoundNucleus(fEx0,(float)dJ0);
}
//**********************************************
/**
 * Initializes the compound nucleus excitation and spin.
 \param dEx0 is the compound nucleus excitation energy in MeV
 \param dJ0 is the compound nucleus spin in hbar
*/
void CNucleus::setCompoundNucleus(double dEx0, double dJ0)
{
  setCompoundNucleus((float)dEx0,(float)dJ0);
}
//*************************************************
/**
 * Returns a pointer to a stable decay product.
 *
 * if no input given, the first or next product is pointed to.
 /param i is the index of the product (0-getNumberOfProducts-1)
 */ 
CNucleus * CNucleus::getProducts(int i/*=-1*/)
{
  
  //returns pointer to array of stable products
  if (i >= 0) iPoint = i;
  else iPoint++;

  return stableProducts[iPoint];
}
//*************************************************
/**
 * Returns a pointer to the parent nucleus,
 *  i.e. the nucleus which emitted this product.
 * This is useful to see if there was secondary decay.
 * If the NULL pointer is returned, then the 
 * initial compound nucleus did not decay presumably
 * because there was not enough excitation energy.
 * Obveriously in this case the final product is the same 
 * as the initial nucleus and there is no parent.
 * 
 */ 
 CNucleus * CNucleus::getParent()
 {
  
  return parent;
 }
//*****************************************
/**
 *Returns a pointer to the heavy daughter nucleus
 *If NULL is returned, then this fragment was stable
 */
CNucleus * CNucleus::getHeavyDaughter()
{
  return daughterHeavy;
}
//*****************************************
/**
 * Returns a pointer to the light daughter nucleus
 * if NULL is returned, then this fragment was stable
 * or a fission event stated, in which case the 
 * heavy daughter pointer will not be NULL
 */
CNucleus * CNucleus::getLightDaughter()
{
  return daughterLight;
} 
//**************************************************
/**
 * Returns the number of stable decay products produced in the 
 *the statistical decay
 */

int CNucleus::getNumberOfProducts()
{
  return iStable;
}
//**************************************************************
/**
 * Causes the nucleus to undergo statistical decay
 */
void CNucleus::decay()
{
  int iProductsOld = iProducts;
  int iStableOld = iStable;
  int tries = 0;
  for (;;)
    {
     recursiveDecay();
     if (abortEvent == 0) break;
     //cout << "aborted event, try again" << endl;

     // valid event could not be generated
     // set arrays points back to NULL and try again

     abortEvent = 0;
     if (iProducts > iProductsOld)
       {
         for (int i=iProducts-1;i>=iProductsOld;i--)
            {
              allProducts[i]->daughterLight = NULL;
              allProducts[i]->daughterHeavy = NULL;
              delete allProducts[i];
              allProducts[i] = NULL;
            }
       }

     if (iStable > iStableOld)
       {
         for (int i=iStable-1;i>=iStableOld;i--)
            {
              stableProducts[i] = NULL;
            }
       }

     iProducts = iProductsOld;
     iStable = iStableOld;
     daughterLight = NULL;
     daughterHeavy = NULL;
     tries++;
     if (tries > 3) 
       {
	 cout << "given up trying to decay this nucleus" << endl;
         cout << " Z = " << iZ << " A= " << iA <<  " Ex = " << fEx <<
              " J= " << fJ << endl;
         abortEvent = 1;
         return;
       }
    }



  multPostLight = 0;
  multPostHeavy = 0;
  multSaddleToScission = 0;
  multPreSaddle = 0;
  for (int i=0;i<iStable;i++)
    {
      if (stableProducts[i]->iZ == 0 && stableProducts[i]->iA == 1)
	{
	  if (stableProducts[i]->origin == 2) multPostLight++;
	  if (stableProducts[i]->origin == 3) multPostHeavy++;
	  if (stableProducts[i]->origin == 1) multSaddleToScission++;
	  if (stableProducts[i]->origin == 0) multPreSaddle++;
	}
    }
}
//************************************************************
bool CNucleus::isSymmetricFission()
/**
 * Returns a true value if a symmetric fission occurred
 * Only use this function for the compound nucleus object, otherwise
 * the output is garbage.
 */
{
  return bSymmetricFission;
}
//**********************************************
bool CNucleus::isAsymmetricFission()
/**
 * Returns a true value if an asymmetyric fission occurred
 * Only use this function for the compound nucleus object, otherwise
 * the output is garbage.
 */
{
  return bAsymmetricFission;
}
//**************************************************************
/**
 * Returns a true value if the event doesn't fission and has an 
 * evaporation residue.
 * Only use this function for the compound nucleus object, otherwise
 * the output is garbage.
 */
bool CNucleus::isResidue()
{
  return bResidue;
}

//*****************************************************************
  /**
   * returns a true value if the fragment does not undergo 
   * statistical decay, for example a particular excited state of 
   * a nucleus. These nuclei are produced in evaporation processes.
   */
bool CNucleus::isNotStatistical()
{
  return notStatistical;
}
//****************************************************************
  /**
   *  returns true if the nucleus is undergoing a saddle to scission
   *  transition.  All symmetric fission events, pass thought this 
   *  stage and some emit light particles during this stage.
   */
bool CNucleus::isSaddleToScission()
{
  return saddleToSciss;
}
//****************************************************************
int CNucleus::getMultPost()
/**
 * Returns the number of neutrons emitted from both fission fragments
 */
{
  return multPostLight + multPostHeavy;
}
//****************************************************************
int CNucleus::getMultPre()
/**
 * Returns the number of neutrons emitted before the scission point
 */
{
  return multPreSaddle + multSaddleToScission;
}
//***************************************************************
int CNucleus::getMultPostLight()
/**
 * returns the multiplicity of neutrons emitted from the lighter
 * fission fragment.
 */
{
  return multPostLight;
}
//*****************************************************************
int CNucleus::getMultPostHeavy()
/**
 * returns the multiplicity of neutrons emitted from the heavier 
 * fission fragment.
 */
{
  return multPostHeavy;
}
//****************************************************************
/**
 * Returns the multiplicity of neutrons emitted before the saddle-point
 */
int CNucleus::getMultPreSaddle()
{
  return multPreSaddle;
}
//************************************************************
/**
 * Returns the multiplicity of neutrons emitted between saddle
 * and scission
 */
int CNucleus::getMultSaddleToScission()
{
  return multSaddleToScission;
}
//************************************************************
/**
 * used to check for conservation of energy.
 *
 * this is a debugging tool. prints out the various contributions 
 * to teh final energy and the energy difference between inital and final 
 * states. The latter shoukd be zero is energy is conserved. 
 */
void CNucleus::energyConservation()
{
  float sumKE = 0.;
  float sumMass = 0.;
  float sumEx = 0.; 
  for (int i=0;i<iStable;i++)
    {
      sumKE += stableProducts[i]->getKE();
      sumMass += stableProducts[i]->getExcessMass();
      sumEx += stableProducts[i]->fEx;
    }
  float Qvalue = fExpMass - sumMass;
  cout << sumKE << " " << Qvalue << " " << sumEx << " " <<
      sumKE-Qvalue+sumEx+sumGammaEnergy 
       << " " <<  fEx - sumKE + Qvalue - sumEx - sumGammaEnergy<< endl;
}
//**********************************************************
float CNucleus::getKE()
/**
 * Returns the fragments kinetic energy in MeV.
 */
{

  float KE = 0.;
  for (int i=0;i<3;i++) KE += pow(velocity[i]/.9784,2);
  KE *= 0.5*(float)iA;
  return KE;
}
//**********************************************************
/**
 * Returns the magnitude of the fragment's velocity in cm/ns
 */
float CNucleus::getVelocity()
{

  float vel = 0.;
  for (int i=0;i<3;i++) vel += pow(velocity[i],2);
  vel = sqrt(vel);
  return vel;
}
//**********************************************************
/**
 * Returns the magnitude of the fragment's momentum in MeV/c
 */
float CNucleus::getMomentum()
{

  return getVelocity()*(float)iA/30.*931.5016;
}
//********************************************************
/**
 * Determines the angles and velocities of the daughter fragments 
 * when the decay was determined by the transition-state formulism
 */
void CNucleus::split(CAngle symmetryCM)
{
  //angles and velocity of fragments after binary decay (not evaporation)
  float A1 = (float)fissionA;
  float A2 = (float)(iA-fissionA);
  float Ared = A1*A2/(A1+A2); 
  //float Ek = yrast.viola((float)fissionZ,A1,
  //			 (float)(iZ-fissionZ),A2);

  float Ek = scission.ekTot;

  float Vrel = sqrt(2.0*Ek/Ared)*0.9794;
  float V1CM = A2/(A1+A2)*Vrel;
  float V2CM = Vrel-V1CM;

  // in the c-o-m frame the velocity vectors are

  daughterLight->velocity[0] = V1CM*sin(symmetryCM.theta)
                                  *cos(symmetryCM.phi);
  daughterLight->velocity[1] = V1CM*sin(symmetryCM.theta)
                                  *sin(symmetryCM.phi);
  daughterLight->velocity[2] = V1CM*cos(symmetryCM.theta);



  // fragment 2 get emitted 180 degrees away.
  symmetryCM.theta = pi - symmetryCM.theta;
  symmetryCM.phi += pi;
  if (symmetryCM.phi > 2.*pi) symmetryCM.phi -= 2.*pi;


  daughterHeavy->velocity[0] = V2CM*sin(symmetryCM.theta)
                                 *cos(symmetryCM.phi);
  daughterHeavy->velocity[1] = V2CM*sin(symmetryCM.theta)
                                 *sin(symmetryCM.phi);
  daughterHeavy->velocity[2] = V2CM*cos(symmetryCM.theta);

  for (int i=0;i<3;i++)
    {
     daughterLight->velocity[i] += velocity[i];
     daughterHeavy->velocity[i] += velocity[i];
    }


     return;
}
//**********************************************************
/**
 *  can be used to initialize another compound nucleus
\param iZ0 is the proton number
\param iA0 is the mass number
\param fEx0 is the excitation energy in MeV
\param fJ0 is the spin in hbar
 */
void CNucleus::setNewIsotope(int iZ0, int iA0, float fEx0, float fJ0)
{
  init(iZ0,iA0);
  setCompoundNucleus(fEx0,fJ0);
} 
//*******************************************************
/**
 * sets the transient time (fission delay) in zs different from the default 
 * value. This is a static function.
 */
 void CNucleus::setTimeTransient(float time)
{
  timeTransient = time;
}
//*******************************************************
/**
 * returns the transient time (fission delay) in zs different from the default 
 * value. This is a static function.
 */
 float CNucleus::getTimeTransient()
{
  return timeTransient;
}
//**************************************************************
/**
 * returns a pointer to the array containing the velocity vector.
 * units are in cm/ns
 */
float* CNucleus::getVelocityVector()
{
  return velocity;
}
//**************************************************************
/**
 * returns a pointer to the arrays containing the momentum vector.
 * units are in MeV/c
 */
float* CNucleus::getMomentumVector()
{
  for (int i=0;i<3;i++) momentum[i] = velocity[i]*(float)iA/30.*931.5016;
  return momentum;
}
//****************************************************
/**
 * sets the fission scale factor
\param factor is the scale factor
*/
void CNucleus::setFissionScaleFactor(float factor)
{
  fissionScaleFactor = factor;
}
//***********************************************
/**
 * returns the fission scale  factor
*/
float CNucleus::getFissionScaleFactor()
{
  return fissionScaleFactor;
}
//**************************************************
/**
 * set the parameter controlling the width of the barrier distribution
\param width - width is \f$ \sqrt(T)* width \f$
 */
void CNucleus::setBarWidth(float width)
{
  CTlBarDist::setBarWidth(width);
  CSigBarDist::setBarWidth(width);
}
//***************************************************
/**
 * returns the paramter controlling the width of the barrier dist
 */
float CNucleus::getBarWidth()
{
  return CTlBarDist::getBarWidth();
}
//****************************************************
/**
 * set the width of the energy bins for integating the Hauser-Feshbach
 *formulism
\param de0 with of the energy bin
*/
void CNucleus::setDeltaE(float de0)
{
  de = de0;
}
//*******************************************************
/**
 * returns the energy bin width used for integrating the Hauser-Feshbach
 * formulism
 */
float CNucleus::getDeltaE()
{
  return de;
}
//********************************************************
/**
 * sets the threshold used to cut out low probability evaporation
 * channels. All channels will be included if set to zero.
\param threshold0 is the threshold
*/
void CNucleus::setThreshold(float threshold0)
{
  threshold = threshold0;
}
//*******************************************************
/**
 * returns the threshold used to cut out low probability evaporation
 * channels
 */
float CNucleus::getThreshold()
{
  return threshold;
}
//***********************************************************
  /**
   * the symmetric fission barrier is modified by adding this quantity
   \param barAdd0 barrier adjestment in MeV
  */
void CNucleus::setAddToFisBarrier(float barAdd0)
{
  barAdd = barAdd0;
}
//***********************************************************
  /**
   * returns the quantity by which the symmerty fission is added to.

  */
float CNucleus::getAddToFisBarrier()
{
  return barAdd;
}
//************************************************************
  /**
   * prints out the values of the statistical model parameters
   */
void CNucleus::printParameters()
{
  CLevelDensity::printParameters();
  CTlBarDist::printParameters();
  CYrast::printParameters();
  cout << "barAdd = " << barAdd << endl;
  cout << "fissionScaleFactor = " << fissionScaleFactor << endl;

  if (fissionMassScission)cout << " viscosity= " << viscosity_scission << endl;
  else cout << " viscosity= " << viscosity_saddle << endl;

  cout << " fission delay = " << timeTransient <<" e-21 s" << endl;
  if (noIMF) cout << " IMF emissions turned off" << endl;
  else cout << " IMF turned on" << endl;
  if (fissionMassScission) 
    cout << " fission mass dist determined at scission" << endl;
  else cout << " fission mass dist determined at saddle" << endl;
  if (BohrWheeler) cout << " Bohr-Wheeler decay width" << endl;
  else cout << " Lestone Decay width" << endl;
  if (barAdd != 0.) cout << barAdd  << " added to barriers" << endl;

  //  cout << "WignerAdd= " << WignerAdd << " WignerScaler= " << WignerScaled <<  endl;
  cout << "WignerAdd= " << WignerAdd << endl;

  if (iHF == 1) cout << "Hauser Feshback for evaporation" << endl;
  else if (iHF ==0)
    {
     cout << "Weisskopf for evaporation with" << endl;
     if (Isig) cout << " parameterized inverse xsections" << endl;
     else cout << 
     " inverse xsection from sum of parameterized transmission coef" << endl;
    }
  else if (iHF ==2) 
    {
     cout << "switches from HF to Weisskopf with" << endl;
     if (Isig) cout << " parameterized inverse xsections" << endl;
     else cout << 
      " inverse xsection from sum of parameterized transmission coef" << endl;
     }
}
//************************************************************
  /**
   * Gives a correction to either the Bohr-Wheeler or Morreto formalism
   * when the tilting angular momentum bearing mode is considered.
   * See J. Lestone PRC 59 (1999) 1540
   * the saddle-point shape is assumed constant as a function of K,
   * the projection of the spin on the symmetry axis.
   \param Usaddle is saddle point excitation energy in Mev
   \param momInertiaEff is the effective moment of inertia for tilting
   \param iAfAn switch to allow af/an value.
  */

float CNucleus::LestoneCorrection( float Usaddle, float momInertiaEff,
                                                   short iAfAn)
{
 float saddleLD = levelDensity.getLogLevelDensitySpherical
    (iA,Usaddle,(float)0.,(float)0.,fJ,fMInertia,iAfAn);
 float saddleTemp = levelDensity.getTemp();

  float K = fJ - roundf(fJ);

  float tot = 0.;
  for (;;)
    {
      if (K > fJ) break;

      if (K < 0.1)  tot = 1.;
      else
	{


      
          float ErotK = kRotate/2.*pow(K,2)/momInertiaEff;
          float U = Usaddle - ErotK; 
          if (U <= 0.) break;
          float LD = levelDensity.
             getLogLevelDensitySpherical(iA,U,(float)0.,(float)0.,
				       fJ,fMInertia,iAfAn);

          float yield = exp(LD-saddleLD)*2.*levelDensity.getTemp()/saddleTemp;

          if (yield < 1.e-3) break;
          tot += yield;
	}
      K++;
    }

  return tot/(2.*fJ+1.);   
}
//************************************************************************
  /**
   * returns a pointer to the Compound nucleus
   */
CNucleus* CNucleus::getCompoundNucleus()
{
  CNucleus * parent = this;
  for(;;)
    {
      CNucleus * parentNew = parent->getParent();
      if (parentNew == NULL) break;
      parent = parentNew;
    }
  return parent;
}


//***********************************************************

/** 
 * Hauser-Feshbach routine to sum of the spin of the residual nucleus
\param Ekvalue if < 0, then loop over Ek, other for the specificed Ekvalue
 */
float CNucleus::S2Loop(float Ekvalue)
{

  lightP->iStore = 0;
  S2 = S2Start;
  float width0 = S2Width(Ekvalue);
  S2 += 1.;
  float width1 = S2Width(Ekvalue);

  float width = width0 + width1;
  if (width <= 0.) return 0.;
  float sign = 1.;
  if (width1 < width0) sign = -1.;
  float gammaMax = max(width0,width1);

  //increase or decrease S2
  if (sign > 0.) S2++;
  else S2 = S2Start - 1;
  for (;;)
    {
      if (S2 < 0.) break;
      float gamma = S2Width(Ekvalue);
      gammaMax = max(gammaMax,gamma);
      width += gamma;
      if (gamma < 0.01*gammaMax) break;
      S2 += sign;

    } 

  //now go the other direction
  if (sign < 0.) S2 = S2Start + 2;
  else S2 = S2Start - 1;
  sign *= -1.;
  for (;;)
    {
      if (S2 < 0.) break;
      float gamma = S2Width(Ekvalue);
      gammaMax = max(gammaMax,gamma);
      width += gamma;
      if (gamma < 0.01*gammaMax) break;
      S2 += sign;

    } 

  return width;
}

//********************************************************************
  /** 
   * calculates the Hauser-Feshbach decay width for a single S2 values, 
   * but integrated over l and ek
   \param Ekvalue if < 0, then loop over Ek, other for the specificed Ekvalue
   */ 
float CNucleus::S2Width(float Ekvalue)
{

      //   consider particle"s spin in determining
      //   the l (orbital angular momentum limits)
      lPlusSMin = fabs(fJ-S2);   //minimum value of (l+s) vector
      lPlusSMax = fJ + S2;    //maximum value of (l+s) vector

      lMax = (int)round(lPlusSMax+lightP->fJ); //maximum value of l vector
      if (lMax > lMaxQuantum) lMax = lMaxQuantum;
      lMin = (int)round(lPlusSMin-lightP->fJ); // minimum value of l 


      if  (lMin < 0) lMin = (int)round(lightP->fJ-lPlusSMax);
      if  (lMin < 0) lMin = 0;  

  
      EYrast2 = yrast.getYrast(lightP->residue.iZ,lightP->residue.iA,S2);

      UMin = fEx - lightP->separationEnergy - EYrast2;
      if (UMin < 0) return 0;

      if (Ekvalue < 0)return EkLoop();
      else return EkWidth(Ekvalue);
}

//********************************************************************
  /**
   * Calculates the Hauser-Feshbach decay width for a single S2 and ek, 
   * but integrated over l
   */
float CNucleus::EkWidth(float ek)
{
     float U2 = UMin-ek;
     float daughterLD = levelDensity.
      getLogLevelDensitySpherical(lightP->residue.iA,U2,lightP->fPair,lightP->fShell
           ,S2,lightP->fMInertia);
     if (daughterLD  == 0 ) return 0.;
     float temp = levelDensity.getTemp();
     float sumTl = getSumTl(ek,temp);
     if (sumTl == 0) return 0.;

     float gamma = de*sumTl*exp(daughterLD-logLevelDensity)/2./pi;
     if (gamma <= 0.) return 0.;

     lightP->storeEvap[lightP->iStore].gamma = gamma;
     lightP->storeEvap[lightP->iStore].spin = (int)S2;
     lightP->storeEvap[lightP->iStore].energy = ek;
     lightP->storeEvap[lightP->iStore].L = EvapL;
     if (lightP->iStore > 0) lightP->storeEvap[lightP->iStore].gamma += 
           lightP->storeEvap[lightP->iStore-1].gamma;
     lightP->iStore++;
     if (lightP->iStore > lightP->nStore)
       {
         cout << " increase nStore in Ek loop" << endl;
         abort();
       }
     /*
     if (lightP->iA == 1 && lightP->iZ == 1)
       {
         cout << fEx - lightP->separationEnergy - ek << " " << S2 << " " << gamma << 
	   " " << daughterLD << " " << sumTl << " " << ek << endl;
       }
     */
     return gamma;

}
//******************************************
/** 
 *Hauser-Feshbach function to sum over kinetic energy for a given S2
 */
float CNucleus::EkLoop()
{

  //find energy at which tl=0.5 approximately
 float ekMin = (float)((lMin+1)*lMin)*20.454/
   lightP->fMInertiaOrbit + Ecoul;
 if( ekMin > UMin) ekMin = UMin;
 //start just below the barrier
 ekMin *= 0.8;
 ekMin = round(ekMin+de/2.0) - de/2.0;
 if (ekMin < de/2.) ekMin = de/2.;


 float ek = ekMin;
  float width0 = EkWidth(ek);
  ek += de;
  float width1 = EkWidth(ek);

  float width = width0 + width1;
  if (width <= 0.) return 0.;
  float sign = 1.;
  if (width1 < width0) sign = -1.;
  float gammaMax = max(width0,width1);

  //increase or decrease S2
  if (sign > 0.) ek+= de;
  else ek = ekMin - de;
  for (;;)
    {
      if (ek < 0.) break;
      float gamma = EkWidth(ek);
      gammaMax = max(gammaMax,gamma);
      width += gamma;
      if (gamma < 0.01*gammaMax) break;
      ek += sign*de;

    } 

  //now go the other direction
  if (sign < 0.) ek = ekMin + 2.*de;
  else ek = ekMin - de;
  sign *= -1.;
  for (;;)
    {
      if (ek < 0.) break;
      float gamma = EkWidth(ek);
      gammaMax = max(gammaMax,gamma);
      width += gamma;
      if (gamma < 0.01*gammaMax) break;
      ek += sign*de;

    } 

  return width;
}

//****************************************************
/**
 * Turns off imf emissions. The default option is
 * IMF emsiion on.
*/
void CNucleus::setNoIMF()
{
  noIMF = 1;
}
//****************************************************
/**
 * Turns on imf emissions. This is the default option.
*/
void CNucleus::setYesIMF()
{
  noIMF = 0;
}
//****************************************************
/**
 * Calculations use Lestone fission width.
 * the default is the BohrWheeler decay width
*/
void CNucleus::setLestone()
{
  BohrWheeler = 0;
}
//****************************************************
/**
 * Calculations use BohrWheeler fission width
 * This is the default. Alternative is setLestone()
*/
void CNucleus::setBohrWheeler()
{
  BohrWheeler = 1;
}
//********************************************
  /**
   * Set one of the parameter sets
  \param isol solution number
  */
void CNucleus::setSolution(int isol)
{
  if (isol == 0)
    {
      setBohrWheeler();
      setTimeTransient(0.);
      CLevelDensity::setAfAn(1.036);
      CLevelDensity::setLittleA(7.3,.00517,.0345,12.);
      CTlBarDist::setBarWidth(1.);

    }
  else if (isol == 1)
    {
      setLestone();
      setTimeTransient(1.);
      CLevelDensity::setAfAn(1.057);
      CLevelDensity::setLittleA(7.3,.00517,.0345,12.);
      CTlBarDist::setBarWidth(1.);
    }
  else 
    {
      cout << "this solution not defined" << endl;
    } 
  
}
//***********************************************
  /**
   * returns the total energy removed by gamma rays in MeV
   */
float CNucleus::getSumGammaEnergy()
{
  return sumGammaEnergy;
}
//*************************************************
  /**
   * returns the time in zs at which the particle was created after 
   * the formation of the CN
   */
float CNucleus::getTime()
{
  return timeSinceStart;
}
//*****************************************************************
  /**
   * determined the spin of the residue (for Weisskopf)
    \param saddle bool true=saddleToScission decay false=normal evaporation
   */
void CNucleus::getSpin(bool saddle)
{
  if (HF && !saddle)
    {
     //choose the residue spin and evaporated particles energy;
     float xran = 1.5;
     while(floor(xran) > 0.5)xran = ran.Rndm();

     int i = 0;
     for(;;)
       {
         float prob = lightP->storeEvap[i].gamma/lightP->width;
         if (prob > xran) break;
         if ( i == lightP->iStore-1) break;
         i++;
       }
     float Ek,Ex;
     int iTry = 0;
     for (;;)
       {
         Ek = lightP->storeEvap[i].energy + (1.-2.*ran.Rndm())*de/2.;
         Ex = fEx - lightP->separationEnergy - Ek;
         if (Ek >= 0. && Ex >= 0.) break;
         if (iTry ==4) 
	   {
	     Ek = lightP->storeEvap[i].energy;
             Ex = fEx - lightP->separationEnergy - Ek;
             break;
   	   }
         iTry++;
       }
  
     EvapEx2 = Ex;
     EvapS2 = lightP->storeEvap[i].spin;
     if (lightP->odd) EvapS2 += 0.5;
     EvapEk = Ek;
     EvapL = lightP->storeEvap[i].L;
     return;  

    }
   

  //now for HF == 0
  //choose the residue spin and evaporated particles energy;
  float xran = ran.Rndm();
  int i = 0;
  for(;;)
    {
      float prob = lightP->storeEvap[i].gamma/lightP->width;
      if (prob > xran) break;
      if ( i == lightP->iStore-1) break;
      i++;
    }

  float Ek,Ex;
  int iTry = 0;
  for (;;)
    {
      Ek = lightP->storeEvap[i].energy + (1.-2.*ran.Rndm())*de/2.;
      Ex = fEx - lightP->separationEnergy - Ek;
      if (Ek >= 0. && Ex >= 0.) break;
      if (iTry ==4) 
	{
	  Ek = lightP->storeEvap[i].energy;
          Ex = fEx - lightP->separationEnergy - Ek;
          break;
	}
      iTry++;
    }

  EvapEx2 = Ex;
  EvapEk = Ek;




  if (saddle)  
    {
     EvapS2 = fJ - lightP->fJ;
     EvapL = 0;
    }
  else 
    {

     lightP->iStore = 0;
     S2Start = roundf(lightP->fMInertia/
              (lightP->fMInertia+lightP->fMInertiaOrbit)*fJ);

     if (lightP->odd) S2Start += 0.5;

     //prepare transmission coeff
     lightP->tlArray->prepare(lightP->residue.iZ);
     float width = S2Loop(Ek);


     //choose the residue spin and evaporated particles energy;
     float xran = 1.5;
     while(floor(xran) > 0.5)xran = ran.Rndm();

     int i = 0;
     for(;;)
       {
         float prob = lightP->storeEvap[i].gamma/width;
         if (prob > xran) break;
         if ( i == lightP->iStore-1) break;
         i++;
       }

     EvapS2 = lightP->storeEvap[i].spin;
     if (lightP->odd) EvapS2 += 0.5;
     EvapL = lightP->storeEvap[i].L;

     
     /*
      float sigmaInverse;
      if (Isig) sigmaInverse = lightP->sigBarDist->getInverseXsec(Ek,0.);
      else sigmaInverse = 
                lightP->tlArray->getInverseXsec(Ek,0.);

      float lmax = sqrt(sigmaInverse);
      int lc = (int)((lmax+.5)*sqrt(ran.Rndm()));
      float jmin = fabs(fJ - (float)lc);
      jmin = jmin - lightP->fJ;
      for (;;)
	{
          if (jmin > 0.) break;
          jmin++;
	}
      float jmax = fJ + (float)lc + lightP->fJ;
      for (;;)
	{
	  if (jmax < 90.) break;
	  jmax--;
	}

      int NN = (int) (jmax - jmin) + 1;
      float prob[NN];
      float U = Ex -  Edef;
      int ii = 0;
      for (float jj=jmin;jj<jmax+1.;jj+=1.)
	{
         prob[ii] =  levelDensity.getLogLevelDensitySpherical
	   (lightP->residue.iA,U,fPairing,fShell,jj,fMInertia2)/(2.*jj+1);
         if (ii > 0) prob[ii] += prob[ii-1];
	 ii++;
	}
      if (ii > NN) 
	{
	  cout << "increase NN" << endl;
         abort();     
	}
      float x = ran.Rndm();
      int i = 0;
      for (;;)
	{
          if (prob[i]/prob[ii-1] > x) break;
          i++;
	}
      
      storeChan[iChannel].S2 = jmin + (float) i;
      storeChan[iChannel].Ek = Ek;
      storeChan[iChannel].L = 0;
     */
      /*
      //find mean and sigma of angular momentum distribution of residue

      float Mreduced = (float)(lightP->iA*lightP->residue.iA)/
                       (float)(lightP->iA+lightP->residue.iA);
      float dist = 1.4*pow((float)lightP->residue.iA,(float)(1./3.)) + 2.;
      float fMInertiaOrbit = Mreduced*pow(dist,2);

      //float fS2 = (fJ - lightP->fJ)*fMInertia2/(fMInertia2+fMInertiaOrbit);
      float fS2 = fJ*fMInertia2/(fMInertia2+fMInertiaOrbit);

      float U = Ex -  Edef;
      levelDensity.getLogLevelDensitySpherical
	(lightP->residue.iA,U,fPair2,fShell2,0.,fMInertia2);
      float temp = levelDensity.getTemp();
      float sigma =  0.1546*sqrt(temp*fMInertiaOrbit*fMInertia2/(
                           fMInertiaOrbit+fMInertia2));

      //now pick the final angular momentum
      float value;
      for(;;)
	{      
         value = ran.Gaus(fS2,sigma);
	 if (value >= 0.) break;
	}
      if (del == 0)storeChan[iChannel].S2 = floor(value+.5);
      else storeChan[iChannel].S2 = floor(value)+.5;
      if (storeChan[iChannel].S2 < 0.) storeChan[iChannel].S2 = del;
      storeChan[iChannel].Ek = Ek;
      storeChan[iChannel].L = 0;
      */

    }
}
//********************************************
  /**
   * Sets the mode use for evaporation
   * 1 = Hauser-Feshach formalism as in the original GEMINI
   * 0 = widths calculated from Weisskopf, then kinetic energy
   *     of particle also from Weisskopf, but spin and orbital
   *     angular momentum from Hauser-Feshbach
   * 2 = Switches bewteen options 0 and 1 dependeing of the ratio of 
   *     rotational to thermal energy. (default)
   *
  \param iHF0 =0,1,2
  */
void CNucleus::setEvapMode(int iHF0/*=2*/)
{
  iHF = iHF0;
}

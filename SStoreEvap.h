#ifndef sstoreEvap_
#define sstoreEvap_

/**
 *!\brief storage
 *
 * this structure store information on an evaporation channel (n,p,..)
 * - the final quantities are Monte-Carlo selected from the many
 * subchannels (SStoreChan)
 */ 





struct SStoreEvap
{
  float gamma; //!< decay width of channel in MeV
  float energy; //!< energy of evaporated particle in MeV
  short unsigned  spin; //!< spin of daughter
  short unsigned L; //!< orbitla angular momentum of evaporated particle
};
#endif

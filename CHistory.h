#ifndef chistory_h
#define chistory_h

#include "CNucleus.h"
#include <map>
#include <limits>
#include <iostream>

// Step types for history
enum HistoryStepType {
	Evaporation=1,
    EvaporationResidue,
	Multifragmentation,
	AsymmetricFissionLight,
	AsymmetricFissionHeavy,
	SymmetricFissionLight,
	SymmetricFissionHeavy,
	SaddleToScission,
	NonStatistical
};



/**
 *!\brief figures out the history of a fragment
 *
 * goes through the event tree structure and determined 
 * the origin of the particle; evaporation products, symmetric,
 * or asymmetric fission products, nonstatistical decay, or 
 * saddle-to-scission emission
 */

class CHistory {

	typedef std::map<CNucleus *, int32_t> HistoryMap;
	
	public:

	/**
	 * constructor
         \param np pointer to compound nucleus
	 */
	CHistory(CNucleus *np) {
        maxEvapZ = evap.maxZ;
		theMap[np]=0;
		tagDaughters(np, 0);
	}

	/**
	 * destructor
	 */
	~CHistory() {};

	/**
	 * returns the histroy of a fragment
        \param p is a pointer to the fragment
	*/
	int32_t getHistory(CNucleus *p) {
		HistoryMap::const_iterator iter = theMap.find(p);
		if(iter != theMap.end())
			return iter->second;
		else {
			std::cerr << "Unknown CNucleus pointer in HistoryMap" << std::endl;
			return 0;
		}
	}

	private:
	void tagDaughters(CNucleus *n, int32_t parentHistory);
    int32_t addToHistory(HistoryStepType steptype, int32_t prevhist);

	HistoryMap theMap;
	int maxEvapZ; //!< maximum Z for evaporation
	static CEvap evap; //!< class for evaporation of light particles
};
#endif // chistory_h

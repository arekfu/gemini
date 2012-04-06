#include "CHistory.h"

CEvap CHistory::evap;

void CHistory::tagDaughters(CNucleus *n, int32_t parentHistory) {
	CNucleus *daughterLight = n->getLightDaughter();
	CNucleus *daughterHeavy = n->getHeavyDaughter();

	if(!daughterLight && !daughterHeavy) {
		return;
	}

	int lightHistory, heavyHistory;

	// Saddle-to-scission transition or gamma decay
	if(!daughterLight) {
/*		if(!n->isSaddleToScission() && daughterHeavy->isSaddleToScission()) {
			// S2S
			heavyHistory = addToHistory(SaddleToScission, parentHistory);
			theMap[daughterHeavy] = heavyHistory;
			tagDaughters(daughterHeavy, heavyHistory);
		} else {*/
			// Gamma decay
			theMap[daughterHeavy] = parentHistory;
			tagDaughters(daughterHeavy, parentHistory);
//		}
		return;
	}


	if(n->isSaddleToScission()) {

		if(daughterLight->iZ <= maxEvapZ) {
			lightHistory = addToHistory(Evaporation, addToHistory(SaddleToScission, parentHistory));
			heavyHistory = addToHistory(SaddleToScission, parentHistory);
		} else {
			lightHistory = addToHistory(SymmetricFissionLight, parentHistory);
			heavyHistory = addToHistory(SymmetricFissionHeavy, parentHistory);
		}
		theMap[daughterLight] = lightHistory;
		theMap[daughterHeavy] = heavyHistory;
	} else {
		if(n->isNotStatistical()) {
			lightHistory = addToHistory(NonStatistical, parentHistory);
			heavyHistory = addToHistory(NonStatistical, parentHistory);
		} else if(daughterLight->iZ <= maxEvapZ) {
			lightHistory = addToHistory(Evaporation, parentHistory);
			heavyHistory = addToHistory(EvaporationResidue, parentHistory);
		} else {
			lightHistory = addToHistory(AsymmetricFissionLight, parentHistory);
			heavyHistory = addToHistory(AsymmetricFissionHeavy, parentHistory);
		}
		theMap[daughterLight] = lightHistory;
		theMap[daughterHeavy] = heavyHistory;
	}
	tagDaughters(daughterLight, lightHistory);
	tagDaughters(daughterHeavy, heavyHistory);
}

//************************************************************
/**
 * add a step to the history variable
 */
int32_t CHistory::addToHistory(HistoryStepType steptype, int32_t prevhist)
{

	int32_t history;

	int32_t maxhist=1;
	for( int i=0; i<std::numeric_limits<int32_t>::digits10; i++)
		maxhist *= 10;

	HistoryStepType last = (HistoryStepType) (prevhist % 10);

	if( last == EvaporationResidue && steptype == Evaporation )
		history = prevhist - EvaporationResidue + Evaporation;
/*	else if( last == SaddleToScission && steptype == Evaporation )
		history = prevhist - SaddleToScission + Evaporation;*/
	else if( last == EvaporationResidue && steptype == EvaporationResidue )
		history = prevhist;
	else if( last == SaddleToScission && steptype == SaddleToScission )
		history = prevhist;
	else if( prevhist < maxhist )
		history = prevhist*10 + steptype;
	else
	{
		history = (prevhist%maxhist)*10 + steptype;
		//      cerr << "Capacity of history variable exceeded!" << endl;
	}

	return history;
}


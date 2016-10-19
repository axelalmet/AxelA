/*

Copyright (c) 2005-2014, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "StochasticTargetRatioBasedCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "PanethCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "Exception.hpp"
#include "Debug.hpp"

StochasticTargetRatioBasedCellCycleModel::StochasticTargetRatioBasedCellCycleModel()
    : AbstractSimpleCellCycleModel(),
      mTargetRatioProbability(DOUBLE_UNSET),
      mMinCellCycleDuration(12.0), //Hours
      mMaxCellCycleDuration(14.0) //Hours
{
}

StochasticTargetRatioBasedCellCycleModel::StochasticTargetRatioBasedCellCycleModel(const StochasticTargetRatioBasedCellCycleModel& rModel)
: AbstractSimpleCellCycleModel(rModel),
  mTargetRatioProbability(rModel.mTargetRatioProbability),
  mMinCellCycleDuration(rModel.mMinCellCycleDuration),
  mMaxCellCycleDuration(rModel.mMaxCellCycleDuration)
{
 /*
  * Set each member variable of the new cell-cycle model that inherits
  * its value from the parent.
  *
  * Note 1: some of the new cell-cycle model's member variables will already
  * have been correctly initialized in its constructor or parent classes.
  *
  * Note 2: one or more of the new cell-cycle model's member variables
  * may be set/overwritten as soon as InitialiseDaughterCell() is called on
  * the new cell-cycle model.
  *
  * Note 3: Only set the variables defined in this class. Variables defined
  * in parent classes will be defined there.
  *
  */

 // No new member variables.
}

AbstractCellCycleModel* StochasticTargetRatioBasedCellCycleModel::CreateCellCycleModel()
{
    return new StochasticTargetRatioBasedCellCycleModel(*this);
}

/**
 * Overridden SetCellCycleDuration Method to add stochastic cell cycle times
 */
void StochasticTargetRatioBasedCellCycleModel::SetCellCycleDuration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mCellCycleDuration = DBL_MAX;
    }
    else
    {
        mCellCycleDuration = mMinCellCycleDuration + (mMaxCellCycleDuration - mMinCellCycleDuration) * p_gen->ranf(); // U[MinCCD,MaxCCD]
    }
}


double StochasticTargetRatioBasedCellCycleModel::GetTargetRatioProbability()
{
	return mTargetRatioProbability;
}

void StochasticTargetRatioBasedCellCycleModel::SetTargetRatioProbability(double targetRatioProbability)
{
	mTargetRatioProbability = targetRatioProbability;
}

/**
 * @return mMinCellCycleDuration
 */
double StochasticTargetRatioBasedCellCycleModel::GetMinCellCycleDuration()
{
	return mMinCellCycleDuration;
}

/**
 * Set mMinCellCycleDuration
 *
 * @param minCellCycleDuration
 */
void StochasticTargetRatioBasedCellCycleModel::SetMinCellCycleDuration(double minCellCycleDuration)
{
	mMinCellCycleDuration = minCellCycleDuration;
}

/**
 * @return mMaxCellCycleDuration
 */
double StochasticTargetRatioBasedCellCycleModel::GetMaxCellCycleDuration()
{
	return mMaxCellCycleDuration;
}

/**
 * Set mMaxCellCycleDuration
 *
 * @param maxCellCycleDuration
 */
void StochasticTargetRatioBasedCellCycleModel::SetMaxCellCycleDuration(double maxCellCycleDuration)
{
	mMaxCellCycleDuration = maxCellCycleDuration;
}

/**
 * Overridden GetAverageTransitCellCycleTime() method.
 * @return time
 */
double StochasticTargetRatioBasedCellCycleModel::GetAverageTransitCellCycleTime()
{
	return 0.5 * (mMinCellCycleDuration + mMaxCellCycleDuration);
}

/**
 * Overridden GetAverageStemCellCycleTime() method.
 * @return time
 */
double StochasticTargetRatioBasedCellCycleModel::GetAverageStemCellCycleTime()
{
	return 0.5 * (mMinCellCycleDuration + mMaxCellCycleDuration);
}

void StochasticTargetRatioBasedCellCycleModel::InitialiseDaughterCell()
{

    /*
     * We set the daughter cell to have a probability mTargetRatioProbability  a stem cell or
     * a paneth cell
     */

	//Get the target ratio probability
	double target_ratio_probability = GetTargetRatioProbability();

	//Get random number
	RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
	double uniform_random_number = p_gen->ranf();

	//Set daughter cell to be a transit cell
	boost::shared_ptr<AbstractCellProperty> p_transit_type =
			mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
	mpCell->SetCellProliferativeType(p_transit_type);

	if (uniform_random_number < target_ratio_probability) // set probability of becoming a Paneth cell
	{
		boost::shared_ptr<AbstractCellProperty> p_wildtype_state =
				mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<WildTypeCellMutationState>();
		mpCell->SetMutationState(p_wildtype_state);
	}
	else
	{
		boost::shared_ptr<AbstractCellProperty> p_paneth_state =
				mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<PanethCellMutationState>();
		mpCell->SetMutationState(p_paneth_state);
	}

	AbstractSimpleCellCycleModel::InitialiseDaughterCell();
}

void StochasticTargetRatioBasedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{

    *rParamsFile << "\t\t\t<TargetRatioProbability>" << mTargetRatioProbability << "</TargetRatioProbability>\n";
    *rParamsFile << "\t\t\t<MinCellCycleDuration>" << mMinCellCycleDuration << "</MinCellCycleDuration>\n";
    *rParamsFile << "\t\t\t<MaxCellCycleDuration>" << mMaxCellCycleDuration << "</MaxCellCycleDuration>\n";


    // Nothing to output, so just call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(StochasticTargetRatioBasedCellCycleModel)

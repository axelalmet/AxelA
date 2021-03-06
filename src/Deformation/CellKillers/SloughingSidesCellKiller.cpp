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
#include "SloughingSidesCellKiller.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "Exception.hpp"

template<unsigned DIM>
SloughingSidesCellKiller<DIM>::SloughingSidesCellKiller(AbstractCellPopulation<DIM>* pCrypt, double sloughLeftSide, double sloughRightSide)
    : AbstractCellKiller<DIM>(pCrypt),
      mCellsRemovedBySloughing(0),
      mSloughLeftSide(sloughLeftSide),
      mSloughRightSide(sloughRightSide)
{

    assert(sloughRightSide - sloughLeftSide > 0.0);
    mSloughLeftSide = sloughLeftSide;
    mSloughRightSide = sloughRightSide;
}

template<unsigned DIM>
unsigned SloughingSidesCellKiller<DIM>::GetNumberCellsRemoved() const
{
    return mCellsRemovedBySloughing;
}

template<unsigned DIM>
double SloughingSidesCellKiller<DIM>::GetSloughLeftSide() const
{
    return mSloughLeftSide;
}

template<unsigned DIM>
double SloughingSidesCellKiller<DIM>::GetSloughRightSide() const
{
    return mSloughRightSide;
}

template<unsigned DIM>
void SloughingSidesCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    switch (DIM)
    {
        case 1:
        {
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
                 cell_iter != this->mpCellPopulation->End();
                 ++cell_iter)
            {
                double x = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter)[0];

                if ( ((x < mSloughLeftSide) || (x > mSloughRightSide))&&(!cell_iter->GetCellProliferativeType()->template IsType<DifferentiatedCellProliferativeType>()) ) //If an epithelial cell has moved past the boundaries
                {
                    cell_iter->Kill();
                    mCellsRemovedBySloughing += 1;
                }
            }
            break;
        }
        case 2:
        {
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
                 cell_iter != this->mpCellPopulation->End();
                 ++cell_iter)
            {
                c_vector<double, 2> location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
                double x = location[0];

                if ( ((x < mSloughLeftSide) || (x > mSloughRightSide))&&(!cell_iter->GetCellProliferativeType()->template IsType<DifferentiatedCellProliferativeType>()) ) //If an epithelial cell has moved past the boundaries
                {
                    cell_iter->Kill();
                    mCellsRemovedBySloughing += 1;
                }
            }
            break;
        }
        case 3:
        {
            EXCEPTION("SloughingSidesCellKiller is not yet implemented in 3D");
            break;
        }
        default:
            // This can't happen
            NEVER_REACHED;
    }
}

template<unsigned DIM>
void SloughingSidesCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CellsRemovedBySloughing>" << mCellsRemovedBySloughing << "</CellsRemovedBySloughing>\n";
    *rParamsFile << "\t\t\t<SloughLeftSide>" << mSloughLeftSide << "</SloughLeftSide>\n";
    *rParamsFile << "\t\t\t<SloughRightSide>" << mSloughRightSide << "</SloughRightSide>\n";


    // Call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class SloughingSidesCellKiller<1>;
template class SloughingSidesCellKiller<2>;
template class SloughingSidesCellKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SloughingSidesCellKiller)

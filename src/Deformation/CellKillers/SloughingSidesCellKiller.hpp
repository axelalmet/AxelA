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

#ifndef SLOUGHINGSIDESCELLKILLER_HPP_
#define SLOUGHINGSIDESCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 *  A cell killer that kills cells if they are outside either of the vertical walls.
 *  The left and right sides are specified in the argument.
 */
template<unsigned DIM>
class SloughingSidesCellKiller : public AbstractCellKiller<DIM>
{
private:

    /** The number of cells sloughed **/
	unsigned mCellsRemovedBySloughing;

    /**
     * The left wall of the domain, non-dimensionalised by cell length.
     */
    double mSloughLeftSide;

    /**
     * The right wall of the domain, non-dimensionalised by cell length.
     */
    double mSloughRightSide;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     * @param number of cells removed by SloughingSidesSides
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<DIM> >(*this);
        archive & mCellsRemovedBySloughing;
    }

public:

    /**
     * Default constructor.
     *
     * @param pCrypt pointer to a cell population
     * @param sloughHeight the height at which to slough from the domain
     * @param sloughSides whether to slough cells at the side of the domain
     * @param sloughWidth the width of the domain (note slough on left and right)
     */
    SloughingSidesCellKiller(AbstractCellPopulation<DIM>* pCrypt,
                        double sloughLeftSide = 0.0,
                        double sloughRightSide = 10.0);

    /**
     * Destructor
     */
    virtual ~SloughingSidesCellKiller(){};

    /*
     * @return mCellsRemovedBySloughingSidesSides
     */

    unsigned GetNumberCellsRemoved() const;

    /**
     * @return mSloughLeftSide.
     */
    double GetSloughLeftSide() const;

    /*
     * @return mSloughRightSide
     */
    double GetSloughRightSide() const;

    /**
     * Loops over cells and kills cells outside boundary.
     */
    virtual void CheckAndLabelCellsForApoptosisOrDeath();

    /**
     * Outputs cell killer parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellKillerParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SloughingSidesCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a SloughingSidesSidesCellKiller.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const SloughingSidesCellKiller<DIM> * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_crypt = t->GetCellPopulation();
    ar << p_crypt;
    double slough_left_side = t->GetSloughLeftSide();
    ar << slough_left_side;
    double slough_right_side = t->GetSloughRightSide();
    ar << slough_right_side;
}

/**
 * De-serialize constructor parameters and initialise a SloughingSidesSidesCellKiller.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, SloughingSidesCellKiller<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_crypt;
    ar >> p_crypt;

    double slough_left_side;
    ar >> slough_left_side;
    double slough_right_side;
    ar >> slough_right_side;

    // Invoke inplace constructor to initialise instance
    ::new(t)SloughingSidesCellKiller<DIM>(p_crypt, slough_left_side, slough_right_side);
}
}
} // namespace ...

#endif /*SLOUGHINGSIDESCELLKILLER_HPP_*/

#ifndef ORGANOIDANOIKISCELLKILLER_HPP_
#define ORGANOIDANOIKISCELLKILLER_HPP_
#include "CheckpointArchiveTypes.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "AbstractCellKiller.hpp"

/*
 * Cell killer that removes any epithelial cell that has popped out from the ring
 * and entered the inner hole of ghost nodes. Adapted from SJD's
 * CrossSectionModelRandomCellKiller
 */

class OrganoidAnoikisCellKiller : public AbstractCellKiller<2>
{
private:

	// Number of cells removed by Anoikis
	unsigned mCellsRemovedByAnoikis;

    std::vector<c_vector<double,3> > mLocationsOfAnoikisCells;

    //Cut off radius for NodeBasedCellPopulations
    double mCutOffRadius;

    // The output file directory for the simulation data that corresponds to the number of cells
    // killed by anoikis
    out_stream mAnoikisOutputFile;

    std::string mOutputDirectory;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<2> >(*this);
        archive & mCellsRemovedByAnoikis;
        archive & mCutOffRadius;
        archive & mOutputDirectory;
    }

public:

    /**
     * Default constructor.
     *
     * @param pCellPopulation pointer to a tissue
     * @param sloughOrifice whether to slough compressed cells at crypt orifice
     */
	OrganoidAnoikisCellKiller(AbstractCellPopulation<2>* pCellPopulation);

	// Destructor
	~OrganoidAnoikisCellKiller();

    void SetOutputDirectory(std::string outputDirectory);

    std::string GetOutputDirectory();

    /*
     * @return mCutOffRadius
     */
    double GetCutOffRadius();

    /*
     * Method to defin mCutOffRadius by
     * cutOffRadius
     */
    void SetCutOffRadius(double cutOffRadius);

    std::set<unsigned> GetNeighbouringNodeIndices(unsigned nodeIndex);

    bool HasCellPoppedUp(unsigned nodeIndex);

    std::vector<c_vector<unsigned,2> > RemoveByAnoikis();

    /**
     *  Loops over and kills cells by anoikis or at the orifice if instructed.
     */
    void CheckAndLabelCellsForApoptosisOrDeath();

    /* After each event of cell killing in CheckAndLabelCellsForApoptosisOrDeath(), the information of whether to kill each cell
     * or not is passed to this method which then increments the member variables corresponding to the total number of cells
     * killed by anoikis or apoptosis through compression
     */
    void SetNumberCellsRemoved(std::vector<c_vector<unsigned,2> > cellsRemoved);

    /* Returns the total number of cells removed by anoikis ([0]) and by compression ([1])
     *
     */
    unsigned GetNumberCellsRemoved();

    /* Storing the x-locations of those epithelial cells that get removed by anoikis
     *
     */
    void SetLocationsOfCellsRemovedByAnoikis(std::vector<c_vector<unsigned,2> > cellsRemoved);

    /* Returns the x-coordinates of those cells removed by anoikis
     *
     */
    std::vector<c_vector<double,3> > GetLocationsOfCellsRemovedByAnoikis();

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
CHASTE_CLASS_EXPORT(OrganoidAnoikisCellKiller)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const OrganoidAnoikisCellKiller * t, const unsigned int file_version)
        {
            const AbstractCellPopulation<2>* const p_cell_population = t->GetCellPopulation();
            ar << p_cell_population;
        }

        template<class Archive>
        inline void load_construct_data(
            Archive & ar, OrganoidAnoikisCellKiller * t, const unsigned int file_version)
        {
            AbstractCellPopulation<2>* p_cell_population;
            ar >> p_cell_population;

            // Invoke inplace constructor to initialise instance
            ::new(t)OrganoidAnoikisCellKiller(p_cell_population);
        }
    }
}

#endif /* ORGANOIDANOIKISCELLKILLER_HPP_ */

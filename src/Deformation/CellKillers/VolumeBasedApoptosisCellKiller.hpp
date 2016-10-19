#ifndef VOLUMEBASEDAPOPTOSISCELLKILLER_HPP_
#define VOLUMEBASEDAPOPTOSISCELLKILLER_HPP_

#include "CheckpointArchiveTypes.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "AbstractCellKiller.hpp"

/*
 * Cell killer class that removes cells once they have become too compressed,
 * as determined by their cell volume. The killer only removes mature cells, in order
 * to maintain a cycling population of epithelial cells
 */

class VolumeBasedApoptosisCellKiller : public AbstractCellKiller<2>
{
private:

	// Number of cells removed by Apoptosis
	unsigned mCellsRemovedByApoptosis;

    std::vector<c_vector<double,3> > mLocationsOfAnoikisCells;

    //Cut off radius for NodeBasedCellPopulations
    double mCutOffRadius;

    //Threshold volume for determining apoptosis
    double mThresholdVolume;

    // The output file directory for the simulation data that corresponds to the number of cells
    // killed by anoikis
    out_stream mApoptosisOutputFile;

    std::string mOutputDirectory;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<2> >(*this);
        archive & mCellsRemovedByApoptosis;
        archive & mCutOffRadius;
        archive & mThresholdVolume;
        archive & mOutputDirectory;
    }

public:

    /**
     * Default constructor.
     *
     * @param pCellPopulation pointer to a tissue
     */
	VolumeBasedApoptosisCellKiller(AbstractCellPopulation<2>* pCellPopulation);

	// Destructor
	~VolumeBasedApoptosisCellKiller();

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

    /*
     * @return mThresholdVlume
     */
    double GetThresholdVolume();

    /*
     * Method to define mThresholdVolume by thresholdVolume
     */
    void SetThresholdVolume(double thresholdVolume);

    /**
     *  Loops over and kills cells by apoptosis
     */
    void CheckAndLabelCellsForApoptosisOrDeath();

    /*
     * Returns number of cells removed by apoptosis
     */
    unsigned GetNumberCellsRemoved();

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
CHASTE_CLASS_EXPORT(VolumeBasedApoptosisCellKiller)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const VolumeBasedApoptosisCellKiller * t, const  unsigned int file_version)
        {
            const AbstractCellPopulation<2>* const p_cell_population = t->GetCellPopulation();
            ar << p_cell_population;
        }

        template<class Archive>
        inline void load_construct_data(
            Archive & ar, VolumeBasedApoptosisCellKiller * t, const unsigned int file_version)
        {
            AbstractCellPopulation<2>* p_cell_population;
            ar >> p_cell_population;

            // Invoke inplace constructor to initialise instance
            ::new(t)VolumeBasedApoptosisCellKiller(p_cell_population);
        }
    }
}

#endif /* VOLUMEBASEDAPOPTOSISCELLKILLER_HPP_ */

/*
 * LAST MODIFIED: 21/12/2014
 * Cell killer to remove excessively-compressived mature epithelial cells
 *
 * Created on: Dec 21 2014
 * Last modified:
 * 			Author: Axel Almet
 */

#include "VolumeBasedApoptosisCellKiller.hpp"
#include "AbstractCellKiller.hpp"
#include "AbstractCellProperty.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "PanethCellMutationState.hpp"

VolumeBasedApoptosisCellKiller::VolumeBasedApoptosisCellKiller(AbstractCellPopulation<2>* pCellPopulation)
    : AbstractCellKiller<2>(pCellPopulation),
    mCellsRemovedByApoptosis(0),
    mCutOffRadius(1.5),
    mThresholdVolume(DOUBLE_UNSET)
{
    // Sets up output file
//	OutputFileHandler output_file_handler(mOutputDirectory + "ApoptosisData/", false);
//	mApoptosisOutputFile = output_file_handler.OpenOutputFile("results.Apoptosis");
}

VolumeBasedApoptosisCellKiller::~VolumeBasedApoptosisCellKiller()
{
//    mApoptosisOutputFile->close();
}

//Method to get mCutOffRadius
double VolumeBasedApoptosisCellKiller::GetCutOffRadius()
{
	return mCutOffRadius;
}

//Method to set mCutOffRadius
void VolumeBasedApoptosisCellKiller::SetCutOffRadius(double cutOffRadius)
{
	mCutOffRadius = cutOffRadius;
}

//Method to get mThresoldVolume
double VolumeBasedApoptosisCellKiller::GetThresholdVolume()
{
	return mThresholdVolume;
}

//Method to set mThresholdVolume
void VolumeBasedApoptosisCellKiller::SetThresholdVolume(double thresholdVolume)
{
	mThresholdVolume = thresholdVolume;
}

/*
 * Cell Killer that kills healthy cells that pop outwards and become detached from
 * the labelled tissue cells, i.e. removal by Apoptosis
 *
 * Also will remove differentiated cells at the orifice if mSloughOrifice is true
 */
void VolumeBasedApoptosisCellKiller::CheckAndLabelCellsForApoptosisOrDeath()
{
	double threshold_volume = GetThresholdVolume();

	if (dynamic_cast<MeshBasedCellPopulation<2>*>(this->mpCellPopulation))
	{
		MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*> (this->mpCellPopulation);
		//    assert(p_tissue->GetVoronoiTessellation()!=NULL);	// This fails during archiving of a simulation as Voronoi stuff not archived yet

    	for (AbstractCellPopulation<2>::Iterator cell_iter = p_tissue->Begin();
    			cell_iter != p_tissue->End();
    			++cell_iter)
    	{
    		//Only consider cells
    		unsigned node_index = p_tissue->GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
    		assert((!p_tissue->IsGhostNode(node_index)));

    		//Only consider epithelial cells
    		if (!cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    		{
    			double cell_age = cell_iter->GetAge();
    			double cell_volume = p_tissue->GetVolumeOfCell(*cell_iter);

    			if ( (cell_volume < threshold_volume)&&(cell_age > 1.0) ) //If the cell is an adult cell and is compressed sufficiently, kill it
    			{
    				cell_iter->Kill(); //Kill cell
    				mCellsRemovedByApoptosis += 1; //Update number of cells removed
    			}
    		}
    	}

	}
	else if (dynamic_cast<NodeBasedCellPopulation<2>*>(this->mpCellPopulation))
	{
		NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*> (this->mpCellPopulation);

    	for (AbstractCellPopulation<2>::Iterator cell_iter = p_tissue->Begin();
    			cell_iter != p_tissue->End();
    			++cell_iter)
    	{
    		//Only consider epithelial cells
    		if (!cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    		{
    			double cell_age = cell_iter->GetAge();
    			double cell_volume = p_tissue->GetVolumeOfCell(*cell_iter);

    			if ( (cell_volume < threshold_volume)&&(cell_age > 1.0) ) //If the cell is an adult cell and is compressed sufficiently, kill it
    			{
//    				std::cout << "Cell has been killed, cell volume = " << cell_volume << ", threshold = " << threshold_volume;
    				cell_iter->Kill(); //Kill cell
    				mCellsRemovedByApoptosis += 1; //Update number of cells removed
    			}
    		}
    	}
	}
}

unsigned VolumeBasedApoptosisCellKiller::GetNumberCellsRemoved()
{
	return mCellsRemovedByApoptosis;
}

void VolumeBasedApoptosisCellKiller::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CellsRemovedByApoptosis>" << mCellsRemovedByApoptosis << "</CellsRemovedByApoptosis> \n";
    *rParamsFile << "\t\t\t<CutOffRadius>" << mCutOffRadius << "</CutOffRadius> \n";

    // Call direct parent class
    AbstractCellKiller<2>::OutputCellKillerParameters(rParamsFile);
}




#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(VolumeBasedApoptosisCellKiller)

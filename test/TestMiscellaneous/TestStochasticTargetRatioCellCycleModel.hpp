/*
 * TestOrganoidRingWithTargetRatio.hpp
 *
 * Created on: April 08 2015
 * Last modified: April 08 2015
 * 		Author: Axel A. Almet
 */

#ifndef TESTSTOCHASTICTARGETRATIOCELLCYCLEMODEL_HPP_
#define TESTSTOCHASTICTARGETRATIOCELLCYCLEMODEL_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files
#include "CellBasedSimulationArchiver.hpp" //Needed if we would like to save/load simulations
#include "AbstractCellBasedTestSuite.hpp" //Needed for cell-based tests: times simulations, generates random numbers and has cell properties
#include "CheckpointArchiveTypes.hpp" //Needed if we use GetIdentifier() method (which we do)
#include "SmartPointers.hpp" //Enables macros to save typing
#include "CellsGenerator.hpp" //Generates cell population
#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "VoronoiDataWriter.hpp" //Allows us to visualise output in Paraview
#include "DifferentiatedCellProliferativeType.hpp" //Stops cells from proliferating
#include "StochasticTargetRatioBasedCellCycleModel.hpp" //Cell cycle model that will randomly assign daughter cells to become Paneth Cells
#include "VolumeTrackingModifier.hpp"
#include "TransitCellProliferativeType.hpp"
#include "PanethCellMutationState.hpp"
#include "CellLabel.hpp"
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "CellProliferativeTypesCountWriter.hpp" //Count cells that can proliferate
#include "OrganoidLinearSpringForce.hpp" //Spring force law to account for different node pairs
#include "OrganoidRingDataTrackingModifier.hpp" //Modifier for all the necessary data

#include "Debug.hpp"

class TestStochasticTargetRatioCellCycleModel : public AbstractCellBasedTestSuite
{
public:
	void TestTissueSlabWithTargetRatioCellCycle() throw(Exception)
	{
		//Simulation time parameters
		double dt = 0.01; //Set dt
		double end_time = 1.0; //Set end time
		double sampling_timestep = end_time/dt; //Set sampling timestep

		//Set all the spring stiffness variables
		double cell_cell_stiffness = 15.0;
		double cell_gel_stiffness = 15.0;
		double gel_gel_stiffness = 15.0;

		//Set the target proportions for stem : paneth cells
		double target_proportion_values [2] = {0.2, 0.8};

		//Set the domain of the model
		unsigned cells_across = 4; //Desired width + a few more layers, to avoid the box collapsing
		unsigned cells_up = 4; //Since height of each cell is 0.5*sqrt(3), we need to specify the no. of cells such that #cells*0.5*sqrt(3) + few more layers \approx desired height
		unsigned ghosts = 3;


		for (unsigned tratio = 0; tratio < 2; tratio ++)
		{

			//Set target proportion for stem cells to Paneth cells.
			double target_proportion = target_proportion_values[tratio];

			//				//Command line argument stuff: Get the index parameter
			//				unsigned index = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-myintval");
			//
			//To be extra careful, we reseed the random number generator
			RandomNumberGenerator::Instance()->Reseed(100);


			HoneycombMeshGenerator generator(cells_across, cells_up, ghosts);
			MutableMesh<2,2>* p_mesh = generator.GetMesh();

			//Get the non-ghost node locations
			std::vector<unsigned> real_indices = generator.GetCellLocationIndices();

			//This is for the purpose of tracking cells. Weird computer science stuff.
			boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
			boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
			boost::shared_ptr<AbstractCellProperty> p_paneth_state = CellPropertyRegistry::Instance()->Get<PanethCellMutationState>();
			boost::shared_ptr<AbstractCellProperty> p_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();

			//Create vector of cells
			std::vector<CellPtr> cells;

			//Create contact-inhibited cell cycle for each cell. Set them all to be differentiated.
			for (unsigned i = 0; i<real_indices.size(); i++)
			{
				//Set stochastic duration based cell cycle
				StochasticTargetRatioBasedCellCycleModel* p_cycle_model = new StochasticTargetRatioBasedCellCycleModel();
				p_cycle_model->SetDimension(2);
				p_cycle_model->SetTargetRatioProbability(target_proportion);
				double birth_time = 12.0*RandomNumberGenerator::Instance()->ranf(); //We would like the birth time to be ~U(0,12) and set in the past
				p_cycle_model->SetBirthTime(-birth_time);

				CellPtr p_cell(new Cell(p_state, p_cycle_model));
				p_cell->SetCellProliferativeType(p_stem_type);
				p_cell->InitialiseCellCycleModel();

				cells.push_back(p_cell);
			}

			//Create cell population
			MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, real_indices);


			//Track mutation states
			cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

			//For each cell in the tissue, we draw a random number and assign cells to be stem cells with
			//a probability of target_proportion, where target_proportion is a probability defined above.

//			for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
//					cell_iter != cell_population.End();
//					++cell_iter)
//			{
//				//Randomly generate number
//				double random_number = RandomNumberGenerator::Instance()->ranf();
//
//				unsigned node_index = cell_population.GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
//
//				if(!cell_population.IsGhostNode(node_index))
//				{
//					if(random_number >= target_proportion) //Assign cells to be Paneth with 1 - target_proportion
//					{
//						cell_iter->SetMutationState(p_paneth_state);
//					}
//				}
//			}


			//Allow output in Paraview (for now)
			cell_population.AddPopulationWriter<VoronoiDataWriter>();

			OffLatticeSimulation<2> simulator(cell_population);

			//Simulation pre-amble

			//Set output directory
			std::stringstream out;
			out << "/TRATIO_" << target_proportion;

			std::string output_directory = "TissueSlabWithTargetRatio/" + out.str() + "/";
			simulator.SetOutputDirectory(output_directory);

			simulator.SetOutputDivisionLocations(true); //Output the division locations and times of cells

			simulator.SetDt(dt);
			simulator.SetSamplingTimestepMultiple(sampling_timestep); //Sample the simulation at every hour
			simulator.SetEndTime(end_time); //Hopefully this is long enough for a steady state

			//Make the pointer for the modifier
			MAKE_PTR(OrganoidRingDataTrackingModifier<2>, p_data_tracking_modifier);
			//File reads as such:
			//[0] = Time, [1] = Ring area, [2] = Total area, [3] = Perimeter,
			//[4] = Stem, [5] = Transit, [6] = Differentiated, [7] = Default,
			//[8] = Healthy, [9] = APC One Hit, [10] = APC Two Hit, [11] = Beta Catenin.
			simulator.AddSimulationModifier(p_data_tracking_modifier);

			//Add linear spring force (modified to have three different spring stiffnesses, depending on the type of pair)
			MAKE_PTR(OrganoidLinearSpringForce<2>, p_spring_force);
			p_spring_force->SetCutOffLength(1.5);
			p_spring_force->SetCellCellSpringStiffness(cell_cell_stiffness); //Default is 15
			p_spring_force->SetCellGelSpringStiffness(cell_gel_stiffness); //Default is 15
			p_spring_force->SetGelGelSpringStiffness(gel_gel_stiffness); //Default is 15
			simulator.AddForce(p_spring_force);

			simulator.Solve();

			//Tidying up
			SimulationTime::Instance()->Destroy();
			SimulationTime::Instance()->SetStartTime(0.0);
		}
	}
};

#endif /* TESTSTOCHASTICTARGETRATIOCELLCYCLEMODEL_HPP_ */

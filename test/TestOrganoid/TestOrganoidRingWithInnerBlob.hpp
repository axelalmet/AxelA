/*
 * TestOrganoidRingWithInnerBlob.hpp
 *
 * Created on: Feb 27 2015
 * Last modified: Dec 21 2014
 * 		Author: Axel A. Almet
 */

#ifndef TESTORGANOIDRINGWITHINNERBLOB_HPP_
#define TESTORGANOIDRINGWITHINNERBLOB_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files
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
//#include "ContactInhibitionCellCycleModel.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "VolumeTrackingModifier.hpp"
#include "TransitCellProliferativeType.hpp"
#include "PanethCellMutationState.hpp"
#include "LumenCellMutationState.hpp"
#include "CellLabel.hpp"
#include "FixedRegionPlaneBoundaryCondition.hpp"
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "CellProliferativeTypesCountWriter.hpp" //Count cells that can proliferate
#include "OrganoidLinearSpringForce.hpp" //Spring force law to account for different node pairs
#include "OrganoidBasementMembraneForce.hpp" //Prototype BM force for ring based off SJD's code.
//#include "OrganoidAnoikisCellKiller.hpp" //Cell killer to remove proliferative cells that have detached from ring.
#include "OrganoidRingDataTrackingModifier.hpp" //Modifier for all the necessary data
#include "PetscSetupAndFinalize.hpp"
#include "OrganoidLumenModifier.hpp"

#include "Debug.hpp"

class TestOrganoidRingWithInnerBlob : public AbstractCellBasedTestSuite
{
public:
	void TestTissueRingWithInnerBlob() throw(Exception)
	{
		//Simulation time parameters
		double dt = 0.01; //Set dt
		double end_time = 1.0; //Set end time

		//Set all the spring stiffness variables
		double cell_cell_stiffness = 30.0;
		double cell_gel_stiffness = 30.0;
		double gel_gel_stiffness = 30.0;

		//Set the BM force parameters. As far as I know, we're only testing BM = 10, TC = 0.4 for Paneth cells
		double bm_force = 10.0; //Set the BM force parameter
		double target_curvature = 0.4; //Set the target curvature

		//Set the domain of the model
		unsigned cells_across = 24; //Desired width + a few more layers, to avoid the box collapsing
		unsigned cells_up = 27; //Since height of each cell is 0.5*sqrt(3), we need to specify the no. of cells such that #cells*0.5*sqrt(3) + few more layers \approx desired height
		unsigned ghosts = 0;

		//Translate mesh 1.5 units left and 1.5 units down, so that we have a layer of fixed cells around the boundary for the BC.
		c_vector<double, 2> translate_left = zero_vector<double>(2);
		translate_left(0) = -1.5;

		c_vector<double, 2> translate_down = zero_vector<double>(2);
		translate_down(1) = -1.5;

		//Define circle by centre and radius
		c_vector<double,2> circle_centre;
		circle_centre(0) = 10.5;
		circle_centre(1) = 10.0;

		double circle_radius = 4.5; //Size of hole
		assert(circle_radius > 0.0); //Just in case someone does something crazy.

		for (unsigned iter = 4; iter < 5; iter++)
		{
			//To be extra careful, we reseed the random number generator
			RandomNumberGenerator::Instance()->Reseed(100*iter);

			double sampling_timestep;
//			if(iter == 1) //Only need one realisation for visualisation---saves memory
//			{
				sampling_timestep = 0.5/dt; //Set sampling time---every hour, for the sake of visualisation and plots
//			}
//			else
//			{
//				sampling_timestep = end_time/dt;
//			}

			HoneycombMeshGenerator generator(cells_across, cells_up, ghosts);
			MutableMesh<2,2>* p_mesh = generator.GetMesh();

			//Translate mesh appropriately
			p_mesh->Translate(translate_left);
			p_mesh->Translate(translate_down);

			//Obtain vector of non-ghost nodes
			std::vector<unsigned> real_indices = generator.GetCellLocationIndices();


			//This is for the purpose of tracking cells. Weird computer science stuff.
			boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
			boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
			boost::shared_ptr<AbstractCellProperty> p_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();
			boost::shared_ptr<AbstractCellProperty> p_lumen_state = CellPropertyRegistry::Instance()->Get<LumenCellMutationState>();

			//Create vector of cells
			std::vector<CellPtr> cells;

			//Create stochastic based cell cycle for each cell. Set them all to be differentiated.
			for (unsigned i = 0; i<real_indices.size(); i++)
			{
				//Create the cell cycle
				StochasticDurationCellCycleModel* p_cycle_model = new StochasticDurationCellCycleModel();
				p_cycle_model->SetDimension(2);
				double birth_time = 12.0*RandomNumberGenerator::Instance()->ranf(); //We would like the birth time to be ~U(0,12) and set in the past
				p_cycle_model->SetBirthTime(-birth_time);

				CellPtr p_cell(new Cell(p_state, p_cycle_model));
				p_cell->SetCellProliferativeType(p_diff_type);
				p_cell->InitialiseCellCycleModel();

				cells.push_back(p_cell);
			}

			//Create cell population
			MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, real_indices);

			//Create blob of proliferative cells
	        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
	             cell_iter != cell_population.End();
	             ++cell_iter)
	        {
	            double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
	            double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

	            if (pow(x-circle_centre[0], 2) + pow(y-circle_centre[1], 2) <= pow(circle_radius,2))
				{
					cell_iter->SetCellProliferativeType(p_stem_type);
				}
	        }

//			//Iterate again and check that proliferative cells are also attached to gel nodes.
//			//If they aren't, set them to be a LumenMutationState and thus terminally differentiated.
//			for (unsigned i = 0; i < real_indices.size(); i++)
//			{
//				unsigned cell_index = real_indices[i];
//				CellPtr cell_iter = cell_population.GetCellUsingLocationIndex(cell_index);
//				double x = cell_population.GetLocationOfCellCentre(cell_iter)[0];
//				double y = cell_population.GetLocationOfCellCentre(cell_iter)[1];
//
//				Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(cell_iter);
//				unsigned node_index = p_node->GetIndex();
//
//				if (cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>() == false)
//				{
//					bool element_contains_gel_nodes = false;
//
//					//Iterate over elements (triangles) containing the node
//					for (Node<2>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
//							iter != p_node->ContainingElementsEnd();
//							++iter)
//					{
//						// Get a pointer to the element
//						Element<2,2>* p_element = cell_population.rGetMesh().GetElement(*iter);
//
//						// Check if its triangulation contains a gel node
//						for (unsigned local_index=0; local_index<3; local_index++)
//						{
//							unsigned nodeGlobalIndex = p_element->GetNodeGlobalIndex(local_index);
//							bool is_ghost_node = cell_population.IsGhostNode(nodeGlobalIndex);
//
//							if (is_ghost_node == false) //Make sure we're not dealing with ghost nodes (otherwise this stuff will fail)
//							{
//								CellPtr p_local_cell = cell_population.GetCellUsingLocationIndex(nodeGlobalIndex);
//								if (p_local_cell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()==true)
//								{
//									element_contains_gel_nodes = true;
//									break; 				// This should break out of the inner for loop
//								}
//							}
//						}
//					}
//
//					if(element_contains_gel_nodes == false)
//					{
//						cell_iter->SetMutationState(p_lumen_state); //Make cell a "lumen" cell
//						cell_iter->GetCellCycleModel()->SetTransitCellG1Duration(DBL_MAX); //Make cell terminally differentiated
//					}
//				}
//			}

			//Allow output in Paraview (for now)
			cell_population.AddPopulationWriter<VoronoiDataWriter>();

			OffLatticeSimulation<2> simulator(cell_population);

			//Simulation pre-amble

			//Set output directory
			std::stringstream out;
			out << "BM_" << bm_force << "_CURV_" << target_curvature;
			std::stringstream run;
			run << 100*iter;
			std::string output_directory = "OrganoidRingWithInnerBlob/" + out.str() + "/" + run.str(); //Set main directory (to be re-used later when running through Paneth cell configs)
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

			//Make pointer for lumen creating modifier
			MAKE_PTR(OrganoidLumenModifier<2>, p_lumen_modifier);
			//Add modifier to simulation
			simulator.AddSimulationModifier(p_lumen_modifier);

			//Add linear spring force (modified to have three different spring stiffnesses, depending on the type of pair)
			MAKE_PTR(OrganoidLinearSpringForce<2>, p_spring_force);
			p_spring_force->SetCutOffLength(1.5);
			p_spring_force->SetCellCellSpringStiffness(cell_cell_stiffness); //Default is 15
			p_spring_force->SetCellGelSpringStiffness(cell_gel_stiffness); //Default is 15
			p_spring_force->SetGelGelSpringStiffness(gel_gel_stiffness); //Default is 15 //Default is 1
			simulator.AddForce(p_spring_force);

			//Add basement membrane force
			MAKE_PTR(OrganoidBasementMembraneForce, p_bm_force);
			p_bm_force->SetBasementMembraneParameter(bm_force); //Equivalent to beta in SJD's papers
			p_bm_force->SetRingCurvature(target_curvature); //This is equivalent to 1/R in SJD's papers
			p_bm_force->SetPositionDependentMultiplier(false);	// Multiplier for basement membrane force that may one day be needed
			simulator.AddForce(p_bm_force);

			//Impose a box on the cells by adding BCs.
			c_vector<double,2> point = zero_vector<double>(2);
			c_vector<double,2> normal = zero_vector<double>(2);

			//Impose a wall at x > 0
			normal(0) = -1.0;
			MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
			simulator.AddCellPopulationBoundaryCondition(p_bc1);

			//Impose a right wall
			point(0) = 20.0;
			normal(0) = 1.0;
			MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal));
			simulator.AddCellPopulationBoundaryCondition(p_bc2);

			//Impose a wall at y > 0
			point(0) = 0.0;
			point(1) = 0.0;
			normal(0) = 0.0;
			normal(1) = -1.0;
			MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal));
			simulator.AddCellPopulationBoundaryCondition(p_bc3);

			//Impose the upper wall
			point(1) = 20.0;
			normal(1) = 1.0;
			MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal));
			simulator.AddCellPopulationBoundaryCondition(p_bc4);

			simulator.Solve();
		}
	}
};


#endif /* TESTORGANOIDRINGWITHINNERBLOB_HPP_ */

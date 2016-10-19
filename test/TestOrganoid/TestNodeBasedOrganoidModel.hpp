/* TestNodeBasedOrganoidModel.hpp
 *
 * A test for implementing SJD's basement membrane model
 * using an overlapping spheres based model, extended to
 * an organoid geometry.
 *
 * Created on: Mar 02 2015
 * Last modified: Mar 02 2016
 * 			Author: Axel A Almet
 */

#ifndef TESTNODEBASEDORGANOIDMODEL_HPP_
#define TESTNODEBASEDORGANOIDMODEL_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files
#include "AbstractCellBasedTestSuite.hpp" //Needed for cell-based tests: times simulations, generates random numbers and has cell properties
#include "CheckpointArchiveTypes.hpp" //Needed if we use GetIdentifier() method (which we do)
#include "SmartPointers.hpp" //Enables macros to save typing
#include "PetscSetupAndFinalize.hpp"
#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "NodeBasedCellPopulation.hpp" //Node-based CellPopulation class
#include "SmartPointers.hpp" //Enables macros to save typing
#include "OrganoidLinearSpringForce.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "UniformlyDistributedCellCycleModel.hpp"
#include "VoronoiDataWriter.hpp" //Allows us to visualise output in Paraview
#include "OrganoidRingDataTrackingModifier.hpp"
#include "PanethCellMutationState.hpp"
#include "FakePetscSetup.hpp" //Forbids tests running in parallel (for now)
#include "OrganoidLinearSpringForce.hpp" //Linear spring force with variable spring constant definition for different cell types
#include "NodeBasedBasementMembraneForceWithNearestNeighbours.hpp" //Basement membrane force that defines curvature using two closest epithelial neighbours
#include "OrganoidAnoikisCellKiller.hpp" //Cell killer to remove cells if they lose contact with basement membrane
#include "VolumeBasedApoptosisCellKiller.hpp"
#include "FixedRegionPlaneBoundaryCondition.hpp" //Region that pins down cells when they move past a certain point
#include "Debug.hpp"

class TestNodeBasedOrganoidModel : public AbstractCellBasedTestSuite
{
public:
	void TestNodeOrganoidModel() throw(Exception)
	{
		//Set simulation time parameters
		double dt = 0.005; //Set dt
		double end_time = 1.0; //Set end time
		double sampling_timestep = end_time/dt; //Set sampling time step

		//Set radius of connectivity
		double radius_of_interaction = 1.5;

		//Set division separation
		double division_separation = 0.1;

		//Set spring stiffness constants
		double cell_cell_stiffness = 15.0;
		double cell_stroma_stiffness = 15.0;
		double stroma_stroma_stiffness = 15.0;

		//Set the basement membrane parameters
		double basement_membrane_stiffness = 4.0;
		double target_curvature = 0.4;

		//Set the number of cells across and down
		unsigned cells_across = 24;
		unsigned cells_down = 24;

		//Set the radius of the initial organoid cyst
		double circle_radius = 2.0;

		//Define the centre of the circle
		c_vector<double, 2> circle_centre;
		circle_centre[0] = 10.0;
		circle_centre[1] = 8.0;

		//Translate mesh 1.5 units left and 1.5 units down, so that we have a layer of fixed cells around the boundary for the BC.
		c_vector<double, 2> translate_left = zero_vector<double>(2);
		translate_left(0) = -1.5;

		c_vector<double, 2> translate_down = zero_vector<double>(2);
		translate_down(1) = -1.5;

		unsigned run = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-myintval");

		for (unsigned index = run; index < run + 10; index++)
		{
			//To be extra careful, we reseed the random number generator
			RandomNumberGenerator::Instance()->Reseed(100*index);

			for (double iter = 0.0; iter < 10.0; iter ++)
			{
				//Set the threshold volume for apoptosis
				double threshold_volume = iter*0.1*0.25*M_PI; //Set it to be 0.1*iter*equilibrium volume (i.e. 0.1area, 0.2area etc)

				//Get the initial honeycomb (equilibrium) mesh
				HoneycombMeshGenerator generator(cells_across, cells_down);
				MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();

				//Translate mesh appropriately
				p_generating_mesh->Translate(translate_left);
				p_generating_mesh->Translate(translate_down);

				//We now define the actual mesh by defining a set of nodes and their corresponding locations
				std::vector<Node<2>*> nodes;

				//Get the indices of all the cells
				std::vector<unsigned> initial_cell_indices = generator.GetCellLocationIndices();

				//Iterate over each cell in the initial mesh
				for (unsigned i = 0; i < initial_cell_indices.size(); i++)
				{
					unsigned cell_index = initial_cell_indices[i]; //Get cell index

					//Get the location
					double x = p_generating_mesh->GetNode(cell_index)->rGetLocation()[0];
					double y = p_generating_mesh->GetNode(cell_index)->rGetLocation()[1];

					//If the node falls outside of the hole, add it to the vector of nodes
					if(pow(x - circle_centre(0), 2.0) + pow(y - circle_centre(1), 2.0) > pow(circle_radius, 2.0))
					{
						Node<2>* cell_node = p_generating_mesh->GetNode(cell_index); //Get node
						nodes.push_back(cell_node); //Add node to vector
					}
				}

				//Construct the actual nde mesh now
				NodesOnlyMesh<2> mesh;
				mesh.ConstructNodesWithoutMesh(nodes, 2.0*radius_of_interaction);

				//Create shared pointers for cell and mutation states
				boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
				boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
				boost::shared_ptr<AbstractCellProperty> p_wildtype_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();

				//Create vector of cells
				std::vector<CellPtr> cells;

				//Create contact-inhibited cell cycle for each cell. Set them all to be differentiated.
				for (unsigned i = 0; i<mesh.GetNumNodes(); i++)
				{
					//Set stochastic duration based cell cycle
					UniformlyDistributedCellCycleModel* p_cycle_model = new UniformlyDistributedCellCycleModel();
					p_cycle_model->SetDimension(2);
					double birth_time = 12.0*RandomNumberGenerator::Instance()->ranf(); //We would like the birth time to be ~U(0,12) and set in the past
					p_cycle_model->SetBirthTime(-birth_time);
					//			p_cycle_model->SetTransitCellG1Duration(DBL_MAX);

					CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
					p_cell->SetCellProliferativeType(p_diff_type); //Make cell differentiated
					p_cell->InitialiseCellCycleModel();

					cells.push_back(p_cell);
				}

				NodeBasedCellPopulation<2> cell_population(mesh, cells);//Create the cell population

				//Set the division separation
				cell_population.SetMeinekeDivisionSeparation(division_separation);

				//We now need t define the 'ring' of epithelial cells (oh boy, total fluke if this works)
				for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
						cell_iter != cell_population.End(); ++cell_iter)
				{
					double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0]; //Get the x-coordinate of the node
					double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1]; //Get the y-coordinate of the node

					if(pow(x - circle_centre(0), 2.0) + pow(y - circle_centre(1), 2.0) < pow(circle_radius + 1.0, 2.0) ) //If the cell lies close enough to the ring then make it an epithelial cell
					{
						cell_iter->SetCellProliferativeType(p_stem_type); //Set cell to be a proliferative type
					}

				}

				//Explicitly request output for Paraview
				cell_population.AddPopulationWriter<VoronoiDataWriter>();

				OffLatticeSimulation<2> simulator(cell_population); //Define the simulation class

				//Simulation pre-amble
				std::stringstream params;
				params << "BM_" << basement_membrane_stiffness << "_CURV_" << target_curvature;
				std::stringstream run;
				run << 100*index;
				std::stringstream threshold;
				threshold << 0.1*iter << "V";
				std::string output_directory = "NodeBasedOrganoidModel/"+ params.str() + "/" + run.str()
											+ "/" + threshold.str(); //Set the output directory
				simulator.SetOutputDirectory(output_directory);

				simulator.SetDt(dt); //Set dt
				simulator.SetSamplingTimestepMultiple(sampling_timestep); //Set when to sample the simulation
				simulator.SetEndTime(end_time); //Set the end time

				//Add basement membrane force with nearest neighbours
				MAKE_PTR(NodeBasedBasementMembraneForceWithNearestNeighbours, p_bm_force);
				p_bm_force->SetBasementMembraneParameter(basement_membrane_stiffness); //Equivalent to beta in SJD's papers
				p_bm_force->SetTargetCurvature(target_curvature); //This is equivalent to 1/R in SJD's papers
				p_bm_force->SetPositionDependentMultiplier(false);	//Multiplier for basement membrane force that may one day be needed
				//		p_bm_force->SetLeftCryptBoundary(left_crypt_boundary);
				//		p_bm_force->SetRightCryptBoundary(right_crypt_boundary);
				p_bm_force->SetCryptGeometry(false);
				p_bm_force->SetCutOffRadius(radius_of_interaction); //Set cut off radius for defining neighbours
				simulator.AddForce(p_bm_force);

				//Set the generalised linear spring force law with variable spring constant
				MAKE_PTR(OrganoidLinearSpringForce<2>, p_variable_constant_force);
				p_variable_constant_force->SetCellCellSpringStiffness(cell_cell_stiffness); //Set stiffness between epithelial cells
				p_variable_constant_force->SetCellGelSpringStiffness(cell_stroma_stiffness); //Set stiffness between epithelial cells and stromal cells
				p_variable_constant_force->SetGelGelSpringStiffness(stroma_stroma_stiffness); //Set stiffness between stromal cells
				p_variable_constant_force->SetCutOffLength(radius_of_interaction);
				p_variable_constant_force->SetMeinekeDivisionRestingSpringLength(division_separation);
				simulator.AddForce(p_variable_constant_force);

				//Add anoikis-based cell killer
				MAKE_PTR_ARGS(OrganoidAnoikisCellKiller, p_anoikis_killer, (&cell_population));
				p_anoikis_killer->SetCutOffRadius(1.5); //Set the cut off radius for the cell killer
				simulator.AddCellKiller(p_anoikis_killer);

				//Add volume-based cell killer
				MAKE_PTR_ARGS(VolumeBasedApoptosisCellKiller, p_apoptosis_killer, (&cell_population));
				p_apoptosis_killer->SetThresholdVolume(threshold_volume); //Set the cut off radius for the cell killer
				simulator.AddCellKiller(p_apoptosis_killer);

				//We want to fix the outer layers of the organoid
				//Impose a box on the cells by adding BCs.
				c_vector<double,2> point = zero_vector<double>(2);
				c_vector<double,2> normal = zero_vector<double>(2);

				//Fix the rows of cells below y = 0.0
				normal(1) = -1.0;
				MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bottom_bc, (&cell_population, point, normal));
				simulator.AddCellPopulationBoundaryCondition(p_bottom_bc);

				//Fix the rows of cells above y = 17.0
				point(1) = 17.0;
				normal(1) = 1.0;
				MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_top_bc, (&cell_population, point, normal));
				simulator.AddCellPopulationBoundaryCondition(p_top_bc);

				//Fox the rows of cells to the left of x = 0.0
				point(1) = 0.0;
				normal(0) = -1.0;
				normal(1) = 0.0;
				MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_left_bc, (&cell_population, point, normal));
				simulator.AddCellPopulationBoundaryCondition(p_left_bc);

				//Fox the rows of cells to the right of x = 20.0
				point(0) = 20.0;
				normal(0) = 1.0;
				MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_right_bc, (&cell_population, point, normal));
				simulator.AddCellPopulationBoundaryCondition(p_right_bc);

				simulator.Solve(); //Solve the simulation

				//			for (AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
				//					cell_iter != simulator.rGetCellPopulation().End(); ++cell_iter)
				//			{
				//				if (!cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>() )
				//				{
				//					double volume = simulator.rGetCellPopulation().GetVolumeOfCell(*cell_iter);
				//					PRINT_VARIABLE(volume);
				//				}
				//			}

				//Tidying up
				SimulationTime::Instance()->Destroy();
				SimulationTime::Instance()->SetStartTime(0.0);

			}
		}
	}
};

#endif /*TESTNODEBASEDORGANOIDMODEL_HPP_*/

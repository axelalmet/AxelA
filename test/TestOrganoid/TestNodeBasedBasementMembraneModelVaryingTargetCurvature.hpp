\/* TestNodeBasedBasementMembraneModel.hpp
 *
 * A test for implementing SJD's basement membrane model
 * using an overlapping spheres based model.
 *
 * Created on: Mar 05 2015
 * Last modified: Feb 29 2016
 * 			Author: Axel A Almet
 */

#ifndef TESTNODEBASEDBASEMENTMEMBRANEMODELVARYINGTARGETCURVATURE_HPP_
#define TESTNODEBASEDBASEMENTMEMBRANEMODELVARYINGTARGETCURVATURE_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files
#include "AbstractCellBasedTestSuite.hpp" //Needed for cell-based tests: times simulations, generates random numbers and has cell properties
#include "CheckpointArchiveTypes.hpp" //Needed if we use GetIdentifier() method (which we do)
#include "SmartPointers.hpp" //Enables macros to save typing
#include "PetscSetupAndFinalize.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp" //Generates mesh with periodic vertical BCs
#include "Cylindrical2dNodesOnlyMesh.hpp"
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "NodeBasedCellPopulation.hpp" //Node-based CellPopulation class
#include "SmartPointers.hpp" //Enables macros to save typing
#include "OrganoidLinearSpringForce.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "UniformDistributedCellCycleModel.hpp"
#include "VoronoiDataWriter.hpp" //Allows us to visualise output in Paraview
#include "PanethCellMutationState.hpp"
#include "FakePetscSetup.hpp" //Forbids tests running in parallel (for now)
#include "OrganoidLinearSpringForce.hpp" //Linear spring force with variable spring constant definition for different cell types
#include "NodeBasedBasementMembraneForceWithNearestNeighbours.hpp" //Basement membrane force that defines curvature using two closest epithelial neighbours
#include "OrganoidAnoikisCellKiller.hpp" //Cell killer to remove cells if they lose contact with BM
#include "FixedRegionPlaneBoundaryCondition.hpp"
//#include "CryptBoxBoundaryCondition.hpp" //Box BC that pins cells on the bottom and applies plane BC on the left and right
#include "SloughingSidesCellKiller.hpp" //Removes cells when they go past the vertical edges
//#include "VolumeBasedApoptosisCellKiller.hpp" //Removes cells if they are sufficiently compressed

#include "Debug.hpp"

class TestNodeBasedBasementMembraneModelVaryingTargetCurvature : public AbstractCellBasedTestSuite
{
public:
	void TestNodeBoxModelVaryingTargetCurvature() throw(Exception)
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

		//Set the crypt boundaries
		double left_crypt_boundary = 5.0;
		double right_crypt_boundary = 15.0;

		//Set the left and right crypt sides for sloughing
		double left_slough_side =  0.25; //Set this boundary to be negative so that cells don't get sloughed on the left.
		double right_slough_side = 19.75;

		//Set the domain width for apoptosis from the sides and the minimum volume
		double region_width = 2.0;
		double minimum_volume = 0.75*M_PI*pow(0.5, 2.0);

		//Define the boundaries
		c_vector<double, 2> left_point = zero_vector<double>(2);
		left_point[0] = left_slough_side;

		c_vector<double, 2> right_point = zero_vector<double>(2);
		right_point[0] = right_slough_side;

		c_vector<double, 2> bottom_point = zero_vector<double>(2);
		bottom_point[1] = 0.5;

		//Set the number of cells across and down for the initial mesh
		unsigned cells_across = 20;
		unsigned cells_down = 10;

		//Set the translation factors for the mesh
		c_vector<double, 2> translate_down = zero_vector<double>(2);
		translate_down(1) = 0.0;

		c_vector<double, 2> translate_left = zero_vector<double>(2);
		translate_left(0) = 0.0;

		unsigned start_run = 0;
		unsigned end_run = 5;

		//Command line argument stuff for the target curvature parameters
//		double target_curvature = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-mydoubleval");
		double target_curvature = 0.3;
		double basement_membrane_parameter = 5.0;

		for (unsigned index = start_run; index < end_run; index++)
		{
			double target_curvature = 0.1 * index;

			//Command line argument stuff: Get the index parameter
//			unsigned index = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-myintval");

			//To be extra careful, we reseed the random number generator
			RandomNumberGenerator::Instance()->Reseed(100*index);

			HoneycombMeshGenerator generator(cells_across, cells_down, 0); //Create mesh
			MutableMesh<2,2>* p_generating_mesh = generator.GetMesh(); //Generate mesh

			//Create some nodes
			std::vector<Node<2>*> nodes;
			//Translate mesh
			p_generating_mesh->Translate(translate_left);
			p_generating_mesh->Translate(translate_down);

			//Construct the nodes only mesh using the MutableMesh
			Cylindrical2dNodesOnlyMesh mesh(20.0);
//			mesh.ConstructNodesWithoutMesh(nodes, 2.0*radius_of_interaction);
			mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 20.0); //Construct mesh

			//Create shared pointers for cell and mutation states
			boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
			boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
			boost::shared_ptr<AbstractCellProperty> p_wildtype_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();

			//Create tissue of cells. Initially we set them all to be differentiated
			std::vector<CellPtr> cells; //Create vector of cells

			for (unsigned i = 0; i<mesh.GetNumNodes(); i++)
			{
				//Set stochastic duration based cell cycle
				UniformlyDistributedCellCycleModel* p_cycle_model = new UniformlyDistributedCellCycleModel();
				p_cycle_model->SetDimension(2);
				double birth_time = 13.0*RandomNumberGenerator::Instance()->ranf(); //We would like the birth time to be ~U(0,12) and set in the past
				p_cycle_model->SetBirthTime(-birth_time);

				CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
				p_cell->SetCellProliferativeType(p_diff_type); //Make cell differentiated
				p_cell->InitialiseCellCycleModel();

				cells.push_back(p_cell);
			}

			NodeBasedCellPopulation<2> cell_population(mesh, cells);//Create the cell population

			//Set the division separation
			cell_population.SetMeinekeDivisionSeparation(division_separation);

			//Create the monolayer of proliferative cells

			//First we find the maximum height of the nodes in the tissue
			double max_height = 0.0;

			for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
					cell_iter != cell_population.End(); ++cell_iter)
			{
				double height = cell_population.GetLocationOfCellCentre(*cell_iter)[1]; //Get the y-coordinate of the node

				if(height > max_height)
				{
					max_height = height;
				}

			}

			//Mark any cell centre that has the same height as max_height as a proliferative cell
			for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
					cell_iter != cell_population.End(); ++cell_iter)
			{
				double height = cell_population.GetLocationOfCellCentre(*cell_iter)[1]; //Get the y-coordinate of the node

				if(height == max_height) //If the y-coordinate is equal to the max_height (it should be as the initial configuration is a grid)
				{
					cell_iter->SetCellProliferativeType(p_stem_type); //Set cell to be a proliferative type
				}

			}

			//Explicitly request output for Paraview
			cell_population.AddPopulationWriter<VoronoiDataWriter>();

			OffLatticeSimulation<2> simulator(cell_population); //Define the simulation class

			//Simulation pre-amble
			std::stringstream params;
			params << "BM_" << basement_membrane_parameter << "_CURV_" << target_curvature;
			std::stringstream run;
			run << 100*index;
			std::string output_directory = "NodeBasedBoxModelVaryingTC/" + run.str() + "/" + params.str(); //Set the output directory
			simulator.SetOutputDirectory(output_directory);

			simulator.SetDt(dt); //Set dt
			simulator.SetSamplingTimestepMultiple(sampling_timestep); //Set when to sample the simulation
			simulator.SetEndTime(end_time); //Set the end time


			//Set the generalised linear spring force law with variable spring constant
			MAKE_PTR(OrganoidLinearSpringForce<2>, p_variable_constant_force);
			p_variable_constant_force->SetCellCellSpringStiffness(cell_cell_stiffness); //Set stiffness between epithelial cells
			p_variable_constant_force->SetCellGelSpringStiffness(cell_stroma_stiffness); //Set stiffness between epithelial cells and stromal cells
			p_variable_constant_force->SetGelGelSpringStiffness(stroma_stroma_stiffness); //Set stiffness between stromal cells
			p_variable_constant_force->SetCutOffLength(radius_of_interaction);
			p_variable_constant_force->SetMeinekeDivisionRestingSpringLength(division_separation);
			simulator.AddForce(p_variable_constant_force);

			//Add basement membrane force with nearest neighbours
			MAKE_PTR(NodeBasedBasementMembraneForceWithNearestNeighbours, p_bm_force);
			p_bm_force->SetBasementMembraneParameter(basement_membrane_parameter); //Equivalent to beta in SJD's papers
			p_bm_force->SetTargetCurvature(target_curvature); //This is equivalent to 1/R in SJD's papers
			p_bm_force->SetPositionDependentMultiplier(false);	//Multiplier for basement membrane force that may one day be needed
			p_bm_force->SetLeftCryptBoundary(left_crypt_boundary);
			p_bm_force->SetRightCryptBoundary(right_crypt_boundary);
			p_bm_force->SetCutOffRadius(radius_of_interaction); //Set cut off radius for defining neighbours
			simulator.AddForce(p_bm_force);

			//Add anoikis-based cell killer
			MAKE_PTR_ARGS(OrganoidAnoikisCellKiller, p_anoikis_killer, (&cell_population));
			p_anoikis_killer->SetCutOffRadius(radius_of_interaction); //Set the cut off radius for the cell killer
			simulator.AddCellKiller(p_anoikis_killer);

			//Add volume-based apoptosis cell killer
//			MAKE_PTR_ARGS(VolumeBasedApoptosisCellKiller<2>, p_apoptosis_killer, (&cell_population, left_slough_side, right_slough_side));
//			p_apoptosis_killer->SetApoptosisRegionWidth(region_width);
//			p_apoptosis_killer->SetMinimumVolumeForApoptosis(minimum_volume);
//			simulator.AddCellKiller(p_apoptosis_killer);

			//Fix the bottom row of stromal cells
			c_vector<double, 2> normal = zero_vector<double>(2);
			normal(1) = -1.0;
			MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bottom_bc, (&cell_population, bottom_point, normal));
			simulator.AddCellPopulationBoundaryCondition(p_bottom_bc);

			//Add a sloughing cell killer for the sides
			MAKE_PTR_ARGS(SloughingSidesCellKiller<2>, p_sloughing_killer, (&cell_population, left_slough_side, right_slough_side)); //Set sloughSlides = true
			simulator.AddCellKiller(p_sloughing_killer);

			simulator.Solve(); //Solve the simulation

			/*
			 * Collect all the relevant data from the simulation.
			 */

			//Get the epithelial cells from the basement membrane force
			std::vector<unsigned> epithelial_indices = p_bm_force->GetEpithelialIndices(simulator.rGetCellPopulation());

			//Need to sort the cells in order of increasing x-coordinate
			std::vector<std::pair<double, unsigned> > x_coordinates_and_indices;

			for (unsigned i = 0; i < epithelial_indices.size(); i++)
			{
				//Get the index of the cell
				unsigned epithelial_index = epithelial_indices[i];

				//Get the x-coordinate
				double x = simulator.rGetCellPopulation().GetNode(epithelial_index)->rGetLocation()[0];

				std::pair<double, unsigned> x_coordinate_and_index = std::make_pair(x, epithelial_index);

				x_coordinates_and_indices.push_back(x_coordinate_and_index);
			}

			//Sort the vector by the x-coordinate
			std::sort(x_coordinates_and_indices.begin(), x_coordinates_and_indices.end());

			//We output the locations of the epithelial cells into a file

			OutputFileHandler results_handler(output_directory + "/", false); //Create output file handler
			std::string monolayer_filename = "epithelialmonolayerdata.dat";
			out_stream epithelial_monolayer_file = results_handler.OpenOutputFile(monolayer_filename);

			//Output the locations of the epithelial cells to the file
			for (unsigned i = 0; i < x_coordinates_and_indices.size(); i++)
			{
				unsigned epithelial_node_index = x_coordinates_and_indices[i].second;

				double x = simulator.rGetCellPopulation().GetNode(epithelial_node_index)->rGetLocation()[0];
				double y = simulator.rGetCellPopulation().GetNode(epithelial_node_index)->rGetLocation()[1];

				*epithelial_monolayer_file << x << "\t" << y << "\n";
			}
//
//			//We now calculate the average horizontal distance, average vertical distance and average gap ratio.
//			double average_horizontal_distance = 0.0;
//			double average_vertical_distance = 0.0;
//			double average_gap_ratio = 0.0;
//
//			//Do a bunch of stuff
//			for (unsigned i = 0; i < x_coordinates_and_indices.size() - 1; i++)
//			{
//				//Get the pair of consecutive indices
//				unsigned epithelial_index_A = x_coordinates_and_indices[i].second;
//				unsigned epithelial_index_B = x_coordinates_and_indices[i+1].second;
//
//				//Get the coordinates of index A
//				double x_A = simulator.rGetCellPopulation().GetNode(epithelial_index_A)->rGetLocation()[0];
//				double y_A = simulator.rGetCellPopulation().GetNode(epithelial_index_A)->rGetLocation()[1];
//
//				//Get the coordinates of index B
//				double x_B = simulator.rGetCellPopulation().GetNode(epithelial_index_B)->rGetLocation()[0];
//				double y_B = simulator.rGetCellPopulation().GetNode(epithelial_index_B)->rGetLocation()[1];
//
//				//Calculate the increment in horizontal and vertical distance and the gap ratio
//				double horizontal_increment = fabs(x_B - x_A);
//				double vertical_increment = fabs(y_B - y_A);
//				double gap_ratio_increment = vertical_increment/horizontal_increment;
//
//				//Update the three quantities
//				average_horizontal_distance += horizontal_increment;
//				average_vertical_distance += vertical_increment;
//				average_gap_ratio += gap_ratio_increment;
//
//			}
//
//			//Aveerage the three quantities
//			average_horizontal_distance /= (epithelial_indices.size() - 1);
//			average_vertical_distance /= (epithelial_indices.size() - 1);
//			average_gap_ratio /= (epithelial_indices.size() - 1);
//
//			//Get the number of epithelial cells in the monolayer
//			unsigned number_epithelial_cells = epithelial_indices.size();
//
//			//Get the number of cells removed by anoikis
//			unsigned number_cells_removed_by_anoikis = p_anoikis_killer->GetNumberCellsRemoved();
//
//			//Get the number of cells removd by sloughing
//			unsigned number_cells_removed_by_sloughing = p_sloughing_killer->GetNumberCellsRemoved();
//
////			//Create results file for model
////			OutputFileHandler results_handler(output_directory + "/", false); //Create output file handler
//
//			//Create file
//			std::string box_model_filename = "boxmodeldata.dat";
//			out_stream box_model_file = results_handler.OpenOutputFile(box_model_filename);
//
//			/*File will have the  following columns:
//			 * [0] = average horizontal distance between epithelial cells, [1] = average vertical distance
//			 * [2] = average gap ratio, [3] = number of epithelial cells, [4] = number of cells removed by anoikis
//			 * [5] = number of cells removed by sloughing
//			 */
//
//			*box_model_file << average_horizontal_distance << "\t" << average_vertical_distance << "\t" << average_gap_ratio
//					<< "\t" << number_epithelial_cells << "\t" << number_cells_removed_by_anoikis << "\t"
//					<< number_cells_removed_by_sloughing << "\n";

			//Tidying up
			SimulationTime::Instance()->Destroy();
			SimulationTime::Instance()->SetStartTime(0.0);

		}
	}
};

#endif /*TESTNODEBASEDBASEMENTMEMBRANEMODELVARYINGTARGETCURVATURE_HPP_*/

/*

 * TestOrganoidRingWithTargetRatio.hpp
 *
 * Created on: Mar 14 2015
 * Last modified: Mar 14 2015
 * 		Author: Axel A. Almet
 */

#ifndef TESTORGANOIDRINGWITHTARGETRATIO_HPP_
#define TESTORGANOIDRINGWITHTARGETRATIO_HPP_

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
#include "FixedRegionPlaneBoundaryCondition.hpp"
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "CellProliferativeTypesCountWriter.hpp" //Count cells that can proliferate
#include "OrganoidLinearSpringForce.hpp" //Spring force law to account for different node pairs
#include "EpithelialLayerBasementMembraneForce.hpp" //Prototype BM force for ring based off SJD's code.
#include "OrganoidAnoikisCellKiller.hpp" //Cell killer to remove proliferative cells that have detached from ring.
#include "OrganoidRingDataTrackingModifier.hpp" //Modifier for all the necessary data

#include "Debug.hpp"

class TestOrganoidRingWithTargetRatio : public AbstractCellBasedTestSuite
{
public:
	void TestTissueRingWithTargetRatio() throw(Exception)
	{
		//Simulation time parameters
		double dt = 0.005; //Set dt
		double end_time = 1.0; //Set end time
		double sampling_timestep = end_time/dt; //Set sampling timestep

		//Set all the spring stiffness variables
		double cell_cell_stiffness = 15.0;
		double cell_gel_stiffness = 15.0;
		double gel_gel_stiffness = 15.0;

		//Set the BM force parameters. As far as I know, we're only testing BM = 10, TC = 0.4 for Paneth cells
		double bm_force = 10.0; //Set the BM force parameter
		double target_curvature = 0.2; //Set the target curvature

		//Set the domain of the model
		unsigned cells_across = 24; //Desired width + a few more layers, to avoid the box collapsing
		unsigned cells_up = 27; //Since height of each cell is 0.5*sqrt(3), we need to specify the no. of cells such that #cells*0.5*sqrt(3) + few more layers \approx desired height
		unsigned ghosts = 4;

		//Translate mesh 1.5 units left and 1.5 units down, so that we have a layer of fixed cells around the boundary for the BC.
		c_vector<double, 2> translate_left = zero_vector<double>(2);
		translate_left(0) = -1.5;

		c_vector<double, 2> translate_down = zero_vector<double>(2);
		translate_down(1) = -1.5;

		//Define circle by centre and radius
		c_vector<double,2> circle_centre;
		circle_centre(0) = 10.5;
		circle_centre(1) = 10.0;

		double mutant_drag_multiplier = 1.0; //Relative adhesiveness of Paneth cells

		double circle_radius = 2.5; //Size of hole
		assert(circle_radius > 0); //Just in case someone does something crazy.

		double ring_radius = circle_radius + 2.0; //Radius of the ring of cells. This isn't the actual radius, just has to be large enough for later
		assert((ring_radius <= cells_across)&&(ring_radius <= cells_up)); //Again, just in case.

		//Command line argument stuff: Get the index parameter
		// unsigned index = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-myintval");
		unsigned index = 1;

		//To be extra careful, we reseed the random number generator
		RandomNumberGenerator::Instance()->Reseed(100*index);

		for (double pstiff = 0.0; pstiff < 1.0; pstiff ++)
		{
			for (double tratio = 0.0; tratio < 1.0; tratio ++)
			{

				//Set the Paneth cell spring stiffness multiplier and the target proportion of
				//stem cells to Paneth cells.
				double paneth_stiffness_multiplier = 4.5 + 0.25*pstiff;
				double target_proportion = 0.3 + 0.05*tratio;

				HoneycombMeshGenerator generator(cells_across, cells_up, ghosts);
				MutableMesh<2,2>* p_mesh = generator.GetMesh();

				//Translate mesh appropriately
				p_mesh->Translate(translate_left);
				p_mesh->Translate(translate_down);

				//Obtain vector of non-ghost nodes
				std::vector<unsigned> initial_real_indices = generator.GetCellLocationIndices();

				//Create another vector that will actually be passed into MeshBasedCellPopulationWithGhostNodes.
				//This defines our hole of ghost nodes in the middle.
				std::vector<unsigned> real_indices;

				for (unsigned i = 0; i < initial_real_indices.size(); i++)
				{
					unsigned cell_index = initial_real_indices[i];
					double x = p_mesh->GetNode(cell_index)->rGetLocation()[0];
					double y = p_mesh->GetNode(cell_index)->rGetLocation()[1];

					if (pow(x-circle_centre[0],2) + pow(y-circle_centre[1],2) > pow(circle_radius,2))
					{
						real_indices.push_back(cell_index);
					}
				}

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
					p_cycle_model->SetTargetRatioProbability(target_proportion);
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

				//Track contact inhibited cells
				cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

				//Create ring of proliferative cells
				for (unsigned i = 0; i < real_indices.size(); i++)
				{
					unsigned cell_index = real_indices[i];
					CellPtr cell_iter = cell_population.GetCellUsingLocationIndex(cell_index);
					double x = cell_population.GetLocationOfCellCentre(cell_iter)[0];
					double y = cell_population.GetLocationOfCellCentre(cell_iter)[1];

					//"un-differentiate" the inner cells attached to the ghost nodes.
					if (pow(x-circle_centre[0],2) + pow(y-circle_centre[1],2) <= pow(ring_radius,2))
					{
						Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(cell_iter);
						unsigned node_index = p_node->GetIndex();

						//Iterate over the elements (triangles) containing the nodes
						for (Node<2>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
								iter != p_node->ContainingElementsEnd();
								++iter)
						{
							bool element_contains_ghost_nodes = false;

							// Get a pointer to the element
							Element<2,2>* p_element = cell_population.rGetMesh().GetElement(*iter);

							// Check whether it's triangulation contains a ghost nod
							for (unsigned local_index=0; local_index<3; local_index++)
							{
								unsigned nodeGlobalIndex = p_element->GetNodeGlobalIndex(local_index);

								if (cell_population.IsGhostNode(nodeGlobalIndex) == true)
								{
									element_contains_ghost_nodes = true;
									break; 				// This should break out of the inner for loop
								}
							}

							if(element_contains_ghost_nodes)
							{
								cell_iter->SetCellProliferativeType(p_stem_type);
							}
						}
					}
				}

				//Iterate again and check that proliferative cells are also attached to gel nodes.
				//If they aren't, kill the cell.
				for (unsigned i = 0; i < real_indices.size(); i++)
				{
					unsigned cell_index = real_indices[i];
					CellPtr cell_iter = cell_population.GetCellUsingLocationIndex(cell_index);
					double x = cell_population.GetLocationOfCellCentre(cell_iter)[0];
					double y = cell_population.GetLocationOfCellCentre(cell_iter)[1];

					//Only consider this inside the ring
					if (pow(x-circle_centre[0],2) + pow(y-circle_centre[1],2) <= pow(ring_radius,2))
					{
						Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(cell_iter);
						unsigned node_index = p_node->GetIndex();

						//Only iterate over the initial ring of transit cells
						if (cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>() == false)
						{
							bool element_contains_gel_nodes = false;

							//Iterate over elements (triangles) containing the node
							for (Node<2>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
									iter != p_node->ContainingElementsEnd();
									++iter)
							{
								// Get a pointer to the element
								Element<2,2>* p_element = cell_population.rGetMesh().GetElement(*iter);

								// Check if its triangulation contains a gel node
								for (unsigned local_index=0; local_index<3; local_index++)
								{
									unsigned nodeGlobalIndex = p_element->GetNodeGlobalIndex(local_index);
									bool is_ghost_node = cell_population.IsGhostNode(nodeGlobalIndex);

									if (is_ghost_node == false) //Make sure we're not dealing with ghost nodes (otherwise this stuff will fail)
									{
										CellPtr p_local_cell = cell_population.GetCellUsingLocationIndex(nodeGlobalIndex);
										if (p_local_cell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()==true)
										{
											element_contains_gel_nodes = true;
											break; 				// This should break out of the inner for loop
										}
									}
								}
							}

							if(element_contains_gel_nodes == false)
							{
								cell_iter->Kill();
							}
						}
					}
				}

				//Randomly assign cells int he ring to be Paneth cells. But first, we need to find
				//which cells make up the ring.
				std::vector<unsigned> cells_in_ring; //Initialise vector

				//Obtain the proliferative cells
				for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
						cell_iter != cell_population.End();
						++cell_iter)
				{

					unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);

					//If the cell is within the previously defined ring (narrows our search down)
					//and it is a stem cell
					if (!cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>() )
					{
						cells_in_ring.push_back(node_index); //Add the angle and node index
					}
				}


				//For each cell in the ring, we draw a random number and assign cells to be stem cells with
				//a probability of target_proportion, where target_proportion is a probability defined above.

				for (unsigned i = 0; i < cells_in_ring.size(); i++)
				{
					unsigned node_index = cells_in_ring[i];

					CellPtr cell = cell_population.GetCellUsingLocationIndex(node_index);

					//Randomly generate number
					double random_number = RandomNumberGenerator::Instance()->ranf();

					if(random_number >= target_proportion) //Assign cells to be Paneth with 1 - target_proportion
					{
						cell->SetMutationState(p_paneth_state);
					}
				}

				//Allow output in Paraview (for now)
				cell_population.AddPopulationWriter<VoronoiDataWriter>();

				//We set Paneth cells to be 1.3 times as adhesive as stem cells
				double normal_damping_constant = cell_population.GetDampingConstantNormal();
				cell_population.SetDampingConstantMutant(mutant_drag_multiplier*normal_damping_constant);

				OffLatticeSimulation<2> simulator(cell_population);

				//Simulation pre-amble

				//Set output directory
				std::stringstream out;
				out << "DRAG_" << mutant_drag_multiplier;
				out << "/BM_" << bm_force << "_CURV_" << target_curvature;
				std::stringstream params;
				params << "PSTIFF_" << paneth_stiffness_multiplier << "/TRATIO_" << target_proportion;
				std::stringstream run;
				run << 100*index;
				std::string output_directory = "OrganoidRingWithTargetRatio/" + out.str() + "/"
						+ run.str() + "/" + params.str();
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
				p_spring_force->SetPanethCellStiffnessMultiplier(paneth_stiffness_multiplier); //Default value is 1
				simulator.AddForce(p_spring_force);

				//Add basement membrane force
				MAKE_PTR(EpithelialLayerBasementMembraneForce, p_bm_force);
				p_bm_force->SetBasementMembraneParameter(bm_force); //Equivalent to beta in SJD's papers
				p_bm_force->SetTargetCurvature(target_curvature); //This is equivalent to 1/R in SJD's papers
				simulator.AddForce(p_bm_force);

				//Add anoikis-based cell killer
				MAKE_PTR_ARGS(OrganoidAnoikisCellKiller, p_anoikis_killer, (&cell_population));
				simulator.AddCellKiller(p_anoikis_killer);

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

				//Impose a wall at y < 0
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

				// //Create results file for anoikis events
				// OutputFileHandler results_handler(output_directory + "/", false); //Create output file handler

				// std::string anoikis_events_filename = "anoikisevents_steadystate.dat";
				// out_stream anoikis_events_file = results_handler.OpenOutputFile(anoikis_events_filename);

				// std::vector<c_vector<double,3> > locations_anoikis_cells = p_anoikis_killer->GetLocationsOfCellsRemovedByAnoikis(); //Get the anoikis events

				// for (unsigned k=0; k<locations_anoikis_cells.size(); k++)
				// {
				// 	*anoikis_events_file << locations_anoikis_cells[k][0] << "\t" << locations_anoikis_cells[k][1] << "\t" << locations_anoikis_cells[k][2] << "\n";
				// }

				// //Get the cells in the ring and the centre of mass

				// std::vector<unsigned> indices_of_cells_in_ring; //Initialise vector

				// c_vector<double,2> centre_of_mass; //Initialise centre of mass

				// unsigned num_epithelial_cells = 0; //Initialise count

				// for (AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
				// 		cell_iter != simulator.rGetCellPopulation().End(); ++cell_iter)
				// {
				// 	//Only consider epithelial cells
				// 	if(!cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>() )
				// 	{
				// 		//Get node inex
				// 		unsigned node_index = simulator.rGetCellPopulation().GetLocationIndexUsingCell(*cell_iter);

				// 		//Get location of cell
				// 		c_vector<double,2> cell_location = simulator.rGetCellPopulation().GetLocationOfCellCentre(*cell_iter);

				// 		indices_of_cells_in_ring.push_back(node_index); //Add node index

				// 		centre_of_mass += cell_location; //Update the centre of mass

				// 		num_epithelial_cells += 1; //Update count of epithelial cells in ring
				// 	}
				// }

				// //Normalise centre_of_mass
				// centre_of_mass /= num_epithelial_cells;

				// //Get the angle of each epithelial cell with respect to the centre of mass. That way we can sort the list
				// //in counter-clockwise order before we create the file.

				// //Initialise vector of angles and the cell indices
				// std::vector<std::pair<double, unsigned> > angles_and_cell_indices;

				// for (unsigned k = 0; k < indices_of_cells_in_ring.size(); k++)
				// {
				// 	unsigned cell_index = indices_of_cells_in_ring[k]; //Get the node index

				// 	CellPtr cell = simulator.rGetCellPopulation().GetCellUsingLocationIndex(cell_index); //Get the cell

				// 	//Get the cell location
				// 	double x = simulator.rGetCellPopulation().GetLocationOfCellCentre(cell)[0];
				// 	double y = simulator.rGetCellPopulation().GetLocationOfCellCentre(cell)[1];

				// 	//Make the pair of the angle and cell
				// 	std::pair<double, unsigned> angle_cell;

				// 	//Get point relative to the circle centre
				// 	double rel_x = x - centre_of_mass[0];
				// 	double rel_y = y - centre_of_mass[1];

				// 	double circle_angle = atan(rel_y/rel_x); //Get initial angle argument

				// 	if (rel_x<0.0) //If the point is in the second quadrant or third quadrant
				// 	{
				// 		circle_angle += M_PI;
				// 	}
				// 	else if ((rel_x>=0.0)&&(rel_y<0.0)) //Fourth quadrant
				// 	{
				// 		circle_angle += 2*M_PI;
				// 	}

				// 	angle_cell = std::make_pair(circle_angle, cell_index);

				// 	angles_and_cell_indices.push_back(angle_cell); //Add the angle and node index
				// }

				// //Sort the vector of cell indices by their angle
				// std::sort(angles_and_cell_indices.begin(), angles_and_cell_indices.end());

				// //Create file to store the location of the cells in the ring and their cell type

				// std::string cells_in_ring_filename = "cellsinringdata.dat";
				// out_stream cells_in_ring_file = results_handler.OpenOutputFile(cells_in_ring_filename);

				// //Create another file to store the areas of the cells in the ring and their cell type

				// std::string areas_of_cells_in_ring_filename = "areascellsinring.dat";
				// out_stream areas_of_cells_in_ring_file = results_handler.OpenOutputFile(areas_of_cells_in_ring_filename);

				// //Add each cell location to the file, along with a 1 if it is a stem cell and 0 if it is a Paneth cell
				// for (unsigned k=0; k<angles_and_cell_indices.size(); k++)
				// {
				// 	std::pair<double, unsigned> angle_and_cell_index = angles_and_cell_indices[k]; //Get pair
				// 	unsigned cell_index = angle_and_cell_index.second; //Get cell index

				// 	CellPtr cell = simulator.rGetCellPopulation().GetCellUsingLocationIndex(cell_index); //Get the cell

				// 	//Get the area
				// 	double cell_area = simulator.rGetCellPopulation().GetVolumeOfCell(cell);

				// 	//Get the cell location
				// 	double x = simulator.rGetCellPopulation().GetLocationOfCellCentre(cell)[0];
				// 	double y = simulator.rGetCellPopulation().GetLocationOfCellCentre(cell)[1];

				// 	unsigned is_stem_cell = 1;

				// 	if(!cell->GetMutationState()->IsType<PanethCellMutationState>() ) //Stem cell if 0, Paneth cell if 1
				// 	{
				// 		is_stem_cell = 0;
				// 	}

				// 	*cells_in_ring_file << x << "\t" << y << "\t" << is_stem_cell << "\n";

				// 	*areas_of_cells_in_ring_file << cell_area << "\t" << is_stem_cell << "\n";
				// }

				//Tidying up
				SimulationTime::Instance()->Destroy();
				SimulationTime::Instance()->SetStartTime(0.0);
			}
		}
	}

};

#endif /* TESTORGANOIDRINGWITHTARGETRATIO_HPP_ */

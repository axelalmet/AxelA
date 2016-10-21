
/*
 * TestOrganoidRing.hpp
 *
 * Created on: Dec 18 2014
 * Last modified: Feb 27 2015
 * 		Author: Axel A. Almet
 */

#ifndef TESTORGANOIDRING_HPP_
#define TESTORGANOIDRING_HPP_

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
#include "UniformCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "PanethCellMutationState.hpp"
#include "CellLabel.hpp"
#include "FixedRegionPlaneBoundaryCondition.hpp"
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "CellProliferativeTypesCountWriter.hpp" //Count cells that can proliferate
#include "OrganoidLinearSpringForce.hpp" //Spring force law to account for different node pairs
#include "EpithelialLayerBasementMembraneForce.hpp" //Prototype BM force for ring based off SJD's code.
#include "OrganoidAnoikisCellKiller.hpp" //Cell killer to remove proliferative cells that have detached from ring.
#include "EpithelialLayerDataTrackingModifier.hpp" //Modifier for all the necessary data
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

class TestOrganoidRing : public AbstractCellBasedTestSuite
{
public:
	void TestTissueRing() throw(Exception)
	{
		for (double BM = 2.0; BM < 3.0; BM++)
		{
			for (double TC = 1.0; TC < 2.0; TC++)
			{

				//Command line argument stuff

				//Get the index parameter
			  //	unsigned index = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-myintval");
			        unsigned index = 1;

				//To be extra careful, we reseed the random number generator
				RandomNumberGenerator::Instance()->Reseed(100*index);

				double dt = 0.01; //Set dt
				double end_time = 1.0; //Set end time
				double sampling_timestep = 0.5/dt; //Set sampling timestep multiple for Visualisations

				//Set all the spring stiffness variables
				double cell_cell_stiffness = 15.0;
				double cell_gel_stiffness = 15.0;
				double gel_gel_stiffness = 15.0;
				double paneth_stiffness_multiplier = 1.0;

				double bm_force = 5.0*BM; //Set the BM force parameter
				double target_curvature = 0.2*TC; //Set the target curvature

				unsigned cells_across = 24; //Desired width + a few more layers, to avoid the box collapsing
				unsigned cells_up = 27; //Since height of each cell is 0.5*sqrt(3), we need to specify the no. of cells such that #cells*0.5*sqrt(3) + few more layers \approx desired height
				unsigned ghosts = 4;

				HoneycombMeshGenerator generator(cells_across, cells_up, ghosts);
				MutableMesh<2,2>* p_mesh = generator.GetMesh();

				//Translate mesh 1.5 units left and 1.5 units down, so that we have a layer of fixed cells around the boundary for the BC.
				c_vector<double, 2> translate_left;
				translate_left(0) = -1.5;
				translate_left(1) = 0.0;

				c_vector<double, 2> translate_down;
				translate_down(0) = 0.0;
				translate_down(1) = -1.5;

				p_mesh->Translate(translate_left);
				p_mesh->Translate(translate_down);

				//Define circle by centre and radius
				c_vector<double,2> circle_centre;
				circle_centre(0) = 10.5;
				circle_centre(1) = 10.0;

				double circle_radius = 2.5; //Size of hole
				assert(circle_radius > 0); //Just in case someone does something crazy.

				double ring_radius = circle_radius + 2.0; //Radius of the ring of cells. This isn't the actual radius, just has to be large enough for later
				assert((ring_radius <= cells_across)&&(ring_radius <= cells_up)); //Again, just in case.

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
				boost::shared_ptr<AbstractCellProperty> p_transit_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
				boost::shared_ptr<AbstractCellProperty> p_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();
				boost::shared_ptr<AbstractCellProperty> p_paneth_type = CellPropertyRegistry::Instance()->Get<PanethCellMutationState>();

				//Create vector of cells
				std::vector<CellPtr> cells;

				//Create contact-inhibited cell cycle for each cell. Set them all to be differentiated.
				for (unsigned i = 0; i<real_indices.size(); i++)
				{
					//Set stochastic duration based cell cycle
					UniformCellCycleModel* p_cycle_model = new UniformCellCycleModel();
					p_cycle_model->SetDimension(2);
					double birth_time = 12.0*RandomNumberGenerator::Instance()->ranf(); //We would like the birth time to be ~U(0,12) and set in the past (this is probably not necessary now)
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


				//Create ring of stem cells
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
						
						//Iterate over the elements (triangles) containing the nodes
						for (Node<2>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
								iter != p_node->ContainingElementsEnd();
								++iter)
						{
							bool element_contains_ghost_nodes = false;

							// Get a pointer to the element
							Element<2,2>* p_element = cell_population.rGetMesh().GetElement(*iter);

							// Check whether its triangulation contains a ghost node
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
								cell_iter->SetCellProliferativeType(p_transit_type);
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

				OffLatticeSimulation<2> simulator(cell_population);

				//Simulation pre-amble

				//Set output directory
				std::stringstream out;
				out << 100*index << "/BM_" << bm_force << "_CURV_" << target_curvature;
				std::string output_directory = "OrganoidRing/" + out.str();
				simulator.SetOutputDirectory(output_directory);

				simulator.SetOutputDivisionLocations(true); //Output the division times and locations

				simulator.SetDt(dt);
				simulator.SetSamplingTimestepMultiple(sampling_timestep); //Sample the simulation at every hour
				simulator.SetEndTime(end_time); //Hopefully this is long enough for a steady state

				//Make the pointer for the modifier
				MAKE_PTR(EpithelialLayerDataTrackingModifier<2>, p_data_tracking_modifier);
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
				p_spring_force->SetPanethCellStiffnessMultiplier(paneth_stiffness_multiplier); //Default is 1
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

				//Create results file for anoikis events
				OutputFileHandler results_handler(output_directory + "/	", false); //Create output file handlere

				std::string anoikis_events_filename = "anoikisevents.dat";
				out_stream anoikis_events_file = results_handler.OpenOutputFile(anoikis_events_filename);

				std::vector<c_vector<double,3> > locations_anoikis_cells = p_anoikis_killer->GetLocationsOfCellsRemovedByAnoikis(); //Get the anoikis events

				for (unsigned k=0; k<locations_anoikis_cells.size(); k++)
				{
					*anoikis_events_file << locations_anoikis_cells[k][0] << "\t" << locations_anoikis_cells[k][1] << "\t" << locations_anoikis_cells[k][2] << "\n";
				}

				//Stuff needed to properly loop over simulations
				SimulationTime::Instance()->Destroy();
				SimulationTime::Instance()->SetStartTime(0.0);
			}
		}
	}
};

#endif /* TESTORGANOIDRING_HPP_ */

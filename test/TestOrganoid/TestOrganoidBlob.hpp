/*
 * TestAxelOrganoid.hpp
 *
 *  Last modified on: Dec 18 2014
 *  Created on: Oct 1, 2014
 *      Author: axelalmet
 */

#ifndef TESTORGANOIDBLOB_HPP_
#define TESTORGANOIDBLOB_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files
#include "AbstractCellBasedTestSuite.hpp" //Needed for cell-based tests: times simulations, generates random numbers and has cell properties
#include "CheckpointArchiveTypes.hpp" //Needed if we use GetIdentifier() method (which we do)
#include "SmartPointers.hpp" //Enables macros to save typing
#include "CellsGenerator.hpp" //Generates cell population
#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "GeneralisedLinearSpringForce.hpp" //Spring forces between cells
#include "VoronoiDataWriter.hpp" //Allows us to visualise output in Paraview
#include "DifferentiatedCellProliferativeType.hpp" //Stops cells from proliferating
#include "AbstractCellProliferativeType.hpp"
#include "ContactInhibitionCellCycleModel.hpp"
#include "VolumeTrackingModifier.hpp"
#include "TransitCellProliferativeType.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "CellProliferativeTypesCountWriter.hpp" //Count cells that can proliferate
#include "AxelLinearSpringForce.hpp" //Prototype spring force law to simulate the Matrigel forces

class TestOrganoidBlob : public AbstractCellBasedTestSuite
{
public:
	void TestSquareTissue() throw(Exception)
	{
		for (int iter = 1; iter <= 1; iter++)
		{
			//Iterate over different gel-gel spring stiffnesses
			for (double i = 0.0; i <= 1.0; i++)
			{
				HoneycombMeshGenerator generator(20,20,4); //Cells across, cells up
				MutableMesh<2,2>* p_mesh = generator.GetMesh(); //Obtains mesh from HoneycombMeshGenerator

				//Get the real nodes
				std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

				//Make the cells differentiated, so that they don't proliferate
				std::vector<CellPtr> cells;

				boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
				boost::shared_ptr<AbstractCellProperty> p_transit_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
				boost::shared_ptr<AbstractCellProperty> p_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();

				//Creates contact inhibition cell cycle for each cell
				for (unsigned i=0; i<location_indices.size(); i++)
				{
					double uniform_random_number = RandomNumberGenerator::Instance()->ranf();
					ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
					p_cycle_model->SetDimension(2);
					p_cycle_model->SetBirthTime(-2.0*(double)i);
					//Set equilibrium volume to a random number between 0.5 and 1
					p_cycle_model->SetQuiescentVolumeFraction(0.5 + 0.5*uniform_random_number);
					p_cycle_model->SetEquilibriumVolume(1.0);

					CellPtr p_cell(new Cell(p_state, p_cycle_model));
					p_cell->SetCellProliferativeType(p_transit_type);
					p_cell->InitialiseCellCycleModel();

					cells.push_back(p_cell);
				}

				MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
				//Track contact inhibited cells
				cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

				//Make all the cells outside a circle of specified radius differentiated.
				for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
						cell_iter != cell_population.End();
						++cell_iter)
				{
					double x = cell_population.GetLocationOfCellCentre(*cell_iter)[0];
					double y = cell_population.GetLocationOfCellCentre(*cell_iter)[1];

					//Start off with circle of cells
					if (pow(x-10,2) + pow(y-8.5,2) > pow(5,2))
					{

						cell_iter->SetCellProliferativeType(p_diff_type);
					}

				}

				//Track the cell proliferative types, so as to see when all the proliferating cells become quiescent.
				//[0] = Stem, [1] = Transit, [2] = Differentiated, [3] = Default.
				cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
				cell_population.AddPopulationWriter<VoronoiDataWriter>();

				//Set the Cell-Gel spring stiffness
				double gel_gel_stiffness = pow(2,i)*10.0;

				//Simulates evolution of population
				OffLatticeSimulation<2> simulator(cell_population);

				//Set the output directory
				std::stringstream out;
				out << "GGStiff_" << gel_gel_stiffness << "_" << iter;
				std::string output_directory = "AxelOrganoid/" + out.str();
				simulator.SetOutputDirectory(output_directory);

				simulator.SetDt(0.001);
				simulator.SetSamplingTimestepMultiple(10);
				simulator.SetEndTime(1.0);

				//Tracks volumes for contact inhibition
				MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
				simulator.AddSimulationModifier(p_modifier);

				//Add linear spring force (modified to have three different spring stiffnesses, depending on the type of pair)
				MAKE_PTR(AxelLinearSpringForce<2>, p_force);
				p_force->SetCutOffLength(1.5);

				p_force->SetCellCellSpringStiffness(40.0); //Default is 15
				p_force->SetCellGelSpringStiffness(15.0); //Default is 15
				p_force->SetGelGelSpringStiffness(gel_gel_stiffness); //Default is 15
				simulator.AddForce(p_force);

				c_vector<double,2> point = zero_vector<double>(2);
				c_vector<double,2> normal = zero_vector<double>(2);
				//Impose a wall at x > 0
				normal(0) = -1.0;
				MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
				simulator.AddCellPopulationBoundaryCondition(p_bc1);

				//Impose a wall at x < 20
				point(0) = 20.0;
				normal(0) = 1.0;
				MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal));
				simulator.AddCellPopulationBoundaryCondition(p_bc2);

				//Impose a wall at y > 0
				point(0) = 0.0;
				point(1) = 0.0;
				normal(0) = 0.0;
				normal(1) = -1.0;
				MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal));
				simulator.AddCellPopulationBoundaryCondition(p_bc3);

				//Impose a wall at y < 20
				point(1) = 20.0;
				normal(1) = 1.0;
				MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal));
				simulator.AddCellPopulationBoundaryCondition(p_bc4);


				simulator.Solve();
				SimulationTime::Instance()->Destroy();
				SimulationTime::Instance()->SetStartTime(0.0);
			}
		}
	}
};



#endif /* TESTORGANOIDBLOB_HPP_ */


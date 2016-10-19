#ifndef TESTDUNNBOXMODELFORJAMES_HPP_
#define TESTDUNNBOXMODELFORJAMES_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files
#include "AbstractCellBasedTestSuite.hpp" //Needed for cell-based tests: times simulations, generates random numbers and has cell properties
#include "CheckpointArchiveTypes.hpp" //Needed if we use GetIdentifier() method (which we do)
#include "SmartPointers.hpp" //Enables macros to save typing
#include "CellsGenerator.hpp" //Generates cell population
#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "CylindricalHoneycombMeshGenerator.hpp" //Generates cylindrical mesh
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "VoronoiDataWriter.hpp" //Allows us to visualise output in Paraview
#include "DifferentiatedCellProliferativeType.hpp" //Stops cells from proliferating
#include "StochasticDurationCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "PanethCellMutationState.hpp"
#include "CellLabel.hpp"
#include "FixedRegionPlaneBoundaryCondition.hpp"
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "OrganoidLinearSpringForce.hpp" //Spring force law with variable stiffness constant
#include "OrganoidLinearSpringForceAccountingForPanethCellSize.hpp"
#include "OrganoidBasementMembraneForce.hpp" //Prototype BM force for ring based off SJD's code.
#include "OrganoidAnoikisCellKiller.hpp" //Cell killer to remove proliferative cells that have detached from ring.
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

class TestDunnBoxModelForJames : public AbstractCellBasedTestSuite
{
public:
	void TestTissueSlab() throw(Exception)
	{

		//Get the index parameter (can change later if want to run multiple simulations)
		unsigned index = 1;

		//To be extra careful, we reseed the random number generator
		RandomNumberGenerator::Instance()->Reseed(100*index);

		double dt = 0.01; //Set dt
		double end_time = 10.0; //Set end time
		double sampling_timestep = 0.5/dt;

		//Set all the spring stiffness variables
		double cell_cell_stiffness = 15.0;
		double cell_gel_stiffness = 15.0;
		double gel_gel_stiffness = 15.0;
		double paneth_stiffness_multiplier = 1.0;

		//Set the basement membrane force parameters
		double bm_stiffness = 10.0;
		double target_curvature = -0.3;

		//Set the left and right boundaries for the non-zero target curvature region
		double left_crypt_boundary = 7.5;
		double right_crypt_boundary = 17.5;

		//Set the number of cells across and down for the array
		unsigned cells_across = 20;
		unsigned cells_up = 6;
		unsigned ghosts = 2; //Set the number of ghost node layers

		//Generate the mesh
		CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, ghosts);
		Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

		//Get the real indices
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        //Create the vector of cells
		//Create shared pointers for cell and mutation states
		boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_wildtype_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();

		//Create tissue of cells. Initially we set them all to be differentiated
		std::vector<CellPtr> cells; //Create vector of cells

		for (unsigned i = 0; i<location_indices.size(); i++)
		{
			//Set stochastic duration based cell cycle
			StochasticDurationCellCycleModel* p_cycle_model = new StochasticDurationCellCycleModel();
			p_cycle_model->SetDimension(2);
			double birth_time = 13.0*RandomNumberGenerator::Instance()->ranf(); //We would like the birth time to be ~U(0,13) and set in the past
			p_cycle_model->SetBirthTime(-birth_time);

			CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
			p_cell->SetCellProliferativeType(p_diff_type); //Make cell differentiated
			p_cell->InitialiseCellCycleModel();

			cells.push_back(p_cell);
		}

		//Create cell population
		MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        //Get the maximum height of the real nodes to define the monolayer
        double max_height = 0.0;

        for (unsigned i = 0; i < location_indices.size(); i++)
        {
			unsigned node_index = location_indices[i];
			double y = p_mesh->GetNode(node_index)->rGetLocation()[1];

			if (y > max_height)
			{
				max_height = y;
			}
        }

        //Any cell with the maximum height is converted to be a proliferative cell

		//Create ring of stem cells
		for (unsigned i = 0; i < location_indices.size(); i++)
		{
			unsigned cell_index = location_indices[i];
			CellPtr cell_iter = cell_population.GetCellUsingLocationIndex(cell_index);
			double y = cell_population.GetLocationOfCellCentre(cell_iter)[1];

			if(y==max_height)
			{
				cell_iter->SetCellProliferativeType(p_stem_type);
			}
		}

		OffLatticeSimulation<2> simulator(cell_population);

		//Set output directory
		std::stringstream out;
		out << "BM_" << bm_stiffness << "_CURV_" << target_curvature << "/" << 100*index;
		std::string output_directory = "DunnBoxModel/" + out.str();
		simulator.SetOutputDirectory(output_directory);

		simulator.SetDt(dt);
		simulator.SetSamplingTimestepMultiple(sampling_timestep); //Sample the simulation at every hour
		simulator.SetEndTime(end_time); //Hopefully this is long enough for a steady state

		//Add linear spring force (modified to have three different spring stiffnesses, depending on the type of pair)
		MAKE_PTR(OrganoidLinearSpringForce<2>, p_spring_force);
		p_spring_force->SetCutOffLength(1.5);
		p_spring_force->SetCellCellSpringStiffness(cell_cell_stiffness); //Default is 15
		p_spring_force->SetCellGelSpringStiffness(cell_gel_stiffness); //Default is 15
		p_spring_force->SetGelGelSpringStiffness(gel_gel_stiffness); //Default is 15
		p_spring_force->SetPanethCellStiffnessMultiplier(paneth_stiffness_multiplier); //Default is 1
		simulator.AddForce(p_spring_force);

		//Add basement membrane force
		MAKE_PTR(OrganoidBasementMembraneForce, p_bm_force);
		p_bm_force->SetBasementMembraneParameter(bm_stiffness); //Equivalent to beta in SJD's papers
		p_bm_force->SetRingCurvature(target_curvature); //This is equivalent to 1/R in SJD's papers
		p_bm_force->SetLeftCryptBoundary(left_crypt_boundary); //Set the left crypt boundary of non-zero target curvature
		p_bm_force->SetRightCryptBoundary(right_crypt_boundary); //Set the right crypt boundary of non-zero target curvature
		p_bm_force->SetPositionDependentMultiplier(false);	// Multiplier for basement membrane force that may one day be needed
		simulator.AddForce(p_bm_force);

		//Add anoikis-based cell killer
		MAKE_PTR_ARGS(OrganoidAnoikisCellKiller, p_anoikis_killer, (&cell_population));
		simulator.AddCellKiller(p_anoikis_killer);

		//Fix the bottom row of cells
		c_vector<double, 2> point, normal;

		point(0) = 0.0;
		point(1) = 0.5;
		normal(0) = 0.0;
		normal(1) = -1.0;
		MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal));
		simulator.AddCellPopulationBoundaryCondition(p_bc3);

		simulator.Solve();
	}
};

#endif /* TESTDUNNBOXMODELFORJAMES_HPP_ */

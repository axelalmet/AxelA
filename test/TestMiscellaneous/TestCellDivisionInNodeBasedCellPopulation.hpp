#ifndef TESTCELLDIVISIONINNODEBASEDCELLPOPULATION_HPP_
#define TESTCELLDIVISIONINNODEBASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files
#include "AbstractCellBasedTestSuite.hpp" //Needed for cell-based tests: times simulations, generates random numbers and has cell properties
#include "CheckpointArchiveTypes.hpp" //Needed if we use GetIdentifier() method (which we do)
#include "SmartPointers.hpp" //Enables macros to save typing
#include "PetscSetupAndFinalize.hpp"
#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "CylindricalHoneycombMeshGenerator.hpp" //Generates mesh with periodic vertical BCs
#include "Cylindrical2dNodesOnlyMesh.hpp"
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "NodeBasedCellPopulation.hpp" //Node-based CellPopulation class
#include "SmartPointers.hpp" //Enables macros to save typing
#include "OrganoidLinearSpringForce.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "VoronoiDataWriter.hpp" //Allows us to visualise output in Paraview
#include "PanethCellMutationState.hpp"
#include "FakePetscSetup.hpp" //Forbids tests running in parallel (for now)
#include "GeneralisedLinearSpringForce.hpp" //Linear spring force model
#include "OrganoidLinearSpringForce.hpp" //Linear spring force with variable spring constant definition for different cell types
#include "RepulsionForce.hpp" //Repulsion force for node-based models
#include "NodeBasedBasementMembraneForceWithNearestNeighbours.hpp" //Basement membrane force that defines curvature using two closest epithelial neighbours
#include "OrganoidAnoikisCellKiller.hpp" //Cell killer to remove cells if they lose contact with BM
#include "FixedRegionPlaneBoundaryCondition.hpp"
//#include "CryptBoxBoundaryCondition.hpp" //Box BC that pins cells on the bottom and applies plane BC on the left and right
#include "SloughingSidesCellKiller.hpp" //Removes cells when they go past the vertical edges
//#include "VolumeBasedApoptosisCellKiller.hpp" //Removes cells if they are sufficiently compressed

#include "Debug.hpp"

class TestCellDivisionInNodeBasedCellPopulation : public AbstractCellBasedTestSuite
{
public:
	void TestCellDivisionInLine() throw(Exception)
	{
		//Set simulation time parameters
		double dt = 0.01; //Set dt
		double end_time = 15.0; //Set end time
		double sampling_timestep = 0.1/dt; //Set sampling time step

		//Set radius of connectivity
		double radius_of_interaction = 1.5;

		//Set division separation
		double division_separation = 1.0;

		//Set spring stiffness constants
		double cell_cell_stiffness = 15.0;
		double cell_stroma_stiffness = 15.0;
		double stroma_stroma_stiffness = 15.0;

		//Set the age of maturity for spring growth after division
		double age_of_maturity = 10.0;

		//Command line argument stuff: Get the index parameter
		//			unsigned index = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-myintval");
		unsigned index = 10;

		//To be extra careful, we reseed the random number generator
		RandomNumberGenerator::Instance()->Reseed(100*index);

//		HoneycombMeshGenerator generator(cells_across, cells_down, 0); //Create mesh
//		MutableMesh<2,2>* p_generating_mesh = generator.GetMesh(); //Generate mesh

		//Create some nodes
		std::vector<Node<2>*> nodes;

		//Create some nodes
		nodes.push_back(new Node<2>(0u, false, 0.0, 1.0));
		nodes.push_back(new Node<2>(1u, false, 1.0, 1.0));
		nodes.push_back(new Node<2>(2u, false, 2.0, 1.0));
		nodes.push_back(new Node<2>(3u, false, 3.0, 1.0));
		nodes.push_back(new Node<2>(4u, false, 4.0, 1.0));
		nodes.push_back(new Node<2>(5u, false, 5.0, 1.0));
		nodes.push_back(new Node<2>(6u, false, 6.0, 1.0));


		//Construct the nodes only mesh using the MutableMesh
        NodesOnlyMesh<2> mesh;
		mesh.ConstructNodesWithoutMesh(nodes, radius_of_interaction);

		//Create shared pointers for cell and mutation states
		boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_wildtype_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();

		//Create tissue of cells. Initially we set them all to be differentiated
		std::vector<CellPtr> cells; //Create vector of cells

		for (unsigned i = 0; i<mesh.GetNumNodes(); i++)
		{
			//Set stochastic duration based cell cycle
			StochasticDurationCellCycleModel* p_cycle_model = new StochasticDurationCellCycleModel();
			p_cycle_model->SetDimension(2);
			double birth_time = 8.0*RandomNumberGenerator::Instance()->ranf(); //We would like the birth time to be ~U(0,12) and set in the past
			p_cycle_model->SetBirthTime(-birth_time);

			if (i != 3) //Only allow the middle to proliferate
			{
				p_cycle_model->SetTransitCellG1Duration(DBL_MAX);
			}

			CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
			p_cell->SetCellProliferativeType(p_stem_type); //Make cell differentiated
			p_cell->InitialiseCellCycleModel();

			cells.push_back(p_cell);
		}

		NodeBasedCellPopulation<2> cell_population(mesh, cells); //Create the cell population

		//Reset the division separation length
		cell_population.SetMeinekeDivisionSeparation(division_separation);

		//Explicitly request output for Paraview
		cell_population.AddPopulationWriter<VoronoiDataWriter>();

		OffLatticeSimulation<2> simulator(cell_population); //Define the simulation class

		//Simulation pre-amble
		std::stringstream run;
		run << 100*index;
		std::string output_directory = "NodeBasedCellDivisionInLine/" + run.str() + "/";
		simulator.SetOutputDirectory(output_directory);

		simulator.SetDt(dt); //Set dt
		simulator.SetSamplingTimestepMultiple(sampling_timestep); //Set when to sample the simulation
		simulator.SetEndTime(end_time); //Set the end time


		//Set the generalised linear spring force law
		MAKE_PTR(OrganoidLinearSpringForce<2>, p_variable_constant_force);
		MAKE_PTR(RepulsionForce<2>, p_repulsion_force);
		p_repulsion_force->SetMeinekeSpringStiffness(cell_cell_stiffness);
		p_repulsion_force->SetCutOffLength(radius_of_interaction);
		p_repulsion_force->SetMeinekeSpringGrowthDuration(age_of_maturity);
		p_repulsion_force->SetMeinekeDivisionRestingSpringLength(division_separation);
		simulator.AddForce(p_repulsion_force);

//		MAKE_PTR(OrganoidLinearSpringForce<2>, p_variable_constant_force);
//		p_variable_constant_force->SetCellCellSpringStiffness(cell_cell_stiffness); //Set stiffness between epithelial cells
//		p_variable_constant_force->SetCellGelSpringStiffness(cell_stroma_stiffness); //Set stiffness between epithelial cells and stromal cells
//		p_variable_constant_force->SetGelGelSpringStiffness(stroma_stroma_stiffness); //Set stiffness between stromal cells
//		p_variable_constant_force->SetCutOffLength(radius_of_interaction);
//		p_variable_constant_force->SetMeinekeSpringGrowthDuration(age_of_maturity);
//		p_variable_constant_force->SetMeinekeDivisionRestingSpringLength(division_separation);
//		simulator.AddForce(p_variable_constant_force);

		simulator.Solve(); //Solve the simulation

	}

	void TestCellDivisionInTissueSlab() throw(Exception)
		{
		//Set simulation time parameters
		double dt = 0.01; //Set dt
		double end_time = 30.0; //Set end time
		double sampling_timestep = 0.1/dt; //Set sampling time step

		//Set radius of connectivity
		double radius_of_interaction = 1.5;

		//Set the initial spring separation after division
		double division_separation = 0.01;

		//Set maturity age for cells (to determine spring growth rate)
		double age_of_maturity = 10.0;

		//Set spring stiffness constants
		double cell_cell_stiffness = 15.0;
		double cell_stroma_stiffness = 15.0;
		double stroma_stroma_stiffness = 15.0;

		//Command line argument stuff: Get the index parameter
		//			unsigned index = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-myintval");
		unsigned index = 10;

		//To be extra careful, we reseed the random number generator
		RandomNumberGenerator::Instance()->Reseed(100*index);


		//Generate mesh to extract the nodes
		HoneycombMeshGenerator generator(5,5);
		MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();

		//Construct the nodes only mesh using the MutableMesh
        NodesOnlyMesh<2> mesh;
		mesh.ConstructNodesWithoutMesh(*p_generating_mesh, radius_of_interaction);

		//Create shared pointers for cell and mutation states
		boost::shared_ptr<AbstractCellProperty> p_diff_type = CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_wildtype_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();

		//Create tissue of cells. Initially we set them all to be differentiated
		std::vector<CellPtr> cells; //Create vector of cells

		for (unsigned i = 0; i<mesh.GetNumNodes(); i++)
		{
			//Set stochastic duration based cell cycle
			StochasticDurationCellCycleModel* p_cycle_model = new StochasticDurationCellCycleModel();
			p_cycle_model->SetDimension(2);
			double birth_time = 8.0*RandomNumberGenerator::Instance()->ranf(); //We would like the birth time to be ~U(0,12) and set in the past
			p_cycle_model->SetBirthTime(-birth_time);

			if (i != 12) //Only allow the middle to proliferate
			{
				p_cycle_model->SetTransitCellG1Duration(DBL_MAX);
			}

			CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
			p_cell->SetCellProliferativeType(p_stem_type); //Make cell differentiated
			p_cell->InitialiseCellCycleModel();

			cells.push_back(p_cell);
		}

		NodeBasedCellPopulation<2> cell_population(mesh, cells); //Create the cell population

		//Reset the division separation length
		cell_population.SetMeinekeDivisionSeparation(division_separation);

		//Explicitly request output for Paraview
		cell_population.AddPopulationWriter<VoronoiDataWriter>();

		OffLatticeSimulation<2> simulator(cell_population); //Define the simulation class

		//Simulation pre-amble
		std::stringstream run;
		run << 100*index;
		std::string output_directory = "NodeBasedCellDivisionInBox/" + run.str() + "/";
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
		p_variable_constant_force->SetMeinekeSpringGrowthDuration(age_of_maturity);
		p_variable_constant_force->SetMeinekeDivisionRestingSpringLength(division_separation);
		simulator.AddForce(p_variable_constant_force);

		simulator.Solve(); //Solve the simulation
		}
};

#endif /*TESTCELLDIVISIONINNODEBASEDCELLPOPULATION_HPP_*/

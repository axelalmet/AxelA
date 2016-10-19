#include "NodeBasedBasementMembraneForceWithNearestNeighbours.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "AbstractCellProperty.hpp"
#include "Debug.hpp"
#include "Exception.hpp"

/*
 * MODIFIED FOR PERSONAL USE BY AXEL ALMET
 * Last modified: Feb 22 2016
 */

/**
 * To avoid warnings on some compilers, C++ style initialization of member
 * variables should be done in the order they are defined in the header file.
 */
NodeBasedBasementMembraneForceWithNearestNeighbours::NodeBasedBasementMembraneForceWithNearestNeighbours()
   :  AbstractForce<2>(),
   mBasementMembraneParameter(DOUBLE_UNSET),
   mTargetCurvature(DOUBLE_UNSET),
   mLeftBoundary(DOUBLE_UNSET),
   mRightBoundary(DOUBLE_UNSET),
   mUsePositionDependentMembraneForce(false),
   mApplyForceToCrypt(true),
   mMembraneForceMultiplier(DOUBLE_UNSET),
   mCutOffRadius(1.5)
{
    // Sets up output file
//	OutputFileHandler output_file_handler("CurvatureData/", false);
//	mMeinekeOutputFile = output_file_handler.OpenOutputFile("results.curvature");
}

NodeBasedBasementMembraneForceWithNearestNeighbours::~NodeBasedBasementMembraneForceWithNearestNeighbours()
{
//    mMeinekeOutputFile->close();
}

void NodeBasedBasementMembraneForceWithNearestNeighbours::SetBasementMembraneParameter(double basementMembraneParameter)
{
	mBasementMembraneParameter = basementMembraneParameter;
}

double NodeBasedBasementMembraneForceWithNearestNeighbours::GetBasementMembraneParameter()
{
	return mBasementMembraneParameter;
}


void NodeBasedBasementMembraneForceWithNearestNeighbours::SetTargetCurvature(double targetCurvature)
{
	mTargetCurvature = targetCurvature;
}


double NodeBasedBasementMembraneForceWithNearestNeighbours::GetTargetCurvature()
{
	return mTargetCurvature;
}

void NodeBasedBasementMembraneForceWithNearestNeighbours::SetLeftCryptBoundary(double leftBoundary)
{
	mLeftBoundary = leftBoundary;
}

double NodeBasedBasementMembraneForceWithNearestNeighbours::GetLeftCryptBoundary()
{
	return mLeftBoundary;
}

void NodeBasedBasementMembraneForceWithNearestNeighbours::SetRightCryptBoundary(double rightBoundary)
{
	mRightBoundary = rightBoundary;
}

double NodeBasedBasementMembraneForceWithNearestNeighbours::GetRightCryptBoundary()
{
	return mRightBoundary;
}

void NodeBasedBasementMembraneForceWithNearestNeighbours::SetCryptGeometry(bool applyForceToCrypt)
{
	mApplyForceToCrypt = applyForceToCrypt;
}

bool NodeBasedBasementMembraneForceWithNearestNeighbours::GetCryptGeometryCheck()
{
	return mApplyForceToCrypt;
}

void NodeBasedBasementMembraneForceWithNearestNeighbours::SetPositionDependentMultiplier(bool usePositionDependentMembraneForce, double MembraneForceMultiplier)
{
	mUsePositionDependentMembraneForce = usePositionDependentMembraneForce;
	mMembraneForceMultiplier = MembraneForceMultiplier;
}


double NodeBasedBasementMembraneForceWithNearestNeighbours::GetPositionDependentMultiplier()
{
	return mMembraneForceMultiplier;
}

double NodeBasedBasementMembraneForceWithNearestNeighbours::GetCutOffRadius()
{
	return mCutOffRadius;
}

void NodeBasedBasementMembraneForceWithNearestNeighbours::SetCutOffRadius(double cutOffRadius)
{
	mCutOffRadius = cutOffRadius;
}


void NodeBasedBasementMembraneForceWithNearestNeighbours::RemoveDuplicates1D(std::vector<unsigned>& rVectorWithDuplicates)
{
    std::sort(rVectorWithDuplicates.begin(), rVectorWithDuplicates.end());
    rVectorWithDuplicates.erase(std::unique(rVectorWithDuplicates.begin(), rVectorWithDuplicates.end()), rVectorWithDuplicates.end());
}

/* Method to return the current coordinates of the crypt orifice and
 * crypt base - these can be used to accurately define the region of the
 * crypt base. (This will be the y-coordinate in 2D, or the z coordinate in 3D)
 * [0] - y or z-coordinate of orifice
 * [1] - y or z-coordinate of base
 */
//Don't think I need this (NodeBased).
//c_vector<double,2> NodeBasedBasementMembraneForceWithNearestNeighbours::GetCryptHeightExtremes(AbstractCellPopulation<2>& rCellPopulation)
//{
//    MeshBasedCellPopulation<2>* p_tissue = static_cast<MeshBasedCellPopulation<2>*>(&rCellPopulation);
//
//    // Create a vector to store the y-coordinates of the lowest point of the crypt base and the highest point of the
//    // crypt orifice
//    c_vector<double,2> height_extremes;
//
//    double max_height = 0.0;
//    double min_height = DBL_MAX;
//
//    double current_height_coordinate;
//
//    // We iterate over all cells in the tissue, and deal only with those that are epithelial cells
//    for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
//         cell_iter != rCellPopulation.End();
//         ++cell_iter)
//    {
//    	boost::shared_ptr<AbstractCellProperty> p_state = cell_iter->GetCellProliferativeType();
//
//	   	// Need these to not be labelled cells
//	   	if ( (p_state->IsType<DifferentiatedCellProliferativeType>()==false) )
//	   	{
//	   		Node<2>* p_node = p_tissue->GetNodeCorrespondingToCell(*cell_iter);
//
//   			current_height_coordinate = p_node->rGetLocation()[1];
//
//	    	if (current_height_coordinate > max_height)
//	    	{
//	    		max_height = current_height_coordinate;
//	    	}
//	    	else if (current_height_coordinate < min_height)
//	    	{
//	    		min_height = current_height_coordinate;
//	    	}
//	    }
//    }
//
//    height_extremes[0] = max_height;
//    height_extremes[1] = min_height;
//
//    return height_extremes;
//}

/*
 * Function to find the curvature along three points, using the method previously described by SJD
 * Method has been adjusted to account for periodic meshes etc, i.e. heavy use of GetVectorFromAtoB
 */
double NodeBasedBasementMembraneForceWithNearestNeighbours::FindParametricCurvature(AbstractCellPopulation<2>& rCellPopulation,
															c_vector<double, 2> leftPoint,
															c_vector<double, 2> centrePoint,
															c_vector<double, 2> rightPoint)
{

	//Get the relevant vectors (all possible differences)
	c_vector<double, 2> left_to_centre = rCellPopulation.rGetMesh().GetVectorFromAtoB(leftPoint, centrePoint);
	c_vector<double, 2> centre_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(centrePoint, rightPoint);
	c_vector<double, 2> left_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(leftPoint, rightPoint);

	// Firstly find the parametric intervals
	double left_s = pow(pow(left_to_centre[0],2) + pow(left_to_centre[1],2), 0.5);
	double right_s = pow(pow(centre_to_right[0],2) + pow(centre_to_right[1],2), 0.5);
//	PRINT_2_VARIABLES(left_s, right_s);

	double sum_intervals = left_s + right_s;

	//Calculate finite difference of first derivatives
	double x_prime = (left_to_right[0])/sum_intervals;
	double y_prime = (left_to_right[1])/sum_intervals;

	//Calculate finite difference of second derivatives
	double x_double_prime = 2*(left_s*centre_to_right[0] - right_s*left_to_centre[0])/(left_s*right_s*sum_intervals);
	double y_double_prime = 2*(left_s*centre_to_right[1] - right_s*left_to_centre[1])/(left_s*right_s*sum_intervals);

	//Calculate curvature using formula
	double curvature = (x_prime*y_double_prime - y_prime*x_double_prime)/pow((pow(x_prime,2) + pow(y_prime,2)), 1.5);

	return curvature;
}

/*
 * Method to find all the epithelial cells that make up the monolayer
 */
std::vector<unsigned> NodeBasedBasementMembraneForceWithNearestNeighbours::GetEpithelialIndices(AbstractCellPopulation<2>& rCellPopulation)
{
	//Create the vector of epithelial cell indices
	std::vector<unsigned> epithelial_indices;

	if (dynamic_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation))
	{
	    NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation);

		for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
				cell_iter != rCellPopulation.End();
				++cell_iter)
		{
			//Get the cell type
			boost::shared_ptr<AbstractCellProperty> p_type = cell_iter->GetCellProliferativeType();

			if(!p_type->IsType<DifferentiatedCellProliferativeType>()) //If we have an epithelial cell
			{
				unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
				epithelial_indices.push_back(node_index);
			}
		}
	}

	return epithelial_indices;
}

/*
 * Method to get the y-coordinates of the crypt base and orifice. This determines where we apply
 * a non-zero target curvature force. First coordinate of the vector will be the y-coordinate of
 * the crypt base, the second will be the y-coordinate of the crypt orifice.
 */
//c_vector<double, 2> NodeBasedBasementMembraneForceWithNearestNeighbours::GetCryptHeightExtremes(AbstractCellPopulation<2>& rCellPopulation)
//{
//	//Initialise vector pair
//	c_vector<double, 2> height_extremes;
//
//	//Get all the epithelial indices (do not have to consider stromal cells)
//	std::vector<unsigned> epithelial_indices = GetEpithelialIndices(rCellPopulation);
//
//	//Initialise the minimum and maximum y-coordinates
//	double min_height = DBL_MAX;
//	double max_height = 0.0;
//
//	//Loop over each index
//	for (unsigned i = 0; i < epithelial_indices.size(); i++)
//	{
//		//Get the epithelial index
//		unsigned current_index = epithelial_indices[i];
//
//		//Get the y-coordinate of the cell
//		double current_height = rCellPopulation.GetNode(current_index)->rGetLocation()[1];
//
//		//If the cell is above the min
//		if (current_height < min_height)
//		{
//			min_height = current_height;
//		}
//		else if (current_height > max_height)
//		{
//			max_height = current_height;
//		}
//	}
//
//	//Add the heights
//	height_extremes[0] = min_height; height_extremes[1] = max_height;
//
//	return height_extremes;
//}

/*
 * Method to get the neighbouring nodes that are epithelial cells
 */
std::vector<unsigned> NodeBasedBasementMembraneForceWithNearestNeighbours::GetNeighbouringEpithelialIndices(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex)
{
	// Create a set of neighbouring node indices
	std::vector<unsigned> neighbouring_epithelial_indices;

	if (dynamic_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation))
	{
		NodeBasedCellPopulation<2>* p_tissue = static_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation);

		//Get cut off radius for defining neighbourhood
		double radius = GetCutOffRadius();

		// Find the indices of the elements owned by this node
		std::set<unsigned> neighbouring_indices = p_tissue->GetNodesWithinNeighbourhoodRadius(nodeIndex, radius);

		// Iterate over these elements
		for (std::set<unsigned>::iterator elem_iter = neighbouring_indices.begin();
				elem_iter != neighbouring_indices.end();
				++elem_iter)
		{
			//Get the cell according to the index
			CellPtr cell_iter = rCellPopulation.GetCellUsingLocationIndex(*elem_iter);

			//Get the cell type
			boost::shared_ptr<AbstractCellProperty> p_type = cell_iter->GetCellProliferativeType();

			//if the cell is not differentiated and thus an epithelial cell, we add it to the vector
			if(!p_type->IsType<DifferentiatedCellProliferativeType>())
			{
				neighbouring_epithelial_indices.push_back(*elem_iter);
			}

		}
	}

    return neighbouring_epithelial_indices;
}

/*
 * Method to get the set of left and right neighbours.
 * We need to use this method to avoid pathological cases when a neighbour can technically be
 * both "left" and "right".
 */
std::pair<std::vector<unsigned>, std::vector<unsigned> > NodeBasedBasementMembraneForceWithNearestNeighbours::GetLeftAndRightNeighbours(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex)
{
	//Initialise pair
	std::pair<std::vector<unsigned>, std::vector<unsigned> > neighbours;

	//Initialise the vectors
	std::vector<std::pair<double, unsigned> > left_neighbours_and_distances, right_neighbours_and_distances;
	std::vector<unsigned> left_neighbours, right_neighbours;

	//Get the set of neighbours
	std::vector<unsigned> neighbouring_epithelial_indices = GetNeighbouringEpithelialIndices(rCellPopulation, nodeIndex);

	//Get the considered node's location
	c_vector<double, 2> epithelial_location = rCellPopulation.GetNode(nodeIndex)->rGetLocation();

	//Get the centre of mass
	c_vector<double, 2> centre_of_mass = GetCentreOfMassOfNeighbours(rCellPopulation, nodeIndex);

	//Get the average line of orientation of the neighbours
	c_vector<double, 2> average_line_of_orientation = GetAveragePlaneOfOrientationOfNeighbours(rCellPopulation, nodeIndex);

	/* Sort the neighbours into a "left" and "right" set based on the sign of the inner product
	* of the average line and the relative vector between the neighbour location and the
	* considered epithelial node
	*/

	//Iterate through the neighbours
	for (unsigned i = 0; i < neighbouring_epithelial_indices.size(); i++)
	{
		//Get the neighbour index
		unsigned neighbour_index = neighbouring_epithelial_indices[i];

		//Get the neighbour location
		c_vector<double, 2> neighbour_location = rCellPopulation.GetNode(neighbour_index)->rGetLocation();

		//Get the relative vector from the centre of mass to the neighbor
		c_vector<double, 2> com_to_neighbour = rCellPopulation.rGetMesh().GetVectorFromAtoB(centre_of_mass, neighbour_location);

		//Get the distance from the epithelial node to the neighbour
		c_vector<double, 2> epithelial_node_to_neighbour = rCellPopulation.rGetMesh().GetVectorFromAtoB(epithelial_location, neighbour_location);
		double distance_to_neighbour = norm_2(epithelial_node_to_neighbour);

		//Calculate the inner product that represents the 'direction' of the neighbour
		double signed_direction = inner_prod(com_to_neighbour, average_line_of_orientation);

		std::pair<double, unsigned> neighbour_index_and_distance = std::make_pair(distance_to_neighbour, neighbour_index);
		if (signed_direction < 0.0) //Negative corresponds to left
		{
			left_neighbours_and_distances.push_back(neighbour_index_and_distance);
		}
		else if (signed_direction > 0.0)  //Positive for right
		{
			right_neighbours_and_distances.push_back(neighbour_index_and_distance);
		}
		else if (signed_direction == 0.0)//If the neighbour is perpendicular (it could happen!), we put it in randomly
		{
			double rand = RandomNumberGenerator::Instance()->ranf();
			if (rand < 0.5)
			{
				left_neighbours_and_distances.push_back(neighbour_index_and_distance);
			}
			else
			{
				right_neighbours_and_distances.push_back(neighbour_index_and_distance);
			}
		}
	}

	//Sort the left and right neighbours by their signed distances in descending order
	std::sort(left_neighbours_and_distances.begin(), left_neighbours_and_distances.end());
	std::sort(right_neighbours_and_distances.begin(), right_neighbours_and_distances.end());

	//Get the vectors of indices for the left and right neighbours

	for (unsigned i = 0; i < left_neighbours_and_distances.size(); i++)
	{
		left_neighbours.push_back(left_neighbours_and_distances[i].second);
	}

	for (unsigned i = 0; i < right_neighbours_and_distances.size(); i++)
	{
		right_neighbours.push_back(right_neighbours_and_distances[i].second);
	}

	//Define the pair
	neighbours = std::make_pair(left_neighbours, right_neighbours);

	return neighbours;
}

/*
 * Method to get the epithelial indices that have a left and right neighbour.
 * Returns a map with pairs <index, <left, right>> for convention.
 */
std::map<unsigned, std::pair<unsigned, unsigned> > NodeBasedBasementMembraneForceWithNearestNeighbours::GetEpithelialIndicesAndTheirLeftAndRightEpithelialNeighbours(AbstractCellPopulation<2>& rCellPopulation)
{
	//Declare the map
	std::map<unsigned, std::pair<unsigned, unsigned> > epithelial_indices_and_their_left_and_right_neighbours;

	//Get the epithelial indices
	std::vector<unsigned> epithelial_indices = GetEpithelialIndices(rCellPopulation);

//	for (unsigned i = 0; i < epithelial_indices.size(); i++)
//	{
//		unsigned epithelial_index = epithelial_indices[i];
//
//		//Get the left and right neighbour set. We do it this way to avoid the weird cases where neighbours
//		//could be either left or right. Inefficient, but it works.
//		std::pair<std::vector<unsigned>, std::vector<unsigned> > neighbours = GetLeftAndRightNeighbours(rCellPopulation, epithelial_index);
//
//		std::vector<unsigned> left_neighbours = neighbours.first;
//		std::vector<unsigned> right_neighbours = neighbours.second;
//
//		//If we have a left and right neighbour
//		if ( (left_neighbours.size() > 0)&&(right_neighbours.size() > 0) )
//		{
//
//			//Make the pair of neighbours
//			std::pair<unsigned, unsigned> left_and_right_neighbour = std::make_pair(left_neighbours[0], right_neighbours[0]);
//
//			//Assign the pair as the value of the node index key
//			epithelial_nodes_and_their_left_and_right_neighbours[epithelial_index] = left_and_right_neighbour;
//		}
//	}

	/* We are going to define the basement membrane using the x-coordinate of epithelial cells.
	 * Doesn't generalise to 3D, but it will give us a confluent layer.This method only applies if
	 * we are modelling a crypt. We have to change how we define neighbours if we are modelling an organoid
	 */

	bool is_force_applied_to_crypt = GetCryptGeometryCheck();
	if (is_force_applied_to_crypt)
	{
		std::vector<std::pair<double, unsigned> > epithelial_indices_and_x_coordinates; //Define vector of pairs, so that we may sort by the x-coordinate

		for (unsigned i = 0; i < epithelial_indices.size(); i++)
		{
			unsigned epithelial_index = epithelial_indices[i]; //Get node index
			double epithelial_x_coordinate = rCellPopulation.GetNode(epithelial_index)->rGetLocation()[0]; //Get x-coordinate of node location

			//Make pair
			std::pair<double, unsigned> x_coordinate_and_index = std::make_pair(epithelial_x_coordinate, epithelial_index);

			epithelial_indices_and_x_coordinates.push_back(x_coordinate_and_index);
		}

		//Sort indices by the x-coodrinate
		std::sort(epithelial_indices_and_x_coordinates.begin(), epithelial_indices_and_x_coordinates.end());

		//We now define the nodes and their left and right neighbours in the layer
		for (unsigned i = 0; i < epithelial_indices_and_x_coordinates.size(); i++)
		{
			unsigned centre_node_index = epithelial_indices_and_x_coordinates[i].second; //Get index of centre node

			//Initialise left and right neighbours
			unsigned left_neighbour_index, right_neighbour_index;

			if (i == 0) //If it is the first index, the 'left' neighbour is the last index due to periodicity
			{
				left_neighbour_index = epithelial_indices_and_x_coordinates[epithelial_indices_and_x_coordinates.size() - 1].second;
				right_neighbour_index = epithelial_indices_and_x_coordinates[1].second;
			}
			else if (i == (epithelial_indices_and_x_coordinates.size() - 1) ) //If it is the last node, the 'right' neighbour is the first index
			{
				left_neighbour_index = epithelial_indices_and_x_coordinates[epithelial_indices_and_x_coordinates.size() - 2].second;
				right_neighbour_index = epithelial_indices_and_x_coordinates[0].second;
			}
			else //Otherwise the left neighbour is the index before and the right is the index after
			{
				left_neighbour_index = epithelial_indices_and_x_coordinates[i - 1].second;
				right_neighbour_index = epithelial_indices_and_x_coordinates[i + 1].second;
			}

			epithelial_indices_and_their_left_and_right_neighbours[centre_node_index] = std::make_pair(left_neighbour_index, right_neighbour_index);
		}
	}
	else //If we are applying the force to an organoid, we need to sort the epithelial cells using polar coordinates (note this assumes tha the initial geomotry is that of an organoid)
	{
		std::vector<std::pair<double, unsigned> > epithelial_indices_and_angles; //Define vector of pairs, so that we may sort by the angles

		//First get centre of mass of epithelial nodes, to define a reference point for calculating the polar angle
		c_vector<double, 2> centre_of_mass = zero_vector<double>(2);

		for (unsigned i = 0; i < epithelial_indices.size(); i++)
		{
			unsigned epithelial_index = epithelial_indices[i]; //Get node index
			c_vector<double, 2> epithelial_location = rCellPopulation.GetNode(epithelial_index)->rGetLocation(); //Get x-coordinate of node location

			centre_of_mass += epithelial_location;
		}

		//Average centre of mass
		centre_of_mass /= epithelial_indices.size();

		for (unsigned k = 0; k < epithelial_indices.size(); k++)
		{
			unsigned cell_index = epithelial_indices[k]; //Get the node index

			CellPtr cell = rCellPopulation.GetCellUsingLocationIndex(cell_index); //Get the cell

			//Get the cell location
			double x = rCellPopulation.GetLocationOfCellCentre(cell)[0];
			double y = rCellPopulation.GetLocationOfCellCentre(cell)[1];

			//Make the pair of the angle and cell
			std::pair<double, unsigned> angle_cell;

			//Get point relative to the circle centre
			double rel_x = x - centre_of_mass[0];
			double rel_y = y - centre_of_mass[1];

			double circle_angle = atan(rel_y/rel_x); //Get initial angle argument

			if (rel_x<0.0) //If the point is in the second quadrant or third quadrant
			{
				circle_angle += M_PI;
			}
			else if ((rel_x>=0.0)&&(rel_y<0.0)) //Fourth quadrant
			{
				circle_angle += 2*M_PI;
			}

			angle_cell = std::make_pair(circle_angle, cell_index);

			epithelial_indices_and_angles.push_back(angle_cell); //Add the angle and node index
		}

		//Sort indices by the angle
		std::sort(epithelial_indices_and_angles.begin(), epithelial_indices_and_angles.end());

		//We now define the nodes and their left and right neighbours in the layer
		for (unsigned i = 0; i < epithelial_indices_and_angles.size(); i++)
		{
			unsigned centre_node_index = epithelial_indices_and_angles[i].second; //Get index of centre node

			//Initialise left and right neighbours
			unsigned left_neighbour_index, right_neighbour_index;

			if (i == 0) //If it is the first index, the 'left' neighbour is the last index due to periodicity
			{
				left_neighbour_index = epithelial_indices_and_angles[epithelial_indices_and_angles.size() - 1].second;
				right_neighbour_index = epithelial_indices_and_angles[1].second;
			}
			else if (i == (epithelial_indices_and_angles.size() - 1) ) //If it is the last node, the 'right' neighbour is the first index
			{
				left_neighbour_index = epithelial_indices_and_angles[epithelial_indices_and_angles.size() - 2].second;
				right_neighbour_index = epithelial_indices_and_angles[0].second;
			}
			else //Otherwise the left neighbour is the index before and the right is the index after
			{
				left_neighbour_index = epithelial_indices_and_angles[i - 1].second;
				right_neighbour_index = epithelial_indices_and_angles[i + 1].second;
			}

			epithelial_indices_and_their_left_and_right_neighbours[centre_node_index] = std::make_pair(left_neighbour_index, right_neighbour_index);
		}
	}

	return epithelial_indices_and_their_left_and_right_neighbours;
}

/*
 * Method to remove pairs of adjacent eptihelial-epithelial triangles of nodes
 * and their nearest neighbours, if they exist
 */
std::map<unsigned, std::pair<unsigned, unsigned> > NodeBasedBasementMembraneForceWithNearestNeighbours::CheckForAdjacentEpithelialTriangles(AbstractCellPopulation<2>& rCellPopulation, std::map<unsigned, std::pair<unsigned, unsigned> > epithelialIndicesAndNeighbours)
{
	//Get list of adjacent triangles
	std::vector<std::pair<c_vector<unsigned, 3>, c_vector<unsigned, 3> > > adjacent_triangle_pairs = GetAdjacentEpithelialTrianglePairs(rCellPopulation, epithelialIndicesAndNeighbours);

	//Loop over each pair
	for (unsigned i =  0; i < adjacent_triangle_pairs.size(); i++)
	{
		std::pair<c_vector<unsigned, 3>, c_vector<unsigned, 3> > epithelial_pair = adjacent_triangle_pairs[i];

	    std::vector<std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> > information_about_quadrilateral = GetRelevantInformationAboutQuadrilateral(rCellPopulation, epithelial_pair);

	    //Preliminary check that we can re-assign edges appropriately
	    unsigned number_of_edges_connected_to_other_triangles = 0;

	    for (unsigned j = 0; j < information_about_quadrilateral.size(); j++)
	    {
	    	if (DoesAnotherNodeHaveTheseNeighbours(epithelialIndicesAndNeighbours, information_about_quadrilateral[j].first))
	    	{
	    		number_of_edges_connected_to_other_triangles += 1;
	    	}
	    }

	    //If there is more than one edge attached to a triangle then the quadrilateral
	    //is attached to two other triangles, preventing us from re-assigning edges appropriately
	    if (number_of_edges_connected_to_other_triangles > 2)
	    {
	    	EXCEPTION("Considered quadrilateral's edges are attached to too many triangles to re-assign the neighbours appropriately.");
	    }
	    else
	    {
		    int iter = 3;
		    //Get the edge with the largest distance that isn't part of a triangle
		    while(DoesAnotherNodeHaveTheseNeighbours(epithelialIndicesAndNeighbours, information_about_quadrilateral[iter].first))
		    {
		    	iter -= 1;
		    }

	    	//Get the edge and the opposite edge to be deleted
	    	std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> node_and_edge_to_be_deleted_A = information_about_quadrilateral[iter];
	    	std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> node_and_edge_to_be_deleted_B = GetOppositeQuadrilateralEdge(information_about_quadrilateral, node_and_edge_to_be_deleted_A);

	    	//Second check that we can re-assign both edges appropriately:
	    	//If one of these edges is attached to a triangle, we have to take the other pair of edges
	    	if ( (DoesAnotherNodeHaveTheseNeighbours(epithelialIndicesAndNeighbours, node_and_edge_to_be_deleted_A.first))||
	    			(DoesAnotherNodeHaveTheseNeighbours(epithelialIndicesAndNeighbours, node_and_edge_to_be_deleted_B.first)) )
	    	{
	    		//Get the first edge that isn't either of the above edges
	    		for (unsigned j = 0; j < information_about_quadrilateral.size(); j++)
	    		{
	    			std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> considered_node_and_neighbours = information_about_quadrilateral[j];

	    			//Get the first edge that doesn't have the same neighbours as the previous two
	    			if( !((considered_node_and_neighbours.first.first==node_and_edge_to_be_deleted_A.first.first)&&(considered_node_and_neighbours.first.second.first==node_and_edge_to_be_deleted_A.first.second.first)
	    					&&(considered_node_and_neighbours.first.second.second==node_and_edge_to_be_deleted_A.first.second.second)) )
	    			{
	    				if 	    					(!((considered_node_and_neighbours.first.first==node_and_edge_to_be_deleted_B.first.first)&&(considered_node_and_neighbours.first.second.first==node_and_edge_to_be_deleted_B.first.second.first)
	    						&&(considered_node_and_neighbours.first.second.second==node_and_edge_to_be_deleted_B.first.second.second)) )
	    				{
	    					//Reassign edge pair
	    					std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> node_and_edge_to_be_deleted_A = considered_node_and_neighbours;
	    					std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> node_and_edge_to_be_deleted_B = GetOppositeQuadrilateralEdge(information_about_quadrilateral, node_and_edge_to_be_deleted_A);
	    					break; //Break out of inner for loop

	    				}

	    			}
	    		}
	    	}

	    	//Check again. If we have the same issue, throw an exception
	    	if ( (DoesAnotherNodeHaveTheseNeighbours(epithelialIndicesAndNeighbours, node_and_edge_to_be_deleted_A.first))||
	    			(DoesAnotherNodeHaveTheseNeighbours(epithelialIndicesAndNeighbours, node_and_edge_to_be_deleted_B.first)) )
	    	{
	    		EXCEPTION("Cannot re-assign edges appropriately.");
	    	}
	    	else
	    	{
	    		//Get node A and its neighbours
	    		unsigned node_A_index = node_and_edge_to_be_deleted_A.first.first;
	    		std::pair<unsigned, unsigned> neighbours_of_node_A = node_and_edge_to_be_deleted_A.first.second;

	    		//Get node B and its neighbours
	    		unsigned node_B_index = node_and_edge_to_be_deleted_B.first.first;
	    		std::pair<unsigned, unsigned> neighbours_of_node_B = node_and_edge_to_be_deleted_B.first.second;

	    		//Re-assign node_A and node_B's neighbours to be the endpoints of the respective edges.
	    		//This ensures the nodes will be 'torqued' in the right direction.
	    		if (neighbours_of_node_A.first == node_B_index)
	    		{
    				//Get the second closest left and right neighbours of A's right neighbour. We do it this way to avoid
    				//the weird case when they may coincide.
    				std::pair<std::vector<unsigned>, std::vector<unsigned> > neighbouring_indices = GetLeftAndRightNeighbours(rCellPopulation, neighbours_of_node_A.second);

    				std::vector<unsigned> left_neighbours = neighbouring_indices.first;
    				std::vector<unsigned> right_neighbours = neighbouring_indices.second;

	    			if((left_neighbours.size() > 1)&&(right_neighbours.size() > 1))
	    			{

	    				unsigned second_nearest_left_neighbour = left_neighbours[1];
	    				unsigned second_nearest_right_neighbour = right_neighbours[1];

	    				//We now check neither of these neighbours are in the original quadrilateral, for stability reasons

	    				//Create a vector of the new indices, so that we can use the 'find' method
	    				std::vector<unsigned> points_in_quadrilateral;
	    				points_in_quadrilateral.push_back(node_A_index);
	    				points_in_quadrilateral.push_back(node_B_index);
	    				points_in_quadrilateral.push_back(neighbours_of_node_A.first);
	    				points_in_quadrilateral.push_back(neighbours_of_node_A.second);
	    				points_in_quadrilateral.push_back(neighbours_of_node_B.first);
	    				points_in_quadrilateral.push_back(neighbours_of_node_B.second);

	    				//If neither the third nor fourth neighbour are in this quadrilateral, we can assign them as neighbours
	    				if( (std::find(points_in_quadrilateral.begin(), points_in_quadrilateral.end(), second_nearest_left_neighbour)==points_in_quadrilateral.end())
	    						&&(std::find(points_in_quadrilateral.begin(), points_in_quadrilateral.end(), second_nearest_right_neighbour)==points_in_quadrilateral.end()))
	    				{
	    					epithelialIndicesAndNeighbours[neighbours_of_node_A.second] = std::make_pair(second_nearest_left_neighbour, second_nearest_right_neighbour);
	    				}
	    				else //Else delete the node from the map
	    				{
	    					if (epithelialIndicesAndNeighbours.find(neighbours_of_node_A.second) != epithelialIndicesAndNeighbours.end())
	    					{
	    						epithelialIndicesAndNeighbours.erase(epithelialIndicesAndNeighbours.find(neighbours_of_node_A.second));
	    					}
	    				}
	    			}
	    			else
	    			{
    					if (epithelialIndicesAndNeighbours.find(neighbours_of_node_A.second) != epithelialIndicesAndNeighbours.end())
    					{
    						epithelialIndicesAndNeighbours.erase(epithelialIndicesAndNeighbours.find(neighbours_of_node_A.second));
    					}
	    			}

	    			//Re-assign the neighbours of A in terms of 'left' and 'right'-ness.

	    			//Get the location of node A and the other neighbour. We check with this neighbour because
	    			//node B may be perpendicular to the average plane (less likely with the other neighbour as it's
	    			//on the acute corner)
	    			c_vector<double, 2> node_A_location = rCellPopulation.GetNode(node_A_index)->rGetLocation();
	    			c_vector<double, 2> second_neighbour_of_A_location = rCellPopulation.GetNode(neighbours_of_node_A.second)->rGetLocation();

	    			c_vector<double, 2> A_to_second_neighbour = rCellPopulation.rGetMesh().GetVectorFromAtoB(node_A_location, second_neighbour_of_A_location);

	    			//Get the average plane of orientation of node A
	    			c_vector<double, 2> average_plane_of_orientation_A = GetAveragePlaneOfOrientationOfNeighbours(rCellPopulation, node_A_index);

	    			//Get the inner product
	    			double signed_direction = inner_prod(average_plane_of_orientation_A, A_to_second_neighbour);

	    			if (signed_direction < 0.0) //If B is the right neighbour of A
	    			{
	    				epithelialIndicesAndNeighbours[node_A_index] = neighbours_of_node_A;
	    			}
	    			else
	    			{
	    				epithelialIndicesAndNeighbours[node_A_index] = std::make_pair(node_B_index, neighbours_of_node_A.second);
	    			}
	    		}
	    		else
	    		{
    				//Get the second closest left and right neighbours of A's left neighbour. We do it this way to avoid
    				//the weird case when they may coincide.
    				std::pair<std::vector<unsigned>, std::vector<unsigned> > neighbouring_indices = GetLeftAndRightNeighbours(rCellPopulation, neighbours_of_node_A.first);

    				std::vector<unsigned> left_neighbours = neighbouring_indices.first;
    				std::vector<unsigned> right_neighbours = neighbouring_indices.second;

	    			if( (left_neighbours.size() > 1)&&(right_neighbours.size() > 1))
	    			{


	    				unsigned second_nearest_left_neighbour = left_neighbours[1];
	    				unsigned second_nearest_right_neighbour = right_neighbours[1];

	    				//We now check neither of these neighbours are in the original quadrilateral, for stability reasons

	    				//Create a vector of the new indices, so that we can use the 'find' method
	    				std::vector<unsigned> points_in_quadrilateral;
	    				points_in_quadrilateral.push_back(node_A_index);
	    				points_in_quadrilateral.push_back(node_B_index);
	    				points_in_quadrilateral.push_back(neighbours_of_node_A.first);
	    				points_in_quadrilateral.push_back(neighbours_of_node_A.second);
	    				points_in_quadrilateral.push_back(neighbours_of_node_B.first);
	    				points_in_quadrilateral.push_back(neighbours_of_node_B.second);

	    				//If neither the third nor fourth neighbour are in this quadrilateral, we can assign them as neighbours
	    				if( (std::find(points_in_quadrilateral.begin(), points_in_quadrilateral.end(), second_nearest_left_neighbour)==points_in_quadrilateral.end())
	    						&&(std::find(points_in_quadrilateral.begin(), points_in_quadrilateral.end(), second_nearest_right_neighbour)==points_in_quadrilateral.end()))
	    				{
	    					epithelialIndicesAndNeighbours[neighbours_of_node_A.first] = std::make_pair(second_nearest_left_neighbour, second_nearest_right_neighbour);
	    				}
	    				else //Else delete the node from the map
	    				{
	    					if (epithelialIndicesAndNeighbours.find(neighbours_of_node_A.first) != epithelialIndicesAndNeighbours.end())
	    					{
	    						epithelialIndicesAndNeighbours.erase(epithelialIndicesAndNeighbours.find(neighbours_of_node_A.first));
	    					}
	    				}
	    			}
	    			else
	    			{
    					if (epithelialIndicesAndNeighbours.find(neighbours_of_node_B.first) != epithelialIndicesAndNeighbours.end())
    					{
    						epithelialIndicesAndNeighbours.erase(epithelialIndicesAndNeighbours.find(neighbours_of_node_A.first));
    					}
	    			}

	    			//Re-assign the neighbours of A in terms of 'left' and 'right'-ness.

	    			//Get the location of node A and the other neighbour. We check with this neighbour because
	    			//node B may be perpendicular to the average plane (less likely with the other neighbour as it's
	    			//on the acute corner)
	    			c_vector<double, 2> node_A_location = rCellPopulation.GetNode(node_A_index)->rGetLocation();
	    			c_vector<double, 2> first_neighbour_of_A_location = rCellPopulation.GetNode(neighbours_of_node_A.first)->rGetLocation();

	    			c_vector<double, 2> A_to_first_neighbour = rCellPopulation.rGetMesh().GetVectorFromAtoB(node_A_location, first_neighbour_of_A_location);

	    			//Get the average plane of orientation of node A
	    			c_vector<double, 2> average_plane_of_orientation_A = GetAveragePlaneOfOrientationOfNeighbours(rCellPopulation, node_A_index);

	    			//Get the inner product
	    			double signed_direction = inner_prod(average_plane_of_orientation_A, A_to_first_neighbour);

	    			if (signed_direction < 0.0) //If B is the right neighbour of A
	    			{
	    				epithelialIndicesAndNeighbours[node_A_index] = neighbours_of_node_A;
	    			}
	    			else
	    			{
	    				epithelialIndicesAndNeighbours[node_A_index] = std::make_pair(node_B_index, neighbours_of_node_A.first);
	    			}

	    		}

	    		//We now do the same for node B's neighbours
	    		if (neighbours_of_node_B.first == node_A_index)
	    		{
    				//Get the second closest left and right neighbours of B's right neighbour. We do it this way to avoid
    				//the weird case when they may coincide.
    				std::pair<std::vector<unsigned>, std::vector<unsigned> > neighbouring_indices = GetLeftAndRightNeighbours(rCellPopulation, neighbours_of_node_B.second);

    				std::vector<unsigned> left_neighbours = neighbouring_indices.first;
    				std::vector<unsigned> right_neighbours = neighbouring_indices.second;

	    			if( (left_neighbours.size() > 1)&&(right_neighbours.size() > 1) )
	    			{
	    				unsigned second_nearest_left_neighbour = left_neighbours[1];
	    				unsigned second_nearest_right_neighbour = right_neighbours[1];

	    				//We now check neither of these neighbours are in the original quadrilateral, for stability reasons

	    				//Create a vector of the new indices, so that we can use the 'find' method
	    				std::vector<unsigned> points_in_quadrilateral;
	    				points_in_quadrilateral.push_back(node_A_index);
	    				points_in_quadrilateral.push_back(node_B_index);
	    				points_in_quadrilateral.push_back(neighbours_of_node_A.first);
	    				points_in_quadrilateral.push_back(neighbours_of_node_A.second);
	    				points_in_quadrilateral.push_back(neighbours_of_node_B.first);
	    				points_in_quadrilateral.push_back(neighbours_of_node_B.second);

	    				//If neither the third nor fourth neighbour are in this quadrilateral, we can assign them as neighbours
	    				if( (std::find(points_in_quadrilateral.begin(), points_in_quadrilateral.end(), second_nearest_left_neighbour)==points_in_quadrilateral.end())
	    						&&(std::find(points_in_quadrilateral.begin(), points_in_quadrilateral.end(), second_nearest_right_neighbour)==points_in_quadrilateral.end()))
	    				{
	    					epithelialIndicesAndNeighbours[neighbours_of_node_B.second] = std::make_pair(second_nearest_left_neighbour, second_nearest_right_neighbour);
	    				}
	    				else //Else delete the node from the map
	    				{
	    					if (epithelialIndicesAndNeighbours.find(neighbours_of_node_B.second) != epithelialIndicesAndNeighbours.end())
	    					{
	    						epithelialIndicesAndNeighbours.erase(epithelialIndicesAndNeighbours.find(neighbours_of_node_B.second));
	    					}
	    				}
	    			}
	    			else
	    			{
    					if (epithelialIndicesAndNeighbours.find(neighbours_of_node_B.second) != epithelialIndicesAndNeighbours.end())
    					{
    						epithelialIndicesAndNeighbours.erase(epithelialIndicesAndNeighbours.find(neighbours_of_node_B.second));
    					}
	    			}

	    			//Re-assign the neighbours of B in terms of 'left' and 'right'-ness.

	    			//Get the location of node B and the other neighbour. We check with this neighbour because
	    			//node A may be perpendicular to the average plane (less likely with the other neighbour as it's
	    			//on the acute corner)
	    			c_vector<double, 2> node_B_location = rCellPopulation.GetNode(node_B_index)->rGetLocation();
	    			c_vector<double, 2> second_neighbour_of_B_location = rCellPopulation.GetNode(neighbours_of_node_B.second)->rGetLocation();

	    			c_vector<double, 2> B_to_second_neighbour = rCellPopulation.rGetMesh().GetVectorFromAtoB(node_B_location, second_neighbour_of_B_location);

	    			//Get the average plane of orientation of node B
	    			c_vector<double, 2> average_plane_of_orientation_B = GetAveragePlaneOfOrientationOfNeighbours(rCellPopulation, node_B_index);

	    			//Get the inner product
	    			double signed_direction = inner_prod(average_plane_of_orientation_B, B_to_second_neighbour);

	    			if (signed_direction < 0.0) //If B is the right neighbour of A
	    			{
	    				epithelialIndicesAndNeighbours[node_B_index] = neighbours_of_node_B;
	    			}
	    			else
	    			{
	    				epithelialIndicesAndNeighbours[node_B_index] = std::make_pair(node_A_index, neighbours_of_node_B.second);
	    			}
	    		}
	    		else
	    		{
    				//Get the second closest left and right neighbours of A's right neighbour. We do it this way to avoid
    				//the weird case when they may coincide.
    				std::pair<std::vector<unsigned>, std::vector<unsigned> > neighbouring_indices = GetLeftAndRightNeighbours(rCellPopulation, neighbours_of_node_B.first);

    				std::vector<unsigned> left_neighbours = neighbouring_indices.first;
    				std::vector<unsigned> right_neighbours = neighbouring_indices.second;

	    			if((left_neighbours.size() > 1)&&(right_neighbours.size() > 1))
	    			{
	    				unsigned second_nearest_left_neighbour = left_neighbours[1];
	    				unsigned second_nearest_right_neighbour = right_neighbours[1];

	    				//We now check neither of these neighbours are in the original quadrilateral, for stability reasons

	    				//Create a vector of the new indices, so that we can use the 'find' method
	    				std::vector<unsigned> points_in_quadrilateral;
	    				points_in_quadrilateral.push_back(node_A_index);
	    				points_in_quadrilateral.push_back(node_B_index);
	    				points_in_quadrilateral.push_back(neighbours_of_node_A.first);
	    				points_in_quadrilateral.push_back(neighbours_of_node_A.second);
	    				points_in_quadrilateral.push_back(neighbours_of_node_B.first);
	    				points_in_quadrilateral.push_back(neighbours_of_node_B.second);

	    				//If neither the third nor fourth neighbour are in this quadrilateral, we can assign them as neighbours
	    				if( (std::find(points_in_quadrilateral.begin(), points_in_quadrilateral.end(), second_nearest_left_neighbour)==points_in_quadrilateral.end())
	    						&&(std::find(points_in_quadrilateral.begin(), points_in_quadrilateral.end(), second_nearest_right_neighbour)==points_in_quadrilateral.end()))
	    				{
	    					epithelialIndicesAndNeighbours[neighbours_of_node_B.first] = std::make_pair(second_nearest_left_neighbour, second_nearest_right_neighbour);
	    				}
	    				else //Else delete the node from the map, if it is already in there
	    				{
	    					if (epithelialIndicesAndNeighbours.find(neighbours_of_node_B.first) != epithelialIndicesAndNeighbours.end())
	    					{
	    						epithelialIndicesAndNeighbours.erase(epithelialIndicesAndNeighbours.find(neighbours_of_node_B.first));
	    					}
	    				}
	    			}
	    			else
	    			{
	    				if (epithelialIndicesAndNeighbours.find(neighbours_of_node_B.first) != epithelialIndicesAndNeighbours.end())
	    				{
	    					epithelialIndicesAndNeighbours.erase(epithelialIndicesAndNeighbours.find(neighbours_of_node_B.first));
	    				}
	    			}

	    			//Re-assign the neighbours of B in terms of 'left' and 'right'-ness.

	    			//Get the location of node B and the other neighbour. We check with this neighbour because
	    			//node A may be perpendicular to the average plane (less likely with the other neighbour as it's
	    			//on the acute corner)
	    			c_vector<double, 2> node_B_location = rCellPopulation.GetNode(node_B_index)->rGetLocation();
	    			c_vector<double, 2> first_neighbour_of_B_location = rCellPopulation.GetNode(neighbours_of_node_B.first)->rGetLocation();

	    			c_vector<double, 2> B_to_first_neighbour = rCellPopulation.rGetMesh().GetVectorFromAtoB(node_B_location, first_neighbour_of_B_location);

	    			//Get the average plane of orientation of node B
	    			c_vector<double, 2> average_plane_of_orientation_B = GetAveragePlaneOfOrientationOfNeighbours(rCellPopulation, node_B_index);

	    			//Get the inner product
	    			double signed_direction = inner_prod(average_plane_of_orientation_B, B_to_first_neighbour);

	    			if (signed_direction < 0.0) //If B is the right neighbour of A
	    			{
	    				epithelialIndicesAndNeighbours[node_B_index] = neighbours_of_node_B;
	    			}
	    			else
	    			{
	    				epithelialIndicesAndNeighbours[node_B_index] = std::make_pair(node_A_index, neighbours_of_node_B.first);
	    			}
	    		}
	    	}
	    }
	}

	return epithelialIndicesAndNeighbours;
}

/*
 * Method to remove any 'triangles of doom' of epithelial-epithelial connections for
 * epithelial nodes and their nearest neighbours. Furthermore, we re-assign the nearest neighbours if
 * an edge has to be removed.
 */
std::map<unsigned, std::pair<unsigned, unsigned> > NodeBasedBasementMembraneForceWithNearestNeighbours::CheckForEpithelialTriangles(AbstractCellPopulation<2>& rCellPopulation, std::map<unsigned, std::pair<unsigned, unsigned> > epithelialIndicesAndNeighbours)
{
	//Get the epithelial triangles
	std::vector<c_vector<unsigned, 3> > nodes_in_triangles = GetEpithelialTriangles(epithelialIndicesAndNeighbours);

	//Iterate through the triangles and re-assign the 'delete' the longest edge from neighbour connections,
	//provided the edge is not also in another triangle
	for (unsigned i = 0; i < nodes_in_triangles.size(); i++)
	{
		//Get the nodes in the triangle and their locations
		c_vector<unsigned, 3> triangle = nodes_in_triangles[i];

		std::vector<std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> > information_about_triangle = GetRelevantInformationAboutTriangle(rCellPopulation, triangle);

		//Find the longest triangle edge that isn't in another triangle
		int iter = 2;

		while(iter >= 0)
		{
			if(DoesAnotherNodeHaveTheseNeighbours(epithelialIndicesAndNeighbours, information_about_triangle[iter].first) == false)
			{
				break;
			}

			iter -= 1;
		}

		//If all the edges are in a triangle, flag the error
		if(iter < 0)
		{
			EXCEPTION("Considered triangle is an internal triangle.");
		}
		else
		{
			std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> edge_to_be_deleted = information_about_triangle[iter];

			//Get the two endpoints of the edge
			unsigned node_A_index = edge_to_be_deleted.first.second.first;
			unsigned node_B_index = edge_to_be_deleted.first.second.second;

			//If node A originally had a left and right neighbour, re-assign its neighbours
			if (epithelialIndicesAndNeighbours.find(node_A_index) != epithelialIndicesAndNeighbours.end())
			{
				std::pair<unsigned, unsigned> neighbours_of_node_A = epithelialIndicesAndNeighbours[node_A_index];

				//Get the left and right neighbour sets
				std::pair<std::vector<unsigned>, std::vector<unsigned> > left_and_right_neighbours_of_node_A = GetLeftAndRightNeighbours(rCellPopulation, node_A_index);
				std::vector<unsigned> left_neighbours_of_node_A = left_and_right_neighbours_of_node_A.first;
				std::vector<unsigned> right_neighbours_of_node_A = left_and_right_neighbours_of_node_A.second;

				//Re-assign the neighbours of A
				if (neighbours_of_node_A.first == node_B_index) //if node B is the left neighbour of A
				{
					//Only re-assign if a second left neighbour is available
					if (left_neighbours_of_node_A.size() > 1)
					{
						//Get the second closest left neighbour
						unsigned second_closest_left_neighbour_of_A = left_neighbours_of_node_A[1];;

						//Re-assign the neighbours of
						epithelialIndicesAndNeighbours[node_A_index] = std::make_pair(second_closest_left_neighbour_of_A,
								neighbours_of_node_A.second);
					}
					else //Delete node A from map
					{
						epithelialIndicesAndNeighbours.erase(epithelialIndicesAndNeighbours.find(node_A_index));
					}
				}
				else if (neighbours_of_node_A.second == node_B_index) //Else if node B is right neighbour of A
				{
					//Only assign if a second right neighbour is available
					if (right_neighbours_of_node_A.size() > 1)
					{
						//Get the second closest right neighbour
						unsigned second_closest_right_neighbour_of_A = right_neighbours_of_node_A[1];;

						//Re-assign the neighbours of
						epithelialIndicesAndNeighbours[node_A_index] = std::make_pair(neighbours_of_node_A.first,
								second_closest_right_neighbour_of_A);
					}
					else //Delete node A from map
					{
						epithelialIndicesAndNeighbours.erase(epithelialIndicesAndNeighbours.find(node_A_index));
					}
				}
			}

			//If node B had left and right neighbours, re-assign them. Note that at least one of A or
			//B has to have a complete pair of neighbours, they were found in a triangle.
			if (epithelialIndicesAndNeighbours.find(node_B_index) != epithelialIndicesAndNeighbours.end())
			{
				std::pair<unsigned, unsigned> neighbours_of_node_B = epithelialIndicesAndNeighbours[node_B_index];

				//Get the left and right neighbours of B
				std::pair<std::vector<unsigned>, std::vector<unsigned> > left_and_right_neighbours_of_node_B = GetLeftAndRightNeighbours(rCellPopulation, node_B_index);
				std::vector<unsigned> left_neighbours_of_node_B = left_and_right_neighbours_of_node_B.first;
				std::vector<unsigned> right_neighbours_of_node_B = left_and_right_neighbours_of_node_B.second;

				//Re-assign the neighbours of B
				if (neighbours_of_node_B.first == node_A_index) //if node B is the left neighbour of A
				{
					//Only assign if a second left neighbour is available
					if (left_neighbours_of_node_B.size() > 1)
					{
						//Get the second closest left neighbour
						unsigned second_closest_left_neighbour_of_B = left_neighbours_of_node_B[1];

						//Re-assign the neighbours of
						epithelialIndicesAndNeighbours[node_B_index] = std::make_pair(second_closest_left_neighbour_of_B,
								neighbours_of_node_B.second);
					}
					else //Delete node A from map
					{
						epithelialIndicesAndNeighbours.erase(epithelialIndicesAndNeighbours.find(node_B_index));
					}
				}
				else if (neighbours_of_node_B.second == node_A_index) //Else if node B is right neighbour of A
				{
					//Only assign if a second right neighbour is available
					if (right_neighbours_of_node_B.size() > 1)
					{
						//Get the second closest right neighbour
						unsigned second_closest_right_neighbour_of_B = right_neighbours_of_node_B[1];

						//Re-assign the neighbours of
						epithelialIndicesAndNeighbours[node_B_index] = std::make_pair(neighbours_of_node_B.first,
								second_closest_right_neighbour_of_B);
					}
					else //Delete node B from map
					{
						epithelialIndicesAndNeighbours.erase(epithelialIndicesAndNeighbours.find(node_B_index));
					}
				}
			}
		}
	}

	return epithelialIndicesAndNeighbours;
}

/*
 * Method to obtain list of epithelial triangles, if they exist
 */
std::vector<c_vector<unsigned, 3> > NodeBasedBasementMembraneForceWithNearestNeighbours::GetEpithelialTriangles(std::map<unsigned, std::pair<unsigned, unsigned> > epithelialIndicesAndNeighbours)
{
	//Declare vector of 'triples' of nodes that are in a triangle
	std::vector<c_vector<unsigned, 3> > nodes_in_triangles;

	//Iterate over each node and its neighbours
	for (std::map<unsigned, std::pair<unsigned, unsigned> >::iterator map_iter = epithelialIndicesAndNeighbours.begin();
			map_iter != epithelialIndicesAndNeighbours.end();
			map_iter++)
	{

		//Get the node
		unsigned epithelial_index = map_iter->first;
		//Get its neighbours
		std::pair<unsigned, unsigned> first_pair_of_neighbours = map_iter->second;
		unsigned first_neighbour = first_pair_of_neighbours.first;
		unsigned second_neighbour = first_pair_of_neighbours.second;

		//Check if node is part of a triangle
		bool is_node_in_triangle = IsNodePartOfATriangle(epithelialIndicesAndNeighbours, epithelial_index);

		//If we do indeed have a triangle
		if(is_node_in_triangle == true)
		{
			//Create vector of node indices
			c_vector<unsigned, 3> epithelial_triangle;

			//I don't know if there's an easier way to sort the node indices
			std::vector<unsigned> triangle_of_nodes;
			triangle_of_nodes.push_back(epithelial_index);
			triangle_of_nodes.push_back(first_neighbour);
			triangle_of_nodes.push_back(second_neighbour);

			//Sort the node indices, as this makes it easier to check if
			//we have already found this triangle
			std::sort(triangle_of_nodes.begin(), triangle_of_nodes.end());

			//Add sorted indices to vector
			for (unsigned i = 0; i < triangle_of_nodes.size(); i++)
			{
				epithelial_triangle[i] = triangle_of_nodes[i];
			}

			//Add triple to vector of triangle, if it hasn't been added previously
			if(!HasThisTriangleAlreadyBeenFound(nodes_in_triangles, epithelial_triangle))
			{
				nodes_in_triangles.push_back(epithelial_triangle);
			}
		}
	}

	return nodes_in_triangles;
}

/*
 * Method to obtain list of adjacent epithelial triangle pairs.
 * List contains non-overlapping triangles, so if a triangle is already
 * part of a pair, it will not be considered again.
 * Structure of vector is < <nodeA, commonEdge>, <nodeB, commonEdge> >
 */
std::vector<std::pair<c_vector<unsigned, 3>, c_vector<unsigned, 3> > > NodeBasedBasementMembraneForceWithNearestNeighbours::GetAdjacentEpithelialTrianglePairs(AbstractCellPopulation<2>& rCellPopulation, std::map<unsigned, std::pair<unsigned, unsigned> > epithelialIndicesAndNeighbours)
{
	//Initialise vector
	std::vector<std::pair<c_vector<unsigned, 3>, c_vector<unsigned, 3> > > adjacent_triangles;

		for (std::map<unsigned, std::pair<unsigned, unsigned> >::iterator map_iter = epithelialIndicesAndNeighbours.begin();
				map_iter != epithelialIndicesAndNeighbours.end();
				map_iter++)
		{
			//Get the node index and its neighbours
			unsigned current_node_index = map_iter->first;
			std::pair<unsigned, unsigned> current_neighbours = map_iter->second;

			std::pair<unsigned, std::pair<unsigned, unsigned> > current_node_and_neighbours = std::make_pair(current_node_index, current_neighbours);
			//If another node has ths neighbours, sweep the map for the other node
			if(DoesAnotherNodeHaveTheseNeighbours(epithelialIndicesAndNeighbours, current_node_and_neighbours))
			{
				for (std::map<unsigned, std::pair<unsigned, unsigned> >::iterator second_map_iter = epithelialIndicesAndNeighbours.begin();
						second_map_iter != epithelialIndicesAndNeighbours.end();
						second_map_iter++)
				{
					unsigned considered_node_index = second_map_iter->first;
					std::pair<unsigned, unsigned> considered_neighbours = second_map_iter->second;

					//If this node is different and has the same neighbours, then add the pair of triples to the vector
					if( (current_node_index != considered_node_index)&&(((current_neighbours.first == considered_neighbours.first)
							&&(current_neighbours.second == considered_neighbours.second))||((current_neighbours.first == considered_neighbours.second)
							&&(current_neighbours.second == considered_neighbours.first))) )
					{
						//Initialise two triangle vectors
						c_vector<unsigned, 3> first_triangle, second_triangle;

						//We need to order the nodes in the right way so that the check torques the nodes in the
						//correct manner
						//Get the location of the current node and its neighbours
						c_vector<double, 2> current_node_location = rCellPopulation.GetNode(current_node_index)->rGetLocation();
						c_vector<double, 2> left_node_location = rCellPopulation.GetNode(current_neighbours.first)->rGetLocation();
						c_vector<double, 2> right_node_location = rCellPopulation.GetNode(current_neighbours.second)->rGetLocation();

						//Get the vectors from the current node to its neighbours
						c_vector<double, 2> current_node_to_left_neighbour = rCellPopulation.rGetMesh().GetVectorFromAtoB(current_node_location, left_node_location);
						c_vector<double, 2> current_node_to_right_neighbour = rCellPopulation.rGetMesh().GetVectorFromAtoB(current_node_location, right_node_location);

						//Get the "signed distance" between the two vectors
						double signed_distance = inner_prod(current_node_to_left_neighbour, current_node_to_right_neighbour);

						//If the signed distance is positive, then the angle is acute and we want to torque the common neighbours
						if (signed_distance > 0.0)
						{
							first_triangle[0] = current_node_index; second_triangle[0] = considered_node_index;
							if (current_neighbours.first < current_neighbours.second)
							{
								first_triangle[1] = current_neighbours.first; first_triangle[2] = current_neighbours.second;
								second_triangle[1] = current_neighbours.first; second_triangle[2] = current_neighbours.second;
							}
							else
							{
								first_triangle[1] = current_neighbours.second; first_triangle[2] = current_neighbours.first;
								second_triangle[1] = current_neighbours.second; second_triangle[2] = current_neighbours.first;
							}
						}
						else
						{
							first_triangle[0] = current_neighbours.first; second_triangle[0] = current_neighbours.second;
							if (current_node_index < considered_node_index)
							{
								first_triangle[1] = current_node_index; first_triangle[2] = considered_node_index;
								second_triangle[1] = current_node_index; second_triangle[2] = considered_node_index;
							}
							else
							{
								first_triangle[1] = considered_node_index; first_triangle[2] = current_node_index;
								second_triangle[1] = considered_node_index; second_triangle[2] = current_node_index;
							}
						}

						//If neither triangle has already been determined to be part of a pair
						if ( (!IsThisTriangleAlreadyInAPair(adjacent_triangles, first_triangle))
								&&(!IsThisTriangleAlreadyInAPair(adjacent_triangles, second_triangle)) )
						{
							//Initialise pair
							std::pair<c_vector<unsigned, 3>, c_vector<unsigned, 3> > triangle_pair;

							//Sort the pair by the first index
							if(current_neighbours.first < current_neighbours.second)
							{
								triangle_pair = std::make_pair(first_triangle, second_triangle);
							}
							else
							{
								triangle_pair = std::make_pair(second_triangle, first_triangle);
							}
							adjacent_triangles.push_back(triangle_pair);
						}
					}
				}
			}
		}

	return adjacent_triangles;
}

/*
 * Method to check if a node shares its neighbours with another node
 */
bool NodeBasedBasementMembraneForceWithNearestNeighbours::DoesAnotherNodeHaveTheseNeighbours(std::map<unsigned, std::pair<unsigned, unsigned> > epithelialIndicesAndNeighbours, std::pair<unsigned, std::pair<unsigned, unsigned> > nodeAndItsNeighbours)
{
	//Initialise boolean
	bool does_another_node_have_these_neighbours = false;

	//Get the node index and the pair of neighbours
	unsigned epithelial_node_index = nodeAndItsNeighbours.first;
	std::pair<unsigned, unsigned> epithelial_neighbour_pair = nodeAndItsNeighbours.second;

	//Iterate over each node and its neighbours
	for (std::map<unsigned, std::pair<unsigned, unsigned> >::iterator map_iter = epithelialIndicesAndNeighbours.begin();
			map_iter != epithelialIndicesAndNeighbours.end();
			map_iter++)
	{
		unsigned current_node_index = map_iter->first;
		std::pair<unsigned, unsigned> current_neighbour_pair = map_iter->second;

		//If there is another node with the same neighbour pair as our concerned node
		if( (current_node_index != epithelial_node_index)&&( ((current_neighbour_pair.first == epithelial_neighbour_pair.first)
				&&(current_neighbour_pair.second == epithelial_neighbour_pair.second))||((current_neighbour_pair.first == epithelial_neighbour_pair.second)
				&&(current_neighbour_pair.second == epithelial_neighbour_pair.first)) ) )
		{
			does_another_node_have_these_neighbours = true;
			break; //Break the for loop (don't need to check any others)
		}

	}

	return does_another_node_have_these_neighbours;
}

/*
 * Method to get each node index, the indices of their neighbours and the lengths between the neighbours
 * for a triangle, i.e. all the information we need to perform the necessary checks
 */
std::vector<std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> > NodeBasedBasementMembraneForceWithNearestNeighbours::GetRelevantInformationAboutTriangle(AbstractCellPopulation<2>& rCellPopulation, c_vector<unsigned, 3> epithelialTriangle)
{
	//Declare vector
	std::vector<std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> > information_about_triangle;

	//Declare vector for lengths
	std::vector<double> lengths;

	//Go through each node index
	for (unsigned i = 0; i < epithelialTriangle.size(); i++)
	{
		unsigned current_index = epithelialTriangle[i];

		//Get the other two neighbours
		std::vector<unsigned> neighbours_in_triangle;

		for (unsigned j = 0; j < epithelialTriangle.size(); j++)
		{
			if( epithelialTriangle[j] != current_index) //Only consider the other indices
			{
				neighbours_in_triangle.push_back(epithelialTriangle[j]);
			}
		}

		std::sort(neighbours_in_triangle.begin(), neighbours_in_triangle.end());

		//Get the distance between the two neighbours
		unsigned first_neighbour = neighbours_in_triangle[0];
		c_vector<double, 2> first_neighbour_location = rCellPopulation.GetNode(first_neighbour)->rGetLocation();

		unsigned second_neighbour = neighbours_in_triangle[1];
		c_vector<double, 2> second_neighbour_location = rCellPopulation.GetNode(second_neighbour)->rGetLocation();

		//Get the length of the edge using GetVectorFromAtoB to account for periodicities etc
		c_vector<double, 2> second_to_first = rCellPopulation.rGetMesh().GetVectorFromAtoB(second_neighbour_location, first_neighbour_location);
		double length_of_edge = norm_2(second_to_first);

		//Create the pair of neighbours
		std::pair<unsigned, unsigned> two_neighbours = std::make_pair(first_neighbour, second_neighbour);

		//Make the pair of the node and its neighbours
		std::pair<unsigned, std::pair<unsigned, unsigned> > node_and_its_neighbours = std::make_pair(current_index, two_neighbours);

		//Make the final pair of the 'triple' of the node and its neighbours and the distance between the neighbours
		std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> information_about_node = std::make_pair(node_and_its_neighbours, length_of_edge);

		//Add information to vector
		information_about_triangle.push_back(information_about_node);

		//Add lengths to vector
		lengths.push_back(length_of_edge);
	}

	//Sort the vector of lengths
	std::sort(lengths.begin(), lengths.end());

	//Declare vector with sorted information
	std::vector<std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> > sorted_information_about_triangle;

	for (unsigned i = 0; i < lengths.size(); i++)
	{
		for (unsigned j = 0; j < information_about_triangle.size(); j++)
		{
			std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> node_information = information_about_triangle[j];

			if( (node_information.second == lengths[i])&&(!HasThisInformationAlreadyBeenCollected(sorted_information_about_triangle, node_information)) )
			{
				sorted_information_about_triangle.push_back(node_information); //If considered node has corresponding distance
			}
		}
	}

	//Return vector
	return sorted_information_about_triangle;
}

/*
 * Method to get each node index, the indices of their neighbours and the lengths between the neighbours
 * for a quadrilateral, i.e. all the information we need to perform the necessary checks.
 * Recall the structure of the triangle pair:
 * < <nodeA, commonEdge>, <nodeB, commonEdge> >
 */
std::vector<std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> > NodeBasedBasementMembraneForceWithNearestNeighbours::GetRelevantInformationAboutQuadrilateral(AbstractCellPopulation<2>& rCellPopulation, std::pair<c_vector<unsigned, 3>, c_vector<unsigned, 3> > epithelialTrianglePair)
{
	//Declare vector
	std::vector<std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> > information_about_quadrilateral;

	//Declare vector for lengths
	std::vector<double> lengths;

	//Get the two triangles
	c_vector<unsigned, 3> first_triangle = epithelialTrianglePair.first;
	c_vector<unsigned, 3> second_triangle = epithelialTrianglePair.second;

	//Obtain the relevant information from the first triangle
	for (unsigned i = 1; i < first_triangle.size(); i++) //Only consider the endpoints of the common edge
	{
		unsigned current_index = first_triangle[i];

		//Get the other two neighbours
		std::vector<unsigned> neighbours_in_triangle;

		for (unsigned j = 0; j < first_triangle.size(); j++)
		{
			if( first_triangle[j] != current_index) //Only consider the other indices
			{
				neighbours_in_triangle.push_back(first_triangle[j]);
			}
		}

		std::sort(neighbours_in_triangle.begin(), neighbours_in_triangle.end()); //Sort the neighbours vector

		//Get the distance between the two neighbours
		unsigned first_neighbour = neighbours_in_triangle[0];
		c_vector<double, 2> first_neighbour_location = rCellPopulation.GetNode(first_neighbour)->rGetLocation();

		unsigned second_neighbour = neighbours_in_triangle[1];
		c_vector<double, 2> second_neighbour_location = rCellPopulation.GetNode(second_neighbour)->rGetLocation();

		//Get length of edge. We do it using GetVectorFromAtoB to account for periodic meshes
		c_vector<double, 2> second_to_first = rCellPopulation.rGetMesh().GetVectorFromAtoB(second_neighbour_location, first_neighbour_location);
		double length_of_edge = norm_2(second_to_first);

		//Create the pair of neighbours
		std::pair<unsigned, unsigned> two_neighbours = std::make_pair(first_neighbour, second_neighbour);

		//Make the pair of the node and its neighbours
		std::pair<unsigned, std::pair<unsigned, unsigned> > node_and_its_neighbours = std::make_pair(current_index, two_neighbours);

		//Make the final pair of the 'triple' of the node and its neighbours and the distance between the neighbours
		std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> information_about_node = std::make_pair(node_and_its_neighbours, length_of_edge);

		//Add information to vector
		information_about_quadrilateral.push_back(information_about_node);

		//Add lengths to vector
		lengths.push_back(length_of_edge);
	}

	//Go through each node index
	for (unsigned i = 1; i < second_triangle.size(); i++)
	{
		unsigned current_index = second_triangle[i];

		//Get the other two neighbours
		std::vector<unsigned> neighbours_in_triangle;

		for (unsigned j = 0; j < second_triangle.size(); j++)
		{
			if( second_triangle[j] != current_index) //Only consider the other indices
			{
				neighbours_in_triangle.push_back(second_triangle[j]);
			}
		}

		std::sort(neighbours_in_triangle.begin(), neighbours_in_triangle.end()); //Sort the neighbours

		//Get the distance between the two neighbours
		unsigned first_neighbour = neighbours_in_triangle[0];
		c_vector<double, 2> first_neighbour_location = rCellPopulation.GetNode(first_neighbour)->rGetLocation();

		unsigned second_neighbour = neighbours_in_triangle[1];
		c_vector<double, 2> second_neighbour_location = rCellPopulation.GetNode(second_neighbour)->rGetLocation();

		//Get the length of the edge
		c_vector<double, 2> second_to_first = rCellPopulation.rGetMesh().GetVectorFromAtoB(second_neighbour_location, first_neighbour_location);
		double length_of_edge = norm_2(second_to_first);

		//Create the pair of neighbours
		std::pair<unsigned, unsigned> two_neighbours = std::make_pair(first_neighbour, second_neighbour);

		//Make the pair of the node and its neighbours
		std::pair<unsigned, std::pair<unsigned, unsigned> > node_and_its_neighbours = std::make_pair(current_index, two_neighbours);

		//Make the final pair of the 'triple' of the node and its neighbours and the distance between the neighbours
		std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> information_about_node = std::make_pair(node_and_its_neighbours, length_of_edge);

		//Add information to vector
		information_about_quadrilateral.push_back(information_about_node);

		//Add lengths to vector
		lengths.push_back(length_of_edge);
	}

	//Sort the vector
	std::sort(lengths.begin(), lengths.end());

	//Declare vector with sorted information
	std::vector<std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> > sorted_information_about_quadrilateral;

	for (unsigned i = 0; i < lengths.size(); i++)
	{
		for (unsigned j = 0; j < information_about_quadrilateral.size(); j++)
		{
			std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> node_information = information_about_quadrilateral[j];

			if( (node_information.second == lengths[i])&&(!HasThisInformationAlreadyBeenCollected(sorted_information_about_quadrilateral, node_information)) )
			{
				sorted_information_about_quadrilateral.push_back(node_information); //If considered node has corresponding distance
			}
		}
	}

	//Return vector
	return sorted_information_about_quadrilateral;
}

/*
 * Method to obtain the opposite edge to a specified edge in a quadrilateral,
 * in order to appropriately re-assign neighbours
 */
std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> NodeBasedBasementMembraneForceWithNearestNeighbours::GetOppositeQuadrilateralEdge(std::vector<std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> > quadrilateralInformation, std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> consideredNodeAndNeighbours)
{
	std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> node_with_opposite_edge;

	for (unsigned i = 0; i < quadrilateralInformation.size(); i++)
	{
		std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> current_node_and_neighbours = quadrilateralInformation[i];

		//The opposite edge will be completely 'disjoint' to the node and neighbours
		if ( (consideredNodeAndNeighbours.first.first != current_node_and_neighbours.first.first)&&(consideredNodeAndNeighbours.first.second.first != current_node_and_neighbours.first.second.first)
				&&(consideredNodeAndNeighbours.first.second.second != current_node_and_neighbours.first.second.second)&&(consideredNodeAndNeighbours.first.second.second != current_node_and_neighbours.first.second.first)
				&&(consideredNodeAndNeighbours.first.second.first != current_node_and_neighbours.first.second.second) )
		{
			node_with_opposite_edge = current_node_and_neighbours;
			break; //Break out of for loop
		}
	}

	return node_with_opposite_edge;
}

/*
 * Method to check if node index and its neighbours form a triangle
 */
bool NodeBasedBasementMembraneForceWithNearestNeighbours::IsNodePartOfATriangle(std::map<unsigned, std::pair<unsigned, unsigned> > epithelialIndicesAndNeighbours, unsigned epithelialNodeIndex)
{
	//Declare boolean
	bool is_node_part_of_triangle = false;

	//Get the node's neighbours (first pair to check)
	std::pair<unsigned, unsigned> first_pair_of_neighbours = epithelialIndicesAndNeighbours[epithelialNodeIndex];

	//Get the indices of the neighbours
	unsigned first_neighbour = first_pair_of_neighbours.first;
	unsigned second_neighbour = first_pair_of_neighbours.second;

	//If either if the neighbours has two neighbours
	if( (epithelialIndicesAndNeighbours.find(first_neighbour) != epithelialIndicesAndNeighbours.end())
			||(epithelialIndicesAndNeighbours.find(second_neighbour) != epithelialIndicesAndNeighbours.end()) )
	{
		if( epithelialIndicesAndNeighbours.find(first_neighbour) != epithelialIndicesAndNeighbours.end() ) //If first neighbour also has two neighbours
		{
			std::pair<unsigned, unsigned> second_pair_of_neighbours = epithelialIndicesAndNeighbours[first_neighbour]; //Get the first neighbour's neighbours

			if( ((second_pair_of_neighbours.first == epithelialNodeIndex)&&(second_pair_of_neighbours.second == second_neighbour))||((second_pair_of_neighbours.first == second_neighbour)&&(second_pair_of_neighbours.second == epithelialNodeIndex)) ) //If the neighbours coincide with the original node and second neighbour
			{
//				//If the nodes coincide, check the second neighbour's neighbours
//				if( epithelialIndicesAndNeighbours.find(second_neighbour) != epithelialIndicesAndNeighbours.end() ) //If the second neighbour also has two neighbours
//				{
//					std::pair<unsigned, unsigned> third_pair_of_neighbours = epithelialIndicesAndNeighbours[second_neighbour]; //Get the second neighbour's neighbours
//
//					if( ((third_pair_of_neighbours.first == epithelialNodeIndex)&&(third_pair_of_neighbours.second == first_neighbour))||((third_pair_of_neighbours.first == first_neighbour)&&(third_pair_of_neighbours.second == epithelialNodeIndex)) ) //If the neighbours coincide with the original node and first neighbour
//					{
						is_node_part_of_triangle = true;
//					}
//				}
			}
		}
		else if (epithelialIndicesAndNeighbours.find(second_neighbour) != epithelialIndicesAndNeighbours.end()) //Else if the second neighbour has two neighbours
		{
			std::pair<unsigned, unsigned> second_pair_of_neighbours = epithelialIndicesAndNeighbours[second_neighbour]; //Get the first neighbour's neighbours

			if( ((second_pair_of_neighbours.first == epithelialNodeIndex)&&(second_pair_of_neighbours.second == first_neighbour))||((second_pair_of_neighbours.first == first_neighbour)&&(second_pair_of_neighbours.second == epithelialNodeIndex)) ) //If the neighbours coincide with the original node and second neighbour
			{
//				//If the nodes coincide, check the second neighbour's neighbours
//				if( epithelialIndicesAndNeighbours.find(second_neighbour) != epithelialIndicesAndNeighbours.end() ) //If the second neighbour also has two neighbours
//				{
//					std::pair<unsigned, unsigned> third_pair_of_neighbours = epithelialIndicesAndNeighbours[second_neighbour]; //Get the second neighbour's neighbours
//
//					if( ((third_pair_of_neighbours.first == epithelialNodeIndex)&&(third_pair_of_neighbours.second == first_neighbour))||((third_pair_of_neighbours.first == first_neighbour)&&(third_pair_of_neighbours.second == epithelialNodeIndex)) ) //If the neighbours coincide with the original node and first neighbour
//					{
						is_node_part_of_triangle = true;
//					}
//				}
			}
		}

	}

	//Return boolean variable
	return is_node_part_of_triangle;
}

/*
 * Method to check if an E-E-E triangle has already been
 * considered to be part of an E-E-E triangle pair. This method ensures
 * the pairs are not overlapping in terms of the triangles.
 */
bool NodeBasedBasementMembraneForceWithNearestNeighbours::IsThisTriangleAlreadyInAPair(std::vector<std::pair<c_vector<unsigned, 3>, c_vector<unsigned, 3> > > foundTrianglePairs, c_vector<unsigned, 3> currentTriangle)
{
	//Declare boolean value
	bool has_this_triangle_already_been_found = false;

	for (unsigned i = 0; i < foundTrianglePairs.size(); i++)
	{
		std::pair<c_vector<unsigned, 3>, c_vector<unsigned, 3> > epithelial_triangle_pair = foundTrianglePairs[i];

		//If either the first or second triangle in the pair is the considered triangle, declare the boolean to be true
		if ( ((epithelial_triangle_pair.first[0]==currentTriangle[0])&&(epithelial_triangle_pair.first[1]==currentTriangle[1])&&(epithelial_triangle_pair.first[2]==currentTriangle[2]))
				||((epithelial_triangle_pair.second[0]==currentTriangle[0])&&(epithelial_triangle_pair.second[1]==currentTriangle[1])&&(epithelial_triangle_pair.second[2]==currentTriangle[2])) )
		{
			has_this_triangle_already_been_found = true;
			break; //Break for loop, as we no longer need to consider other triangle pairs
		}

	}

	return has_this_triangle_already_been_found; //Return boolean
}
/*
 * Method to check if an E-E-E triangle has been found previously
 */
bool NodeBasedBasementMembraneForceWithNearestNeighbours::HasThisTriangleAlreadyBeenFound(std::vector<c_vector<unsigned, 3> > foundTriangles, c_vector<unsigned, 3> currentTriangle)
{
	bool has_triangle_already_been_found = false;

	//Iterate through each 'triple' of nodes
	for (unsigned i = 0; i < foundTriangles.size(); i++)
	{
		c_vector<unsigned, 3> considered_triangle = foundTriangles[i];

		if( (considered_triangle[0] == currentTriangle[0])&&(considered_triangle[1] == currentTriangle[1])
				&&(considered_triangle[2] == currentTriangle[2]) ) //If each index is the same
		{
			has_triangle_already_been_found = true;
			break; //As soon as we find one, we don't need to check the others
		}
	}

	return has_triangle_already_been_found;
}

/*
 * Criteria check to properly sort information vector for triangles
 * and quadrilaterals.
 */

bool NodeBasedBasementMembraneForceWithNearestNeighbours::HasThisInformationAlreadyBeenCollected(std::vector<std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> > currentInformationCollection, std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> consideredInformation)
{
	bool has_this_already_been_collected = false;

	//Get all the pieces of the considered information
	unsigned considered_node = consideredInformation.first.first;
	std::pair<unsigned, unsigned> considered_neighbours = consideredInformation.first.second;
	double considered_length = consideredInformation.second;

	for (unsigned i = 0; i < currentInformationCollection.size(); i++)
	{
		//Get current information
		std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> current_information_piece = currentInformationCollection[i];

		unsigned current_node = current_information_piece.first.first;
		std::pair<unsigned, unsigned> current_neighbours = current_information_piece.first.second;
		double current_length = current_information_piece.second;

		if( (current_node == considered_node)&&(current_neighbours.first == considered_neighbours.first)
				&&(current_neighbours.second == considered_neighbours.second)&&(considered_length == current_length) )
		{
			has_this_already_been_collected = true;
			break;
		}
	}

	return has_this_already_been_collected;
}

/*
 * Method to calculate the centre of mass of a considered index's  set of neighbouring indices.
 * We have to do it this weird way to account for periodic meshes, but one can prove (I did) that
 * this gives the same centre of mass. Essentially, we calculate the average 'relative vector' from the node
 * to its neighbours and add that to the reference point.
 */
c_vector<double, 2> NodeBasedBasementMembraneForceWithNearestNeighbours::GetCentreOfMassOfNeighbours(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex)
{
	//Create vector
	c_vector<double,2> centre_of_mass = zero_vector<double>(2);

	//Get the considered node's location
	c_vector<double,2> reference_location = rCellPopulation.GetNode(nodeIndex)->rGetLocation();

	//Get the set of neighbouring indices

	std::vector<unsigned> neighbouringIndices = GetNeighbouringEpithelialIndices(rCellPopulation, nodeIndex);

	//Get number of nodes
	unsigned num_nodes = neighbouringIndices.size();

	//Iterate over nodes
	for (unsigned i = 0; i < neighbouringIndices.size(); i++)
	{
		unsigned neighbour_index = neighbouringIndices[i];

		//Get location of node index
		c_vector<double,2> neighbour_location = rCellPopulation.GetNode(neighbour_index)->rGetLocation();

		//Get the relative position between the reference point and the neighbouring point.
		//This is the main reason we need to define the function this way: to use GetVectorFromAtoB
		c_vector<double,2> reference_to_neighbour = rCellPopulation.rGetMesh().GetVectorFromAtoB(reference_location, neighbour_location);

		centre_of_mass += reference_to_neighbour; //Add to centre of mass
	}

	centre_of_mass /= num_nodes; //Average vectors

	//Add the reference point
	centre_of_mass += reference_location;

	return centre_of_mass;
}

/*
 * Method to calculate the 'orientation' of the relative vector
 */
double NodeBasedBasementMembraneForceWithNearestNeighbours::GetRelativeAngleWithRespectToVerticalAxis(AbstractCellPopulation<2>& rCellPopulation, c_vector<double, 2> firstPoint, c_vector<double, 2> secondPoint)
{
	//Get relative vector
	c_vector<double, 2> relative_vector = rCellPopulation.rGetMesh().GetVectorFromAtoB(secondPoint, firstPoint);

	//Get x and y coordinates
	double x = relative_vector[0];
	double y = relative_vector[1];

	double angle = atan(x/y); //Calculate initial angle with respect to vertical axis

	//We now adjust for the quadrant that relative_vector sits in. For the 1st and 2nd quadrant,
	//the angle is correct so we only have to adjust for the 3rd and 4th quadrant.

	if ( (x<0.0)&&(y<0.0) ) //3rd quadrant
	{
		angle -= M_PI;
	}
	else if ( (x>0.0)&&(y<0.0) ) //4th quadrant
	{
		angle += M_PI;
	}
	else if ( (x<0.0)&&(y==0.0)) //If we're sitting on the x-axis, in the 2nd and 3rd quadrant
	{
		angle = -0.5*M_PI;
	}
	else if ( (x>0.0)&&(y==0.0)) //If we're sitting on the x-axis, in the 1st and 4th quadrant
	{
		angle = 0.5*M_PI;
	}
	else if ( (x==0.0)&&(y==0.0)) //If the two points coincide
	{
		angle = 0.0;
	}

	return angle;
}

/*
 * Method to calculate the average oientation of a given node index's epithelial neighbours.
 * Returns a vector (cos(theta), sin(theta)), where theta is the average angle of orientation.
 */
c_vector<double, 2> NodeBasedBasementMembraneForceWithNearestNeighbours::GetAveragePlaneOfOrientationOfNeighbours(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex)
{
	//Initialise the orientation vector
	c_vector<double, 2> average_orientation;

	//Initialise the average angle of orientation
	double angle_of_orientation = 0.0;

	//Get the neighbour set of the node
	std::vector<unsigned> neighbouring_epithelial_indices = GetNeighbouringEpithelialIndices(rCellPopulation, nodeIndex);

	//Get the number of neighbours
	unsigned num_epithelial_neighbours = neighbouring_epithelial_indices.size();

	//Gte the centre of mass of the neighbours
	c_vector<double, 2> centre_of_mass_of_neighbours = GetCentreOfMassOfNeighbours(rCellPopulation, nodeIndex);

	for (unsigned i = 0; i < num_epithelial_neighbours; i++)
	{
		//Get the index
		unsigned neighbour_index = neighbouring_epithelial_indices[i];

		//Get the location
		c_vector<double, 2> neighbour_location = rCellPopulation.GetNode(neighbour_index)->rGetLocation();

		//Get relative angle from centre of mass to neighbour location
		double theta = GetRelativeAngleWithRespectToVerticalAxis(rCellPopulation, neighbour_location, centre_of_mass_of_neighbours);

		angle_of_orientation += theta; //Update angle of orientation
	}

	//Average the angle
	angle_of_orientation /= num_epithelial_neighbours;

	//Define the vector
	average_orientation[0] = cos(angle_of_orientation);
	average_orientation[1] = sin(angle_of_orientation);

	return average_orientation;
}

/*
 * Method to calculate the orthogonal projection of firstPoint onto secondPoint
 */
c_vector<double, 2> NodeBasedBasementMembraneForceWithNearestNeighbours::GetOrthogonalProjection(AbstractCellPopulation<2>& rCellPopulation, c_vector<double, 2> firstPoint, c_vector<double, 2> secondPoint)
{
	double second_norm = norm_2(secondPoint);
	c_vector<double,2> unit_second = secondPoint/second_norm; //Get the unit vector of secondPoint

	double scalar_proj = inner_prod(firstPoint, unit_second); //Get scalar projection

	c_vector<double,2> vector_proj = scalar_proj*unit_second; //Get vector projection

	//Get orthogonal projection
	c_vector<double,2> orthog_proj = rCellPopulation.rGetMesh().GetVectorFromAtoB(vector_proj, firstPoint);

	return orthog_proj;
}

//Method to calculate the force due to the basement membrane on an epithelial cell
c_vector<double, 2> NodeBasedBasementMembraneForceWithNearestNeighbours::CalculateForceDueToBasementMembrane(AbstractCellPopulation<2>& rCellPopulation, std::map<unsigned, std::pair<unsigned, unsigned > > epithelialIndicesAndNeighbours, unsigned nodeIndex)
{
	//Get the left and right crypt boundaries and the target curvature
	double left_boundary = GetLeftCryptBoundary();
	double right_boundary = GetRightCryptBoundary();
	double target_curvature = GetTargetCurvature();

	//Get the neighbours of the node index
	std::pair<unsigned, unsigned> epithelial_neighbours = epithelialIndicesAndNeighbours[nodeIndex];
	unsigned left_node_index = epithelial_neighbours.first;
	unsigned right_node_index = epithelial_neighbours.second;

	//Get the location of all three points for the curvature calculation
	c_vector<double, 2> centre_point = rCellPopulation.GetNode(nodeIndex)->rGetLocation();
	c_vector<double, 2> left_point = rCellPopulation.GetNode(left_node_index)->rGetLocation();
	c_vector<double, 2> right_point = rCellPopulation.GetNode(right_node_index)->rGetLocation();

	double basement_membrane_parameter = GetBasementMembraneParameter(); //Get the basement membrane stiffness

	double curvature = FindParametricCurvature(rCellPopulation, left_point, centre_point, right_point);

	/*
	 * We need to make sure the curvature isn't too large. Degenerate cases occur during proliferation events, when newly-born
	 * cells may get compressed, resulting in a massive curvature
	 */
	if (fabs(curvature) > 5.8647) //This is the curvature value for the degenerate case \/, when arc lengths are 1 and the left and right cells are 0.1 apart, mimicking a proliferation event
	{
		if (curvature < 0.0)
		{
			curvature = -5.8647;
		}
		else
		{
			curvature = 5.8647;
		}
	}

	//Get the unit vectors from the centre points to its left and right neighbours
	c_vector<double, 2> centre_to_left = rCellPopulation.rGetMesh().GetVectorFromAtoB(centre_point, left_point);
	centre_to_left /= norm_2(centre_to_left); //Normalise vector

	c_vector<double, 2> centre_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(centre_point, right_point);
	centre_to_right /= norm_2(centre_to_right); //Normalise vector

	/* Define direction using a formula derived (essentially solve for the vector w such that u.w = v.w, where w is the force direction
	* and u = unit vector to the left and v = unit vector pointing to the right)
	*/
	c_vector<double, 2> left_to_right = rCellPopulation.rGetMesh().GetVectorFromAtoB(centre_to_left, centre_to_right);
	assert(norm_2(left_to_right) != 0.0);

	//Define force direction
	c_vector<double, 2> force_direction;
	force_direction(0) = left_to_right[1];
	force_direction(1) = -left_to_right[0];

	force_direction /= norm_2(left_to_right);

	/* We now ensure the vector is pointing in the appropriate direction
	 * (it will always point "down" and "outwards" initially).
	 * */
	if (target_curvature > 0.0) //If the force has overshot the target curvature, we need to reverse the force direction
	{
		if ( (curvature > 0.0)&&(curvature - target_curvature > 0.0) ) //If points look like V and the 'v' is too pointy, we send it away from the CoM
		{
			force_direction *= -1.0;
		}

	}
	else if (target_curvature < 0.0) //Similar situation but with "/\"
	{
		if ( (curvature < 0.0)&&(curvature - target_curvature < 0.0) )
		{
			force_direction *= -1.0;
		}
	}
	else //Reverse the force direction if we get a "V"
	{
		if (curvature > 0.0)
		{
			force_direction *= -1.0;
		}
	}

	//Initialise the force vector
	c_vector<double, 2> force_due_to_basement_membrane;

	//Again, teh geometry of the model alters how we apply target curvature
	bool is_force_applied_to_crypt = GetCryptGeometryCheck();

	//If we are considering a crypt geometry
	if(is_force_applied_to_crypt)
	{
		//If the cell falls in the region of non-zero target curvature
		if( (centre_point[0] > left_boundary)&&(centre_point[0] < right_boundary) )
		{
			force_due_to_basement_membrane = basement_membrane_parameter*(fabs(curvature - target_curvature) )*force_direction;
		}
		else
		{
			//We take the absolute value of the local curvature as the force direction vector already accounts for which way
			//the epithelial node 'should' go
			force_due_to_basement_membrane = basement_membrane_parameter*fabs(curvature)*force_direction;
		}
	}
	else //Else we are modelling organoid
	{
		force_due_to_basement_membrane = basement_membrane_parameter*(fabs(curvature - target_curvature) )*force_direction;
	}

	if (norm_2(force_due_to_basement_membrane) > 2.0/0.005) //Hacked check to make sure the force isn't excessively large
	{
		PRINT_2_VARIABLES(left_point[0], left_point[1]);
		PRINT_2_VARIABLES(centre_point[0], centre_point[1]);
		PRINT_2_VARIABLES(right_point[0], right_point[1]);
		PRINT_VARIABLE(curvature);
	}
	return force_due_to_basement_membrane;
}

//Method overriding the virtual method for AbstractForce. The crux of what really needs to be done.
void NodeBasedBasementMembraneForceWithNearestNeighbours::AddForceContribution(AbstractCellPopulation<2>& rCellPopulation)
{
	//Get the left and right crypt boundaries and the target curvature
	double left_boundary = GetLeftCryptBoundary();
	double right_boundary = GetRightCryptBoundary();
	double target_curvature = GetTargetCurvature();
	double basement_membrane_stiffness = GetBasementMembraneParameter();

	//Get the epithelial nodes that have a left and right nearest neighbour
	std::map<unsigned, std::pair<unsigned, unsigned> > epithelial_nodes_and_their_neighbours = GetEpithelialIndicesAndTheirLeftAndRightEpithelialNeighbours(rCellPopulation);

	/*
	 * Need to check two things:
	 * 1) the map of nodes and their neighbours defines a confluent monolayer
	 * 2) every epithelial index in the tissue population is mapped to a pair of neighbours in the map
	 */

	//Get the epithelial indices
	std::vector<unsigned> epithelial_indices = GetEpithelialIndices(rCellPopulation);

	//If the sizes are different, we have not accounted for every epithelial layer in the population
	if (epithelial_indices.size() != epithelial_nodes_and_their_neighbours.size())
	{
		EXCEPTION("Not every epithelial cell has been accounted for in defined monolayer");
	}

	//Check confluence of monolayer
	for (std::map<unsigned, std::pair<unsigned, unsigned> >::iterator map_iter = epithelial_nodes_and_their_neighbours.begin();
			map_iter != epithelial_nodes_and_their_neighbours.end();
			map_iter++)
	{

		unsigned centre_node = map_iter->first;

		std::pair<unsigned, unsigned> centre_nodes_neighbours = map_iter->second;
		unsigned left_node = centre_nodes_neighbours.first;
		unsigned right_node = centre_nodes_neighbours.second;

		//Get the neighbours of the left and right neighbours
		std::pair<unsigned, unsigned> left_nodes_neighbours = epithelial_nodes_and_their_neighbours[left_node];
		std::pair<unsigned, unsigned> right_nodes_neighbours = epithelial_nodes_and_their_neighbours[right_node];

		/*
		 * If the right neighbour of the left node is NOT the centre node, OR
		 * the left neighbour of the right node is NOT the centre node, the layer is not confluent.
		 * If this condition doesn't hold for any of the epithelial nodes, then the layer is confluent.
		 */

		if (left_nodes_neighbours.second != centre_node)
		{
			//We would like to know where the discontinuity is: we print the indices of the discontinuous triple
			//and the locations of the relevant nodes
			PRINT_5_VARIABLES(left_node, left_nodes_neighbours.second, centre_node, right_nodes_neighbours.first, right_node);
			PRINT_5_VARIABLES(rCellPopulation.GetNode(left_node)->rGetLocation()[0], rCellPopulation.GetNode(left_nodes_neighbours.second)->rGetLocation()[0],
					rCellPopulation.GetNode(centre_node)->rGetLocation()[0], rCellPopulation.GetNode(right_nodes_neighbours.first)->rGetLocation()[0],
					rCellPopulation.GetNode(right_node)->rGetLocation()[0]);
			PRINT_5_VARIABLES(rCellPopulation.GetNode(left_node)->rGetLocation()[1], rCellPopulation.GetNode(left_nodes_neighbours.second)->rGetLocation()[1],
					rCellPopulation.GetNode(centre_node)->rGetLocation()[1], rCellPopulation.GetNode(right_nodes_neighbours.first)->rGetLocation()[1],
					rCellPopulation.GetNode(right_node)->rGetLocation()[1]);

			//Flag the error
			EXCEPTION("Defined monolayer is not confluent.");
		}
		else if (right_nodes_neighbours.first != centre_node)
		{
			//We would like to know where the discontinuity is: we print the indices of the discontinuous triple
			//and the locations of the relevant nodes
			PRINT_5_VARIABLES(left_node, left_nodes_neighbours.second, centre_node, right_nodes_neighbours.first, right_node);
			PRINT_5_VARIABLES(rCellPopulation.GetNode(left_node)->rGetLocation()[0], rCellPopulation.GetNode(left_nodes_neighbours.second)->rGetLocation()[0],
					rCellPopulation.GetNode(centre_node)->rGetLocation()[0], rCellPopulation.GetNode(right_nodes_neighbours.first)->rGetLocation()[0],
					rCellPopulation.GetNode(right_node)->rGetLocation()[0]);
			PRINT_5_VARIABLES(rCellPopulation.GetNode(left_node)->rGetLocation()[1], rCellPopulation.GetNode(left_nodes_neighbours.second)->rGetLocation()[1],
					rCellPopulation.GetNode(centre_node)->rGetLocation()[1], rCellPopulation.GetNode(right_nodes_neighbours.first)->rGetLocation()[1],
					rCellPopulation.GetNode(right_node)->rGetLocation()[1]);

			//Flag the error
				EXCEPTION("Defined monolayer is not confluent.");
		}


	}

//	//Check for E-E-E-E 'quadrilaterals' and remove them
//	std::map<unsigned, std::pair<unsigned, unsigned> > epithelial_nodes_and_their_neighbours_without_quadrilaterals = CheckForAdjacentEpithelialTriangles(rCellPopulation, epithelial_nodes_and_their_neighbours);
//
//	//Check for E-E-E triangles and remove them
//	std::map<unsigned, std::pair<unsigned, unsigned> > epithelial_nodes_and_their_neighbours_without_triangles = CheckForEpithelialTriangles(rCellPopulation, epithelial_nodes_and_their_neighbours_without_quadrilaterals);

	//Iterate over each node and its neighbours
	for (std::map<unsigned, std::pair<unsigned, unsigned> >::iterator map_iter = epithelial_nodes_and_their_neighbours.begin();
			map_iter != epithelial_nodes_and_their_neighbours.end();
			map_iter++)
	{

		//Get the node and its neighbours
		unsigned centre_node_index = map_iter->first;

		std::pair<unsigned, unsigned> epithelial_neighbours = map_iter->second;
		unsigned left_node_index = epithelial_neighbours.first;
		unsigned right_node_index = epithelial_neighbours.second;

		//Calculate the forces exerted on the left, centre and right nodes by the basement membrane
		c_vector<double, 2> force_on_centre_node = CalculateForceDueToBasementMembrane(rCellPopulation, epithelial_nodes_and_their_neighbours, centre_node_index);

		assert(norm_2(force_on_centre_node) < 2.0/0.005); //0.005 is the imposed dt, while 2.0 is the absolute movement threshold

		rCellPopulation.GetNode(centre_node_index)->AddAppliedForceContribution(force_on_centre_node);
	}
}

void NodeBasedBasementMembraneForceWithNearestNeighbours::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile <<  "\t\t\t<BasementMembraneParameter>"<<  mBasementMembraneParameter << "</BasementMembraneParameter> \n" ;
	*rParamsFile <<  "\t\t\t<TargetCurvature>" << mTargetCurvature << "</TargetCurvature> \n";
	*rParamsFile <<  "\t\t\t<LeftBoundary>"<<  mLeftBoundary << "</LeftBoundary> \n" ;
	*rParamsFile <<  "\t\t\t<RightBoundary>"<<  mRightBoundary << "</RightBoundary> \n" ;
	*rParamsFile <<  "\t\t\t<UsePositionDependentMembraneForce>"<<  mUsePositionDependentMembraneForce << "</UsePositionDependentMembraneForce> \n" ;
	*rParamsFile <<  "\t\t\t<MembraneForceMultiplier>"<<  mMembraneForceMultiplier << "</MembraneForceMultiplier> \n" ;

	// Call direct parent class
	AbstractForce<2>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(NodeBasedBasementMembraneForceWithNearestNeighbours)

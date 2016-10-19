#ifndef NODEBASEDBASEMENTMEMBRANEFORCEWITHNEARESTNEIGHBOURS_HPP_
#define NODEBASEDBASEMENTMEMBRANEFORCEWITHNEARESTNEIGHBOURS_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "MeshBasedCellPopulation.hpp"

#include "RandomNumberGenerator.hpp"


#include <cmath>
#include <list>
#include <fstream>

/**
 * MODIFIED FOR PERSONAL USE BY AXEL ALMET
 * A force class that defines the force due to the basement membrane.
 */

class NodeBasedBasementMembraneForceWithNearestNeighbours : public AbstractForce<2>
{
    friend class TestCrossSectionModelInteractionForce;

private :

    /** Parameter that multiplies the curvature to give the basement membrane force */
    double mBasementMembraneParameter;

    /** Target curvature for the ring of cells (NodeBased) */
    double mTargetCurvature;

    /** x-coordinate that encloses the region in which to apply a non-zero target curvature */
    double mLeftBoundary;

    /** x-coordinate that encloses the region in which to apply a non-zero target curvature */
    double mRightBoundary;

    /** Make the basement membrane force dependent on the position of a cell up the crypt */
    bool mUsePositionDependentMembraneForce;

    /** Boolean to check whether force is applied to crypt-like epithelium or organoid */
    bool mApplyForceToCrypt;

    /** The multiplication factor for the basement membrane parameter */
    double mMembraneForceMultiplier;

    /** The cut off radius for defining neighbouring nodes */
    double mCutOffRadius;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractForce<2> >(*this);
        archive & mBasementMembraneParameter;
        archive & mTargetCurvature;
        archive & mLeftBoundary;
        archive & mRightBoundary;
        archive & mApplyForceToCrypt;
        archive & mUsePositionDependentMembraneForce;
        archive & mMembraneForceMultiplier;
        archive & mCutOffRadius;
    }

public :

    /**
     * Constructor.
     */
	NodeBasedBasementMembraneForceWithNearestNeighbours();

    /**
     * Destructor.
     */
    ~NodeBasedBasementMembraneForceWithNearestNeighbours();

    /* Set method for Basement Membrane Parameter
     */
    void SetBasementMembraneParameter(double basementMembraneParameter);

    /* Get method for Basement Membrane Parameter
     */
    double GetBasementMembraneParameter();

    /* Value of curvature at crypt base */
    void SetTargetCurvature(double targetCurvature = 0.0);

    /*
     * Get method for Target Curvature Parameters
     */
    double GetTargetCurvature();

    /*
     * Set method for left crypt boundary
     */
    void SetLeftCryptBoundary(double leftBoundary);

    /*
     * Get method for Left crypt boundary parameter
     */
    double GetLeftCryptBoundary();

    /*
     * Set method for right crypt boundary parameter
     */
    void SetRightCryptBoundary(double rightBoundary);

    /*
     * Get method for right crypt boundary
     */
    double GetRightCryptBoundary();

    /*
     * Set method for geometry-dependent basement membrane for application, i.e.crypt or organoid
     */
    void SetCryptGeometry(bool applyForceToCrypt = true);

    /*
     * Check to determine whether or not basement membrane force is to be applied to a crypt or organoid
     */
    bool GetCryptGeometryCheck();

    /* Set method for position-dependent basement membrane force multiplier (i.e. if you want to apply a different
     * basement membrane parameter in the crypt base, or at the orifice)
     * @param usePositionDependentMembraneForce whether to multiply the basement membrane force by a factor
     * @param membraneForceMultiplier the multiplication factor for the basement membrane force
     */
    void SetPositionDependentMultiplier(bool usePositionDependentMembraneForce = false, double membraneForceMultiplier = 1.0);

    /* Get method for basement membrane force strength multiplier
     */
    double GetPositionDependentMultiplier();

    /*Get method for cut off radius */
    double GetCutOffRadius();

    /*Set method for cut off radius */
    void SetCutOffRadius(double cutOffRadius);

    /* Removing duplicated entries of a vector
     */
    void RemoveDuplicates1D(std::vector<unsigned>& rVectorWithDuplicates);

    double FindParametricCurvature(AbstractCellPopulation<2>& rCellPopulation,
    								c_vector<double, 2> leftPoint,
									c_vector<double, 2> centrePoint,
									c_vector<double, 2> rightPoint);

    /*
     * Method to get the indices of the monolayer
     */
    std::vector<unsigned> GetEpithelialIndices(AbstractCellPopulation<2>& rCellPopulation);

    /*
     * Get method for the crypt height extremes, adapted from Dunn et al. (2012).
     * Returns a vector with two components, first component is the y-coordinate of
     * the crypt base, second component is the y-coordinate of the crypt collar.
     */
    c_vector<double, 2> GetCryptHeightExtremes(AbstractCellPopulation<2>& rCellPopulation);

    /*
     * Method to return the nodes connected to a particular node within a defined cut-off
     * radius
     */
    std::vector<unsigned> GetNeighbouringEpithelialIndices(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex);

    /* Finding the y-coordinates of the crypt orifice and crypt base
     * The first entry of the resulting vector is the orifice, the second is the base
     */
//    c_vector<double,2> GetCryptHeightExtremes(AbstractCellPopulation<2>& rCellPopulation);

    /*
     * Method to return the set of left and right neighbours.
     * Require this to avoid any weird issues where a neighbour can technically be
     * both 'left and right'.
     */
    std::pair<std::vector<unsigned>, std::vector<unsigned> > GetLeftAndRightNeighbours(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex);

    /*
     * Method to get the epithelial nodes and their left and right neighbours, if they have any
     */
    std::map<unsigned, std::pair<unsigned, unsigned> > GetEpithelialIndicesAndTheirLeftAndRightEpithelialNeighbours(AbstractCellPopulation<2>& rCellPopulation);

    /*
     * Method to remove pairs of adjacent eptihelial-epithelial triangles of nodes
     * and their nearest neighbours, if they exist
     */
    std::map<unsigned, std::pair<unsigned, unsigned> > CheckForAdjacentEpithelialTriangles(AbstractCellPopulation<2>& rCellPopulation, std::map<unsigned, std::pair<unsigned, unsigned> > epithelialIndicesAndNeighbours);

    /*
     * Method to remove epithelial-epithelial triangles of nodes and their
     * nearest neighbours and to re-assign nearest neighbours if such triangles exist
     */
    std::map<unsigned, std::pair<unsigned, unsigned> > CheckForEpithelialTriangles(AbstractCellPopulation<2>& rCellPopulation, std::map<unsigned, std::pair<unsigned, unsigned> > epithelialIndicesAndNeighbours);

    /*
     * Method to obtain list of epithelial triangles, if they exist
     */
    std::vector<c_vector<unsigned, 3> > GetEpithelialTriangles(std::map<unsigned, std::pair<unsigned, unsigned> > epithelialIndicesAndNeighbours);

    /*
     * Method to obtain list of adjacent epithelial triangle pairs.
     * List contains non-overlapping triangles, so if a triangle is already
     * part of a pair, it will not be considered again.
     * Structure of vector is < <nodeA, commonEdge>, <nodeB, commonEdge> >.
     */
    std::vector<std::pair<c_vector<unsigned, 3>, c_vector<unsigned, 3> > > GetAdjacentEpithelialTrianglePairs(AbstractCellPopulation<2>& rCellPopulation, std::map<unsigned, std::pair<unsigned, unsigned> > epithelialIndicesAndNeighbours);

    /*
     * Method to check if more than one node has the same pair of nearest neighbours
     */
    bool DoesAnotherNodeHaveTheseNeighbours(std::map<unsigned, std::pair<unsigned, unsigned> > epithelialIndicesAndNeighbours, std::pair<unsigned, std::pair<unsigned, unsigned> > nodeAndItsNeighbours);

    /*
     * Method to obtain all relevant information about an E-E-E triangle:
     * each vertex, its neighbours and the lengths between the neighbours
     */
    std::vector<std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> > GetRelevantInformationAboutTriangle(AbstractCellPopulation<2>& rCellPopulation, c_vector<unsigned, 3> epithelialTriangle);

    /*
     * Method to obtain all relevant information about a quadrilateral formed
     * by two adjacent E-E-E triangles: each vertex, its neighbours and the lengths
     * between the neighbours
     */
    std::vector<std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> > GetRelevantInformationAboutQuadrilateral(AbstractCellPopulation<2>& rCellPopulation, std::pair<c_vector<unsigned, 3>, c_vector<unsigned, 3> > epithelialTrianglePair);

    /*
     * Method to obtain the opposite edge to a specified edge in a quadrilateral,
     * in order to appropriately re-assign neighbours
     */
    std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> GetOppositeQuadrilateralEdge(std::vector<std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> > quadrilateralInformation, std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> consideredNodeAndNeighbours);

    /*
     * Method to see if the considered node and its nearest neighbours form
     * a triangle of connections
     */
    bool IsNodePartOfATriangle(std::map<unsigned, std::pair<unsigned, unsigned > > epithelialIndicesAndNeighbours, unsigned epithelialNodeIndex);

    /*
     * Method to check if an E-E-E triangle has been found
     * previously in a pair
     */
    bool IsThisTriangleAlreadyInAPair(std::vector<std::pair<c_vector<unsigned, 3>, c_vector<unsigned, 3> > > foundTrianglePairs, c_vector<unsigned, 3> currentTriangle);

    /*
     * Method to check if an E-E-E triangle has been found previously
     */
    bool HasThisTriangleAlreadyBeenFound(std::vector<c_vector<unsigned, 3> > foundTriangles, c_vector<unsigned, 3> currentTriangle);

    /*
     * Criteria check to properly sort information vector for triangles
     * and quadrilaterals.
     */
    bool HasThisInformationAlreadyBeenCollected(std::vector<std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> > currentInformationCollection, std::pair<std::pair<unsigned, std::pair<unsigned, unsigned> >, double> consideredInformation);

    /*
     * Method to calculate centre of mass of a given index's neighbours with respect to the considered node
     */
    c_vector<double, 2> GetCentreOfMassOfNeighbours(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex);

    /*
     * Method to calculate the angle of the relative vector firstPoint - secondPoint with respect to the vertical y-axis
     */
    double GetRelativeAngleWithRespectToVerticalAxis(AbstractCellPopulation<2>& rCellPopulation, c_vector<double, 2> firstPoint, c_vector<double, 2> secondPoint);

    /*
     * Method to calculate the average orientation of a given epithelial node's neighbour set
     */
    c_vector<double, 2> GetAveragePlaneOfOrientationOfNeighbours(AbstractCellPopulation<2>& rCellPopulation, unsigned nodeIndex);

    /*
     * Method to calculate the unit orthogonal projection of firstPoint onto secondPoint
     */
    c_vector<double, 2> GetOrthogonalProjection(AbstractCellPopulation<2>& rCellPopulation, c_vector<double, 2> firstPoint, c_vector<double, 2> secondPoint);

    /*
     * Method to calculate the force due to basement membrane on an epithelial cell
     */
    c_vector<double, 2> CalculateForceDueToBasementMembrane(AbstractCellPopulation<2>& rCellPopulation, std::map<unsigned, std::pair<unsigned, unsigned > > epithelialIndicesAndNeighbours, unsigned nodeIndex);

    /**
     * Overridden AddForceContribution method.
     *
     * @param rCellPopulation reference to the tissue
     */
    void AddForceContribution(AbstractCellPopulation<2>& rCellPopulation);

    /**
     * Outputs force Parameters to file
	 *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(NodeBasedBasementMembraneForceWithNearestNeighbours)

#endif /*NODEBASEDBASEMENTMEMBRANEFORCE_HPP_*/

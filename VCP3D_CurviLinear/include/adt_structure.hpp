/*!
 * \file adt_structure.hpp
 * \brief Headers of the subroutines for carrying out geometrical searches using an
 *        alternating digital tree (ADT).
 *        The subroutines and functions are in the <i>adt_structure.cpp</i> file.
 * \author E. van der Weide
 * \version 7.0.2 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

/*! 
 * \class CADTComparePointClass
 * \brief  Functor, used for the sorting of the points when building an ADT.
 * \author E. van der Weide
 */
class CADTComparePointClass {
private:
  const su2double      *pointCoor;      /*!< \brief Pointer to the coordinates of the points. */
  const unsigned short splitDirection;  /*!< \brief Split direction used in the sorting. */
  const unsigned short nDim;            /*!< \brief Number of spatial dimensions stored in the coordinates. */

public:
  /*!
   * \brief Constructor of the class. The member variables are initialized.
   * \param[in] coor      Pointer to the coordinates of the points.
   * \param[in] splitDir  Direction that must be used to sort the coordinates.
   * \param[in] nDimADT   Number of spatial dimensions of the ADT and coordinates.
   */
  CADTComparePointClass(const su2double      *coor,
                        const unsigned short splitDir,
                        const unsigned short nDimADT);
  /*!
   * \brief Destructor, nothing to be done.
   */
  ~CADTComparePointClass();

  /*!
   * \brief Operator used for the sorting of the points.
   * \param[in] p0  Index of the first point to be compared.
   * \param[in] p1  Index of the second point to be compared.
   */
  bool operator()(const unsigned long p0,
                  const unsigned long p1) const;
private:
  /*!
   * \brief Default constructor of the class, disabled.
   */
  CADTComparePointClass();
};

/*! 
 * \class CADTNodeClass
 * \brief  Class for storing the information needed in a node of an ADT.
 * \author E. van der Weide
 */
class CADTNodeClass {
public:
  bool          childrenAreTerminal[2];  /*!< \brief Whether or not the child leaves are terminal. */
  unsigned long children[2];             /*!< \brief Child leaves. If childrenAreTerminal is true the children
                                                     contain the point ID's or bounding box ID's. Note that it
                                                     is allowed that one child is termimal and the other is not. */
  unsigned long centralNodeID;           /*!< \brief ID of a node, which is near the center of the leaf. */

  su2double *xMin;  /*!< \brief The minimum coordinates of this leaf. It points to a position in the large
                                vector, which contains the coordinates of all leaves. */
  su2double *xMax;  /*!< \brief The maximum coordinates of this leaf. It points to a position in the large
                                vector, which contains the coordinates of all leaves. */

  /*!
   * \brief Constructor of the class. Nothing to be done.
   */
  CADTNodeClass();

  /*!
   * \brief Destructor of the class. Nothing to be done.
   */
  ~CADTNodeClass();

  /*!
   * \brief Copy constructor of the class.
   * \param[in] other  Object from which the data must be copied.
   */
  CADTNodeClass(const CADTNodeClass &other);

  /*!
   * \brief Assignment operator.
   * \param[in] other  Object from which the data must be copied.
   */
  CADTNodeClass& operator=(const CADTNodeClass &other);

private:

  /*!
   * \brief Copy function, which copies the data from the given object.
   * \param[in] other  Object from which the data must be copied.
   */
  void Copy(const CADTNodeClass &other);
};

/*! 
 * \class CADTBaseClass
 * \brief  Base class for storing an ADT in an arbitrary number of dimensions.
 * \author E. van der Weide
 */
class CADTBaseClass {
protected:
  unsigned long nLeaves;    /*!< \brief Number of leaves in the ADT. */
  unsigned short nDimADT;   /*!< \brief Number of dimensions of the ADT. */
  bool           isEmpty;   /*!< \brief Whether or not the ADT is empty. */

  std::vector<CADTNodeClass> leaves; /*!< \brief Vector, which contains all the leaves of the ADT. */

private:
  std::vector<su2double> coorMinLeaves; /*!< \brief Vector, which contains all the minimum coordinates
                                               of the leaves. */
  std::vector<su2double> coorMaxLeaves; /*!< \brief Vector, which contains all the maximum coordinates
                                               of the leaves. */
protected:
  /*!
   * \brief Constructor of the class. Nothing to be done.
   */
  CADTBaseClass();

  /*!
   * \brief Destructor of the class. Nothing to be done.
   */
  virtual ~CADTBaseClass();  

  /*!
   * \brief Function, which builds the ADT of the given coordinates.
   * \param[in] nDim    Number of dimensions of the ADT.
   * \param[in] nPoints Number of points present in the ADT.
   * \param[in] coor    Coordinates of the points.
   */
  void BuildADT(unsigned short  nDim,
                unsigned long   nPoints,
                const su2double *coor);
public:
  /*!
   * \brief Function, which returns whether or not the ADT is empty.
   * \return  Whether or not the ADT is empty.
   *
   */
  bool IsEmpty(void) const;

private:
  /*!
   * \brief Copy constructor of the class, disabled.
   */
  CADTBaseClass(const CADTBaseClass &);

  /*!
   * \brief Assignment operator, disabled.
   */
  CADTBaseClass& operator=(const CADTBaseClass &);
};

/*! 
 * \class CADTElemClass
 * \brief  Class for storing an ADT of (linear) elements in an arbitrary number of dimensions.
 * \author E. van der Weide
 * \version 7.0.2 "Blackbird"
 */
class CADTElemClass : public CADTBaseClass {
private:
  unsigned short nDim; /*!< \brief Number of spatial dimensions. */

  std::vector<su2double>     coorPoints;    /*!< \brief Vector, which contains the coordinates
                                                   of the points in the ADT. */
  std::vector<su2double>     BBoxCoor;      /*!< \brief Vector, which contains the coordinates
                                                   of the bounding boxes of the elements. */

  std::vector<unsigned short> elemVTK_Type; /*!< \brief Vector, which the type of the elements
                                                   using the VTK convention. */
  std::vector<unsigned long> nDOFsPerElem;  /*!< \brief Vector, which contains the number of DOFs
                                                   of the elements in the ADT in cumulative
                                                   storage format. */
  std::vector<unsigned long> elemConns;     /*!< \brief Vector, which contains the connectivities
                                                   of the elements in the ADT. */
  std::vector<unsigned short> localMarkers; /*!< \brief Vector, which contains the marker ID's
                                                   of the elements in the ADT. */
  std::vector<unsigned long> localElemIDs;  /*!< \brief Vector, which contains the local element ID's
                                                    of the elements in the ADT. */
  std::vector<int>           ranksOfElems;  /*!< \brief Vector, which contains the ranks
                                                    of the elements in the ADT. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in]     val_nDim     Number of spatial dimensions of the problem.
   * \param[in]     val_coor     Coordinates of the local points to be stored in the ADT.
   * \param[in,out] val_connElem Local connectivity of the elements to be stored in the ADT.
                                 In parallel mode the connectivities are corrected for the
                                 offset in node numbers of the rank.
   * \param[in]     val_VTKElem  Type of the elements using the VTK convention.
   * \param[in]     val_markerID Markers of the local elements.
   * \param[in]     val_elemID   Local element IDs of the elements. 
   * \param[in]     globalTree   Whether or not a global tree must be built. If false
                                 a local ADT is built.
   */
  CADTElemClass(unsigned short              val_nDim,
                std::vector<su2double>      &val_coor,
                std::vector<unsigned long>  &val_connElem,
                std::vector<unsigned short> &val_VTKElem,
                std::vector<unsigned short> &val_markerID,
                std::vector<unsigned long>  &val_elemID);

  /*!
   * \brief Destructor of the class. Nothing to be done.
   */
  ~CADTElemClass();

  /*!
   * \brief Function, which determines the element that contains the given
            coordinate.
   * \param[in]  coor             Coordinate which the element must contain.
   * \param[out] markerID         Local marker ID of the element containing the coordinate.
   * \param[out] elemID           Local element ID of the element containing the coordinate.
   * \param[out] rankID           Rank on which element containing the coordinate is stored.
   * \param[out] parCoor          Parametric coordinates of coor inside the element,
                                  which contains the coordinate.
   * \param[out] weightsInterpol  Interpolation weigts of of coor inside the element,
                                  which contains the coordinate.
   * \return                      True if an element is found, false if not.
   */
  bool DetermineContainingElement(const su2double *coor,
                                  unsigned short  &markerID,
                                  unsigned long   &elemID,
                                  int             &rankID,
                                  su2double       *parCoor,
                                  su2double       *weightsInterpol);

private:

  /*!
   * \brief Function, which checks whether or not the given coordinate is
            inside the given element.
   * \param[in]  elemID          ID of the element for which the containment
                                 must be checked.
   * \param[in]  coor            Coordinate for which the containment must be checked.
   * \param[out] parCoor         Parametric coordinates of coor if it is inside the
                                 given element.
   * \param[out] weightsInterpol Interpolation weights of coor if it is inside the
                                 given element.
   * \return                     True if coor is inside the element and false otherwise.
   */
  bool CoorInElement(const unsigned long elemID,
                     const su2double     *coor,
                     su2double           *parCoor,
                     su2double           *weightsInterpol);

  /*!
   * \brief Function, which checks whether or not the given coordinate is
            inside the given quadrilateral.
   * \param[in]  elemID          ID of the quadrilateral for which the containment
                                 must be checked.
   * \param[in]  coor            Coordinate for which the containment must be checked.
   * \param[out] parCoor         Parametric coordinates of coor if it is inside the
                                 given quadrilateral.
   * \param[out] weightsInterpol Interpolation weights of coor if it is inside the
                                 given quadrilateral.
   * \return                     True if coor is inside the quadrilateral and false otherwise.
   */
  bool CoorInQuadrilateral(const unsigned long elemID,
                           const su2double     *coor,
                           su2double           *parCoor,
                           su2double           *weightsInterpol);

  /*!
   * \brief Function, which checks whether or not the given coordinate is
            inside the given triangle.
   * \param[in]  elemID          ID of the triangle for which the containment
                                 must be checked.
   * \param[in]  coor            Coordinate for which the containment must be checked.
   * \param[out] parCoor         Parametric coordinates of coor if it is inside the
                                 given triangle.
   * \param[out] weightsInterpol Interpolation weights of coor if it is inside the
                                 given triangle.
   * \return                     True if coor is inside the triangle and false otherwise.
   */
  bool CoorInTriangle(const unsigned long elemID,
                      const su2double     *coor,
                      su2double           *parCoor,
                      su2double           *weightsInterpol);

  /*!
   * \brief Function, which checks whether or not the given coordinate is
            inside the given hexahedron.
   * \param[in]  elemID          ID of the hexahedron for which the containment
                                 must be checked.
   * \param[in]  coor            Coordinate for which the containment must be checked.
   * \param[out] parCoor         Parametric coordinates of coor if it is inside the
                                 given hexahedron.
   * \param[out] weightsInterpol Interpolation weights of coor if it is inside the
                                 given hexahedron.
   * \return                     True if coor is inside the hexahedron and false otherwise.
   */
  bool CoorInHexahedron(const unsigned long elemID,
                        const su2double     *coor,
                        su2double           *parCoor,
                        su2double           *weightsInterpol);

  /*!
   * \brief Function, which checks whether or not the given coordinate is
            inside the given prism.
   * \param[in]  elemID          ID of the prism for which the containment
                                 must be checked.
   * \param[in]  coor            Coordinate for which the containment must be checked.
   * \param[out] parCoor         Parametric coordinates of coor if it is inside the
                                 given prism.
   * \param[out] weightsInterpol Interpolation weights of coor if it is inside the
                                 given prism.
   * \return                     True if coor is inside the prism and false otherwise.
   */
  bool CoorInPrism(const unsigned long elemID,
                   const su2double     *coor,
                   su2double           *parCoor,
                   su2double           *weightsInterpol);

  /*!
   * \brief Function, which checks whether or not the given coordinate is
            inside the given pyramid.
   * \param[in]  elemID          ID of the pyramid for which the containment
                                 must be checked.
   * \param[in]  coor            Coordinate for which the containment must be checked.
   * \param[out] parCoor         Parametric coordinates of coor if it is inside the
                                 given pyramid.
   * \param[out] weightsInterpol Interpolation weights of coor if it is inside the
                                 given pyramid.
   * \return                     True if coor is inside the pyramid and false otherwise.
   */
  bool CoorInPyramid(const unsigned long elemID,
                     const su2double     *coor,
                     su2double           *parCoor,
                     su2double           *weightsInterpol);

  /*!
   * \brief Function, which checks whether or not the given coordinate is
            inside the given tetrahedron.
   * \param[in]  elemID          ID of the tetrahedron for which the containment
                                 must be checked.
   * \param[in]  coor            Coordinate for which the containment must be checked.
   * \param[out] parCoor         Parametric coordinates of coor if it is inside the
                                 given tetrahedron.
   * \param[out] weightsInterpol Interpolation weights of coor if it is inside the
                                 given tetrahedron.
   * \return                     True if coor is inside the tetrahedron and false otherwise.
   */
  bool CoorInTetrahedron(const unsigned long elemID,
                         const su2double     *coor,
                         su2double           *parCoor,
                         su2double           *weightsInterpol);

  /*!
   * \brief Function, which provides an initial guess for the parametric coordinates
            of the given point inside a hexahedron by splitting it into tetrahedra.
   * \param[in]  xRelC   Coordinates of the point to be investigated relative to
                         vertex 0 of the hexahedron.
   * \param[in]  xRel    Coordinates of the vertices of the hexahedron relative to
                         vertex 0.
   * \param[out] parCoor Initial guess of the parametric coordinates.
   * \return             True if the initial guess is within the hexahedron and
                         false otherwise.
   */
  bool InitialGuessContainmentHexahedron(const su2double xRelC[3],
                                         const su2double xRel[8][3],
                                         su2double       *parCoor);

  /*!
   * \brief Function, which provides an initial guess for the parametric coordinates
            of the given point inside a prism by splitting it into tetrahedra.
   * \param[in]  xRelC   Coordinates of the point to be investigated relative to
                         vertex 0 of the prism.
   * \param[in]  xRel    Coordinates of the vertices of the prism relative to
                         vertex 0.
   * \param[out] parCoor Initial guess of the parametric coordinates.
   * \return             True if the initial guess is within the prism and
                         false otherwise.
   */
  bool InitialGuessContainmentPrism(const su2double xRelC[3],
                                    const su2double xRel[6][3],
                                    su2double       *parCoor);

  /*!
   * \brief Function, which provides an initial guess for the parametric coordinates
            of the given point inside a pyramid by splitting it into tetrahedra.
   * \param[in]  xRelC   Coordinates of the point to be investigated relative to
                         vertex 0 of the pyramid.
   * \param[in]  xRel    Coordinates of the vertices of the pyramid relative to
                         vertex 0.
   * \param[out] parCoor Initial guess of the parametric coordinates.
   * \return             True if the initial guess is within the pyramid and
                         false otherwise.
   */
  bool InitialGuessContainmentPyramid(const su2double xRelC[3],
                                      const su2double xRel[5][3],
                                      su2double       *parCoor);

  /*!
   * \brief Default constructor of the class, disabled.
   */
  CADTElemClass();

  /*!
   * \brief Copy constructor of the class, disabled.
   */
  CADTElemClass(const CADTElemClass &);

  /*!
   * \brief Assignment operator, disabled.
   */
  CADTElemClass& operator=(const CADTElemClass &);
};

#include "adt_structure.inl"

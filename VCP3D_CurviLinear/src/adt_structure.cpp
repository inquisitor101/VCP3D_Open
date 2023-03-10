/*!
 * \file adt_structure.cpp
 * \brief Main subroutines for for carrying out geometrical searches using an
 *        alternating digital tree (ADT).
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

#include "VCP3D_CurviLinear.hpp"

/* Define the tolerance to decide whether or not a point is inside an element. */
const su2double tolInsideElem   =  (su2double) 1.e-10;
const su2double paramLowerBound = -one - tolInsideElem;
const su2double paramUpperBound =  one + tolInsideElem;

CADTComparePointClass::CADTComparePointClass(const su2double      *coor,
                                             const unsigned short splitDir,
                                             const unsigned short nDimADT)
  : pointCoor(coor),
    splitDirection(splitDir),
    nDim(nDimADT) {}

void CADTNodeClass::Copy(const CADTNodeClass &other) {

  childrenAreTerminal[0] = other.childrenAreTerminal[0];
  childrenAreTerminal[1] = other.childrenAreTerminal[1];

  children[0] = other.children[0];
  children[1] = other.children[1];

  centralNodeID = other.centralNodeID;

  xMin = other.xMin;
  xMax = other.xMax;
}

void CADTBaseClass::BuildADT(unsigned short  nDim,
                             unsigned long   nPoints,
                             const su2double *coor) {

  /*---  Determine the number of leaves. It can be proved that nLeaves equals
         nPoints-1 for an optimally balanced tree. Take the exceptional case of
         nPoints == 1 into account and return if the tree is empty. ---*/
  nDimADT = nDim;
  isEmpty = false;
  nLeaves = nPoints -1;
  if(nPoints <= 1) ++nLeaves;
  if(nLeaves == 0) {isEmpty = true; return;}

  /*--- Allocate the memory for the leaves of the ADT and the minimum and
        maximum coordinates of the leaves. Note that these coordinates are
        stored in one vector, rather than that memory is allocated for the
        individual leaves. ---*/
  leaves.resize(nLeaves);
  coorMinLeaves.resize(nDim*nLeaves);
  coorMaxLeaves.resize(nDim*nLeaves);

  /*--- Define the vectors, which control the subdivision of the leaves. ---*/
  unsigned long nn = (nPoints+1)/2;
  std::vector<unsigned long> pointIDs(nPoints), pointIDsNew(nPoints);
  std::vector<unsigned long> nPointIDs(nn+1),   nPointIDsNew(nn+1);
  std::vector<unsigned long> curLeaf(nn),       curLeafNew(nn);

  /*--------------------------------------------------------------------------*/
  /*---                 Building of the actual ADT                         ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Initialize the arrays pointIDs, nPointIDs and curLeaf such that all
        points belong to the root leaf. Also set the counters nLeavesToDivide
        and nLeavesTot.  ---*/
  nPointIDs[0] = 0; nPointIDs[1] = nPoints;
  curLeaf[0] = 0;

  for(unsigned long i=0; i<nPoints; ++i) pointIDs[i] = i;

  unsigned long nLeavesToDivide = 1, nLeavesTot = 1;

  /*--- Loop to subdivide the leaves. The division is such that the ADT is
        optimally balanced.  ---*/
  for(;;) {

    /* Criterion to exit the loop. */
    if(nLeavesToDivide == 0) break;

    /* Initializations for the next round of subdivisions. */
    unsigned long nLeavesToDivideNew = 0;
    nPointIDsNew[0] = 0;

    /*--- Loop over the current number of leaves to be divided. ---*/
    for(unsigned long i=0; i<nLeavesToDivide; ++i) {

      /* Store the number of points present in the leaf in nn and the
         current leaf number in mm. */
      nn = nPointIDs[i+1] - nPointIDs[i];
      unsigned long mm = curLeaf[i];

      /*--- Set the pointers for the coordinates of the leaf to the correct
            locations in the vectors coorMinLeaves and coorMaxLeaves and
            determine the bounding box coordinates of this leaf. ---*/
      leaves[mm].xMin = coorMinLeaves.data() + nDim*mm;
      leaves[mm].xMax = coorMaxLeaves.data() + nDim*mm;

      unsigned long ll = nDim*pointIDs[nPointIDs[i]];
      for(unsigned short l=0; l<nDim; ++l)
        leaves[mm].xMin[l] = leaves[mm].xMax[l] = coor[ll+l];

      for(unsigned long j=(nPointIDs[i]+1); j<nPointIDs[i+1]; ++j) {
        ll = nDim*pointIDs[j];
        for(unsigned short l=0; l<nDim; ++l) {
          leaves[mm].xMin[l] = std::min(leaves[mm].xMin[l], coor[ll+l]);
          leaves[mm].xMax[l] = std::max(leaves[mm].xMax[l], coor[ll+l]);
        }
      }

      /*--- Determine the split direction for this leaf. The splitting is done
            in such a way that isotropy is reached as quickly as possible.
            Hence the split direction is the direction of the largest dimension
            of the leaf. ---*/
      unsigned short splitDir= 0;
      su2double distMax = -1.0;
      for(unsigned short l=0; l<nDim; ++l) {
        const su2double dist = leaves[mm].xMax[l] - leaves[mm].xMin[l];
        if(dist > distMax) {distMax = dist; splitDir = l;}
      }

      /* Sort the points of the current leaf in increasing order. The sorting
         is based on the coordinate in the split direction, for which the
         functor CADTComparePointClass is used. */
      std::sort(pointIDs.data() + nPointIDs[i], pointIDs.data() + nPointIDs[i+1],
                CADTComparePointClass(coor, splitDir, nDim));

      /* Determine the index of the node, which is approximately central
         in this leave. */
      leaves[mm].centralNodeID = pointIDs[nPointIDs[i] + nn/2];

      /*--- Determine the situation of the leaf. It is either a terminal leaf
            or a leaf that must be subdivided. ---*/
      if(nn <= 2) {

        /* Terminal leaf. Store the ID's of the points as children and
           indicate that the children are terminal. */
        leaves[mm].children[0] = pointIDs[nPointIDs[i]];
        leaves[mm].children[1] = pointIDs[nPointIDs[i+1]-1];

        leaves[mm].childrenAreTerminal[0] = true;
        leaves[mm].childrenAreTerminal[1] = true;
      }
      else {

        /* The leaf must be divided. Determine the number of points in the
           left leaf. This number is at least 2. The actual number stored in kk
           is this number plus an offset. Also initialize the counter nfl, which
           is used to store the bounding boxes in the arrays for the new round. */
        unsigned long kk  = (nn+1)/2 + nPointIDs[i];
        unsigned long nfl = nPointIDsNew[nLeavesToDivideNew];

        /* Copy the ID's of the left points into pointIDsNew. Also update the
           corresponding entry in nPointIDsNew. */
        for(unsigned long k=nPointIDs[i]; k<kk; ++k)
          pointIDsNew[nfl++] = pointIDs[k];

        nPointIDsNew[nLeavesToDivideNew+1] = nfl;

        /* Store the leaf info in the tree and in the leafs for next round and
           update nLeavesToDivideNew and nLeavesTot. */
        leaves[mm].children[0]            = nLeavesTot;
        leaves[mm].childrenAreTerminal[0] = false;

        curLeafNew[nLeavesToDivideNew] = nLeavesTot;
        ++nLeavesToDivideNew;
        ++nLeavesTot;

        /*--- The right leaf will only be created if it has more than one point
              in it, i.e. if the original leaf has more than three points.
              If the new leaf only has one point in it, it is not created.
              Instead, the point is stored in the current leaf. ---*/
        if(nn == 3)
        {
          /* Only three points present in the current leaf. The right leaf is
             not created and the last point is stored as the second child of
             the current leaf. */
          leaves[mm].children[1]            = pointIDs[nPointIDs[i+1]-1];
          leaves[mm].childrenAreTerminal[1] = true;
        }
        else {

          /* More than 3 points are present and thus the right leaf is created.
             Copy the ID's from pointIDs into pointIDsNew and update the
             counters for the new round. */
          unsigned long nfr = nPointIDsNew[nLeavesToDivideNew];

          for(unsigned long k=kk; k<nPointIDs[i+1]; ++k)
            pointIDsNew[nfr++] = pointIDs[k];

          nPointIDsNew[nLeavesToDivideNew+1] = nfr;

          /* Store the leaf info in the tree and in the leaves for next round and
             update nLeavesToDivideNew and nLeavesTot. */
          leaves[mm].children[1]            = nLeavesTot;
          leaves[mm].childrenAreTerminal[1] = false;

          curLeafNew[nLeavesToDivideNew] = nLeavesTot;
          ++nLeavesToDivideNew;
          ++nLeavesTot;
        }
      }
    }

    /* Set the data for the next round. */
    nLeavesToDivide = nLeavesToDivideNew;

    for(unsigned long i=0; i<=nLeavesToDivide; ++i)           nPointIDs[i] = nPointIDsNew[i];
    for(unsigned long i=0; i< nLeavesToDivide; ++i)           curLeaf[i]   = curLeafNew[i];
    for(unsigned long i=0; i<nPointIDs[nLeavesToDivide]; ++i) pointIDs[i]  = pointIDsNew[i];
  }
}

CADTElemClass::CADTElemClass(unsigned short              val_nDim,
                             std::vector<su2double>      &val_coor,
                             std::vector<unsigned long>  &val_connElem,
                             std::vector<unsigned short> &val_VTKElem,
                             std::vector<unsigned short> &val_markerID,
                             std::vector<unsigned long>  &val_elemID) {

  /* Copy the dimension of the problem into nDim. */
  nDim = val_nDim;

  /*--- A local tree must be built. Copy the data from the arguments into the
        member variables and set the ranks to the rank of this processor. ---*/
  coorPoints   = val_coor;
  elemConns    = val_connElem;
  elemVTK_Type = val_VTKElem;
  localMarkers = val_markerID;
  localElemIDs = val_elemID;

  ranksOfElems.assign(elemVTK_Type.size(), rank);

 /*--- Determine the values of the vector nDOFsPerElem, which contains the
        number of DOFs per element in cumulative storage format. ---*/
  const unsigned long nElem = elemVTK_Type.size();
  nDOFsPerElem.resize(nElem+1);

  nDOFsPerElem[0] = 0;
  for(unsigned long i=0; i<nElem; ++i) {

    unsigned short nDOFs = 0;
    switch( elemVTK_Type[i] ) {
      case LINE:          nDOFs = 2; break;
      case TRIANGLE:      nDOFs = 3; break;
      case QUADRILATERAL: nDOFs = 4; break;
      case TETRAHEDRON:   nDOFs = 4; break;
      case PYRAMID:       nDOFs = 5; break;
      case PRISM:         nDOFs = 6; break;
      case HEXAHEDRON:    nDOFs = 8; break;
    }

    nDOFsPerElem[i+1] = nDOFsPerElem[i] + nDOFs;
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Create the bounding boxes of the elements. The coordinates ---*/
  /*---         of these bounding boxes can be interpreted as a point in a ---*/
  /*---         higher dimensional space. The ADT is built as a tree of    ---*/
  /*---         these points in this higher dimensional space.             ---*/
  /*--------------------------------------------------------------------------*/

  /* Allocate the memory for the bounding boxes of the elements. */
  BBoxCoor.resize(2*nDim*nElem);

  /*--- Loop over the elements to determine the minimum and maximum coordinates
        of the elements. These coordinates completely define the bounding box,
        which can be represented as a point in 2*nDim dimensional space. ---*/
  for(unsigned long i=0; i<nElem; ++i) {

    /* Set the pointers for the minimum and maximum coordinates for this
       bounding box. Initialize these coordinates to the first coordinate
       of the corresponding element. */
    su2double *BBMin = BBoxCoor.data() + 2*nDim*i;
    su2double *BBMax = BBMin + nDim;

    unsigned long    j  = nDOFsPerElem[i];
    const su2double *xP = coorPoints.data() + nDim*elemConns[j];

    for(unsigned short k=0; k<nDim; ++k)
      BBMin[k] = BBMax[k] = xP[k];

    /* Loop over the other vertices of the element to determine the
       minimum and maximum coordinates of the bounding box. */
    for(j=(nDOFsPerElem[i]+1); j<nDOFsPerElem[i+1]; ++j) {
      xP = coorPoints.data() + nDim*elemConns[j];

      for(unsigned short k=0; k<nDim; ++k) {
        BBMin[k] = std::min(BBMin[k], xP[k]);
        BBMax[k] = std::max(BBMax[k], xP[k]);
      }
    }

    /* Add a tolerance to the bounding box size, such that the overlap check in
       the tree traversal does not go wrong due to round off error. */
    for(unsigned short k=0; k<nDim; ++k) {
      const su2double lenScale = BBMax[k] - BBMin[k];
      const su2double tol      = std::max(1.e-25, 1.e-6*lenScale);

      BBMin[k] -= tol;
      BBMax[k] += tol;
    }
  }

  /* Build the ADT of the bounding boxes. */
  BuildADT(2*nDim, nElem, BBoxCoor.data());
}

bool CADTElemClass::DetermineContainingElement(const su2double *coor,
                                               unsigned short  &markerID,
                                               unsigned long   &elemID,
                                               int             &rankID,
                                               su2double       *parCoor,
                                               su2double       *weightsInterpol) {

  // Define the vectors to store the frontLeaves and reserve some memory.
  std::vector<unsigned long> frontLeaves, frontLeavesNew;
  frontLeaves.reserve(100);
  frontLeavesNew.reserve(100);

  /* Start at the root leaf of the ADT, i.e. initialize frontLeaves such that
     it only contains the root leaf. Make sure to wipe out any data from a
     previous search. */
  frontLeaves.clear();
  frontLeaves.push_back(0);

  /* Infinite loop of the tree traversal. */
  for(;;) {

    /* Initialize the new front, i.e. the front for the next round, to empty. */
    frontLeavesNew.clear();

    /* Loop over the leaves of the current front. */
    for(unsigned long i=0; i<frontLeaves.size(); ++i) {

      /* Store the current leaf a bit easier in ll and loop over its children. */
      const unsigned long ll = frontLeaves[i];
      for(unsigned short mm=0; mm<2; ++mm) {

        /* Determine whether this child contains a node or a leaf
           of the next level of the ADT. */
        unsigned long kk = leaves[ll].children[mm];
        if( leaves[ll].childrenAreTerminal[mm] ) {

          /* Child contains a bounding box. Check if the coordinate is
             inside the bounding box. */
          const su2double *coorBBMin = BBoxCoor.data() + nDimADT*kk;
          const su2double *coorBBMax = coorBBMin + nDim;

          bool coorIsInside = true;
          for(unsigned short k=0; k<nDim; ++k) {
            if(coor[k] < coorBBMin[k]) coorIsInside = false;
            if(coor[k] > coorBBMax[k]) coorIsInside = false;
          }

          if( coorIsInside ) {

            /* Coordinate is inside the bounding box. Check if it
               is also inside the corresponding element. If so,
               set the required information and return true. */
            if( CoorInElement(kk, coor, parCoor, weightsInterpol) ) {
              markerID = localMarkers[kk];
              elemID   = localElemIDs[kk];
              rankID   = ranksOfElems[kk];
              return true;
            }
          }
        }
        else {

          /* Child contains a leaf. If the coordinate is inside the leaf
             store the leaf for the next round. */
          const su2double *coorBBMin = leaves[kk].xMin;
          const su2double *coorBBMax = leaves[kk].xMax + nDim;

          bool coorIsInside = true;
          for(unsigned short k=0; k<nDim; ++k) {
            if(coor[k] < coorBBMin[k]) coorIsInside = false;
            if(coor[k] > coorBBMax[k]) coorIsInside = false;
          }

          if( coorIsInside ) frontLeavesNew.push_back(kk);
        }
      }
    }

    /*--- End of the loop over the current front. Copy the data from
          frontLeavesNew to frontLeaves for the next round. If the new front
          is empty the entire tree has been traversed and a break can be made
          from the infinite loop. ---*/
    frontLeaves = frontLeavesNew;
    if(frontLeaves.size() == 0) break;
  }

  /* If this point is reached, no element is found that contains the coordinate
     and false must be returned. */
  return false;
}

bool CADTElemClass::CoorInElement(const unsigned long elemID,
                                  const su2double     *coor,
                                  su2double           *parCoor,
                                  su2double           *weightsInterpol) {

  /*--- Make a distinction between the element types. ---*/
  switch( elemVTK_Type[elemID] ) {

    case TRIANGLE:      return CoorInTriangle(elemID, coor, parCoor, weightsInterpol);
    case QUADRILATERAL: return CoorInQuadrilateral(elemID, coor, parCoor, weightsInterpol);
    case TETRAHEDRON:   return CoorInTetrahedron(elemID, coor, parCoor, weightsInterpol);
    case PYRAMID:       return CoorInPyramid(elemID, coor, parCoor, weightsInterpol);
    case PRISM:         return CoorInPrism(elemID, coor, parCoor, weightsInterpol);
    case HEXAHEDRON:    return CoorInHexahedron(elemID, coor, parCoor, weightsInterpol);

    default:
      /* This should not happen. */
      Terminate("CADTElemClass::CoorInElement", __FILE__, __LINE__,
                "This should not happen");
      return false;
  }
}

bool CADTElemClass::CoorInTriangle(const unsigned long elemID,
                                   const su2double     *coor,
                                   su2double           *parCoor,
                                   su2double           *weightsInterpol) {

  /* Determine the indices of the three vertices of the triangle,
     multiplied by nDim (which is 2). This gives the position in the
     coordinate array where the coordinates of these points are stored. */
  unsigned long i0 = nDOFsPerElem[elemID];
  unsigned long i1 = i0 + 1, i2 = i0 + 2;
  i0 = nDim*elemConns[i0]; i1 = nDim*elemConns[i1]; i2 = nDim*elemConns[i2];

  /* Determine the coordinates relative to the vertex 0. */
  const su2double xc = coor[0] - coorPoints[i0];
  const su2double yc = coor[1] - coorPoints[i0+1];

  const su2double x1 = coorPoints[i1]   - coorPoints[i0];
  const su2double y1 = coorPoints[i1+1] - coorPoints[i0+1];

  const su2double x2 = coorPoints[i2]   - coorPoints[i0];
  const su2double y2 = coorPoints[i2+1] - coorPoints[i0+1];

  /* The triangle is parametrized by X-X0 = (r+1)*(X1-X0)/2 + (s+1)*(X2-X0)/2,
     r, s >= -1, r+s <= 0. As this is a containment search, the number of
     dimesions is 2. As a consequence, the parametric coordinates r and s
     can be solved easily. Note that X0 is 0 in the above expression, because
     the coordinates are relative to node 0. */
  const su2double detInv = 2.0/(x1*y2 - x2*y1);
  parCoor[0] = detInv*(xc*y2 - yc*x2) - 1.0;
  parCoor[1] = detInv*(yc*x1 - xc*y1) - 1.0;

  /* Check if the point resides within the triangle and compute the
     interpolation weights if it is. */
  bool coorIsInside = false;
  if((parCoor[0] >= paramLowerBound) && (parCoor[1] >= paramLowerBound) &&
     ((parCoor[0]+parCoor[1]) <= tolInsideElem)) {
    coorIsInside = true;

    weightsInterpol[0] = -0.5*(parCoor[0] + parCoor[1]);
    weightsInterpol[1] =  0.5*(parCoor[0] + 1.0);
    weightsInterpol[2] =  0.5*(parCoor[1] + 1.0);
  }

  /* Return the value of coorIsInside. */
  return coorIsInside;
}

bool CADTElemClass::CoorInQuadrilateral(const unsigned long elemID,
                                        const su2double     *coor,
                                        su2double           *parCoor,
                                        su2double           *weightsInterpol) {

  /* Definition of the maximum number of iterations in the Newton solver
     and the tolerance level. */
  const unsigned short maxIt = 50;
  const su2double tolNewton  = 1.e-10;

  /* Determine the indices of the four vertices of the quadrilatral,
     multiplied by nDim (which is 2). This gives the position in the
     coordinate array where the coordinates of these points are stored. */
  unsigned long i0 = nDOFsPerElem[elemID];
  unsigned long i1 = i0 + 1, i2 = i0 + 2, i3 = i0+3;
  i0 = nDim*elemConns[i0]; i1 = nDim*elemConns[i1];
  i2 = nDim*elemConns[i2]; i3 = nDim*elemConns[i3];

  /* Determine the coordinates relative to the vertex 0. */
  const su2double xc = coor[0] - coorPoints[i0];
  const su2double yc = coor[1] - coorPoints[i0+1];

  const su2double x1 = coorPoints[i1]   - coorPoints[i0];
  const su2double y1 = coorPoints[i1+1] - coorPoints[i0+1];

  const su2double x2 = coorPoints[i2]   - coorPoints[i0];
  const su2double y2 = coorPoints[i2+1] - coorPoints[i0+1];

  const su2double x3 = coorPoints[i3]   - coorPoints[i0];
  const su2double y3 = coorPoints[i3+1] - coorPoints[i0+1];

  /* The parametrization of the quadrilatral is nonlinear, which requires an
     iterative algorithm. Especially for highly skewed quadrilaterals, this
     could lead to convergence problems. Hence the quadrilateral is split into
     linear triangles to check if the point is actually within the quad.
     First check the triangle i0-i1-i3. See CoorInTriangle for more details
     on this test. */
  su2double detInv = 2.0/(x1*y3 - x3*y1);
  parCoor[0] = detInv*(xc*y3 - yc*x3) - 1.0;
  parCoor[1] = detInv*(yc*x1 - xc*y1) - 1.0;

  bool coorIsInside = false;
  if((parCoor[0] >= paramLowerBound) && (parCoor[1] >= paramLowerBound) &&
     ((parCoor[0]+parCoor[1]) <= tolInsideElem)) coorIsInside = true;

  /* Check triangle i2-i3-i1 if the coordinate is not inside i0-i1-i3. */
  if( !coorIsInside ) {

    /* Define the coordinates w.r.t. vertex 2 using the numbering i2-i3-i1. */
    const su2double xxc = xc - x2, yyc = yc - y2;
    const su2double xx1 = x3 - x2, yy1 = y3 - y2;
    const su2double xx3 = x1 - x2, yy3 = y1 - y2;

    /* Check if the coordinate is inside this triangle. */
    detInv = 2.0/(xx1*yy3 - xx3*yy1);
    parCoor[0] = detInv*(xxc*yy3 - yyc*xx3) - 1.0;
    parCoor[1] = detInv*(yyc*xx1 - xxc*yy1) - 1.0;

    if((parCoor[0] >= paramLowerBound) && (parCoor[1] >= paramLowerBound) &&
       ((parCoor[0]+parCoor[1]) <= tolInsideElem)) coorIsInside = true;

    /* Convert the parametric coordinates to the ones used by the
       quadrilatral i0-i1-i2-i3. They serve as initial guess below. */
    parCoor[0] = -parCoor[0];
    parCoor[1] = -parCoor[1];
  }

  /* If the coordinate is in neither triangle, return false. */
  if( !coorIsInside ) return false;

  /* The coordinate is inside the quadrilatral and an initial guess has been
     obtained by splitting the quad into two triangles. Carry out a Newton
     algorithm to obtain the true parametric coordinates.
     The quadrilateral is parametrized by
     X = {(1-r)*(1-s)*X0 + (1+r)*(1-s)*X1 + (1+r)*(1+s)*X2 + (1-r)*(1+s)*X3}/4,
     -1 <= r,s <= 1. As the coordinates are relative to X0, the first term
     drops from this equation. The nonlinear set of equations can be written as
     V0 - r*V1 - s*V2 - r*s*V3 = 0, where V0 = xc - (x1+x2+x3)/4,
     V1 = (x1+x2-x3)/4, V2 = (x2+x3-X1)/4, V3 = (x2-x1-x3)/4. First
     construct the vectors V0, V1, V2 and V3. */
  const su2double V0x = xc - 0.25*(x1+x2+x3), V0y = yc - 0.25*(y1+y2+y3);
  const su2double V1x =      0.25*(x1+x2-x3), V1y =      0.25*(y1+y2-y3);
  const su2double V2x =      0.25*(x2+x3-x1), V2y =      0.25*(y2+y3-y1);
  const su2double V3x =      0.25*(x2-x1-x3), V3y =      0.25*(y2-y1-y3);

  /* Loop over the maximum number of iterations. */
  unsigned short itCount;
  for(itCount=0; itCount<maxIt; ++itCount) {

    /* Compute the values of the nonlinear functions. */
    const su2double f0 = V0x - parCoor[0]*V1x - parCoor[1]*V2x - parCoor[0]*parCoor[1]*V3x;
    const su2double f1 = V0y - parCoor[0]*V1y - parCoor[1]*V2y - parCoor[0]*parCoor[1]*V3y;

    /* Compute the negative of the Jacobian matrix. */
    const su2double a00 = V1x + parCoor[1]*V3x, a01 = V2x + parCoor[0]*V3x;
    const su2double a10 = V1y + parCoor[1]*V3y, a11 = V2y + parCoor[0]*V3y;

    /* Compute the update of the parametric coordinates. */
    detInv = 1.0/(a00*a11 - a01*a10);
    const su2double dr = detInv*(f0*a11 - f1*a01);
    const su2double ds = detInv*(f1*a00 - f0*a10);

    /* Update the parametric coordinates. Note that the negative of the
       Jacobian is used, so the update must be added. */
    parCoor[0] += dr;
    parCoor[1] += ds;

    /* Check for convergence. */
    if(fabs(dr) <= tolNewton && fabs(ds) <= tolNewton) break;
  }

  /* Terminate if the Newton algorithm did not converge. */
  if(itCount == maxIt)
    Terminate("CADTElemClass::CoorInQuadrilateral", __FILE__, __LINE__,
              "Newton did not converge");

  /* Check if the parametric coordinates are inside the quadrilateral.
     If not, something is seriously wrong, because the inside test has been
     done already with the triangles. */
  if(parCoor[0] < paramLowerBound || parCoor[0] > paramUpperBound ||
     parCoor[1] < paramLowerBound || parCoor[1] > paramUpperBound)
    Terminate("CADTElemClass::CoorInQuadrilateral", __FILE__, __LINE__,
              "Point not inside the quadrilateral.");

  /* Compute the interpolation weights. */
  const su2double omr = 0.5*(1.0-parCoor[0]), opr = 0.5*(1.0+parCoor[0]);
  const su2double oms = 0.5*(1.0-parCoor[1]), ops = 0.5*(1.0+parCoor[1]);

  weightsInterpol[0] = omr*oms;
  weightsInterpol[1] = opr*oms;
  weightsInterpol[2] = opr*ops;
  weightsInterpol[3] = omr*ops;

  /* Return true, because the search was successful. */
  return true;
}

bool CADTElemClass::CoorInTetrahedron(const unsigned long elemID,
                                      const su2double     *coor,
                                      su2double           *parCoor,
                                      su2double           *weightsInterpol) {

  /* Determine the indices of the four vertices of the tetrahedron,
     multiplied by nDim (which is 3). This gives the position in the
     coordinate array where the coordinates of these points are stored. */
  unsigned long i0 = nDOFsPerElem[elemID];
  unsigned long i1 = i0 + 1, i2 = i0 + 2, i3 = i0+3;
  i0 = nDim*elemConns[i0]; i1 = nDim*elemConns[i1];
  i2 = nDim*elemConns[i2]; i3 = nDim*elemConns[i3];

  /* Determine the coordinates relative to the vertex 0. */
  const su2double xc = coor[0] - coorPoints[i0];
  const su2double yc = coor[1] - coorPoints[i0+1];
  const su2double zc = coor[2] - coorPoints[i0+2];

  const su2double x1 = coorPoints[i1]   - coorPoints[i0];
  const su2double y1 = coorPoints[i1+1] - coorPoints[i0+1];
  const su2double z1 = coorPoints[i1+2] - coorPoints[i0+2];

  const su2double x2 = coorPoints[i2]   - coorPoints[i0];
  const su2double y2 = coorPoints[i2+1] - coorPoints[i0+1];
  const su2double z2 = coorPoints[i2+2] - coorPoints[i0+2];

  const su2double x3 = coorPoints[i3]   - coorPoints[i0];
  const su2double y3 = coorPoints[i3+1] - coorPoints[i0+1];
  const su2double z3 = coorPoints[i3+2] - coorPoints[i0+2];

  /* The tetrahedron is parametrized by
     X-X0 = (r+1)*(X1-X0)/2 + (s+1)*(X2-X0)/2 + (t+1)*(X3-X0)/2,
     r, s, t >= -1, r+s+t <= -1. As a consequence, the parametric coordinates
     r, s and t can be solved easily. Note that X0 is 0 in the above expression,
     because the coordinates are relative to node 0. */
  const su2double detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  parCoor[0] =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  parCoor[1] = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  parCoor[2] =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* Check if the point resides within the tetrahedron and compute the
     interpolation weights if it is. */
  bool coorIsInside = false;
  if((parCoor[0] >= paramLowerBound) && (parCoor[1] >= paramLowerBound) &&
     (parCoor[2] >= paramLowerBound) &&
     ((parCoor[0]+parCoor[1]+parCoor[2]) <= paramLowerBound)) {
    coorIsInside = true;

    parCoor[0] = -0.5*(parCoor[0] + parCoor[1] + parCoor[2] + 1.0);
    parCoor[1] =  0.5*(parCoor[0] + 1.0);
    parCoor[2] =  0.5*(parCoor[1] + 1.0);
    parCoor[3] =  0.5*(parCoor[2] + 1.0);
  }

  /* Return the value of coorIsInside. */
  return coorIsInside;
}

bool CADTElemClass::CoorInPyramid(const unsigned long elemID,
                                  const su2double     *coor,
                                  su2double           *parCoor,
                                  su2double           *weightsInterpol) {

  /* Definition of the maximum number of iterations in the Newton solver
     and the tolerance level. */
  const unsigned short maxIt = 50;
  const su2double tolNewton  = 1.e-10;

  /* Determine the indices of the five vertices of the pyramid,
     multiplied by nDim (which is 3). This gives the position in the
     coordinate array where the coordinates of these points are stored. */
  unsigned long i0 = nDOFsPerElem[elemID];
  unsigned long i1 = i0 + 1, i2 = i0 + 2, i3 = i0+3, i4 = i0+4;
  i0 = nDim*elemConns[i0]; i1 = nDim*elemConns[i1];
  i2 = nDim*elemConns[i2]; i3 = nDim*elemConns[i3];
  i4 = nDim*elemConns[i4];

  /* Determine the coordinates relative to the vertex 0. */
  su2double xRel[5][3], xc[3];
  xc[0] = coor[0] - coorPoints[i0];
  xc[1] = coor[1] - coorPoints[i0+1];
  xc[2] = coor[2] - coorPoints[i0+2];

  xRel[0][0] = xRel[0][1] = xRel[0][2] = 0.0;

  xRel[1][0] = coorPoints[i1]   - coorPoints[i0];
  xRel[1][1] = coorPoints[i1+1] - coorPoints[i0+1];
  xRel[1][2] = coorPoints[i1+2] - coorPoints[i0+2];

  xRel[2][0] = coorPoints[i2]   - coorPoints[i0];
  xRel[2][1] = coorPoints[i2+1] - coorPoints[i0+1];
  xRel[2][2] = coorPoints[i2+2] - coorPoints[i0+2];

  xRel[3][0] = coorPoints[i3]   - coorPoints[i0];
  xRel[3][1] = coorPoints[i3+1] - coorPoints[i0+1];
  xRel[3][2] = coorPoints[i3+2] - coorPoints[i0+2];

  xRel[4][0] = coorPoints[i4]   - coorPoints[i0];
  xRel[4][1] = coorPoints[i4+1] - coorPoints[i0+1];
  xRel[4][2] = coorPoints[i4+2] - coorPoints[i0+2];

  /* Obtain an initial guess of the parametric coordinates by splitting the
     pyramid into tetrahedra. If this approach is not successful, this means
     that the point is not inside the true pyramid and false can be returned. */
  if( !InitialGuessContainmentPyramid(xc, xRel, parCoor) ) return false;

  /* The pyramid is parametrized by X-X0 = (Xi-X0)*li, where the sum runs over
     i = 0..4, although i = 0 does not give a contribution. The Lagrangian
     interpolation functions for the pyramid are given by:
     l0 = (1-t-2*r)*(1-t-2*s)/(8*(1-t)), l1 = (1-t+2*r)*(1-t-2*s)/(8*(1-t))
     l2 = (1-t+2*r)*(1-t+2*s)/(8*(1-t)), l3 = (1-t-2*r)*(1-t+2*s)/(8*(1-t)),
     l4 = (1+t)/2.
     The boundaries are -1 <= t <= 1, (t-1)/2 <= r,s <= (1-t)/2.
     As all coordinates are taken relative to vertex 0, X0 drops out and the
     nonlinear set of equations can be written as
     V0 - V1*r - V2*s - V3*t - V4*r*s/(1-t)  = 0, where
     V0 = xc - (4*x4+x1+x2+x3)/8, V1 = (x1+x2-x3)/4, V2 = (x2+x3-X1)/4,
     V3 = (4*X4-x1-x2-x3)/8, V4 = (x2-x1-x3)/2. First construct these vectors. */
  const su2double V0x = xc[0] - 0.5*xRel[4][0] - 0.125*(xRel[1][0]+xRel[2][0]+xRel[3][0]);
  const su2double V0y = xc[1] - 0.5*xRel[4][1] - 0.125*(xRel[1][1]+xRel[2][1]+xRel[3][1]);
  const su2double V0z = xc[2] - 0.5*xRel[4][2] - 0.125*(xRel[1][2]+xRel[2][2]+xRel[3][2]);

  const su2double V1x = 0.25*(xRel[1][0]+xRel[2][0]-xRel[3][0]);
  const su2double V1y = 0.25*(xRel[1][1]+xRel[2][1]-xRel[3][1]);
  const su2double V1z = 0.25*(xRel[1][2]+xRel[2][2]-xRel[3][2]);

  const su2double V2x = 0.25*(xRel[2][0]+xRel[3][0]-xRel[1][0]);
  const su2double V2y = 0.25*(xRel[2][1]+xRel[3][1]-xRel[1][1]);
  const su2double V2z = 0.25*(xRel[2][2]+xRel[3][2]-xRel[1][2]);

  const su2double V3x = 0.5*xRel[4][0] - 0.125*(xRel[1][0]+xRel[2][0]+xRel[3][0]);
  const su2double V3y = 0.5*xRel[4][1] - 0.125*(xRel[1][1]+xRel[2][1]+xRel[3][1]);
  const su2double V3z = 0.5*xRel[4][2] - 0.125*(xRel[1][2]+xRel[2][2]+xRel[3][2]);

  const su2double V4x = 0.5*(xRel[2][0]-xRel[1][0]-xRel[3][0]);
  const su2double V4y = 0.5*(xRel[2][1]-xRel[1][1]-xRel[3][1]);
  const su2double V4z = 0.5*(xRel[2][2]-xRel[1][2]-xRel[3][2]);

  /* Loop over the maximum number of iterations. */
  unsigned short itCount;
  for(itCount=0; itCount<maxIt; ++itCount) {

    /* Compute the values of the nonlinear functions. */
    su2double oneMinT = 1.0 - parCoor[2];
    if(fabs(oneMinT) < 1.e-10) {
      if(oneMinT < 0.0) oneMinT = -1.e-10;
      else              oneMinT =  1.e-10;
    }
    const su2double oneMinTInv = 1.0/oneMinT;

    const su2double f0 = V0x - parCoor[0]*V1x - parCoor[1]*V2x - parCoor[2]*V3x
                       - parCoor[0]*parCoor[1]*V4x*oneMinTInv;
    const su2double f1 = V0y - parCoor[0]*V1y - parCoor[1]*V2y - parCoor[2]*V3y
                       - parCoor[0]*parCoor[1]*V4y*oneMinTInv;
    const su2double f2 = V0z - parCoor[0]*V1z - parCoor[1]*V2z - parCoor[2]*V3z
                       - parCoor[0]*parCoor[1]*V4z*oneMinTInv;

    /* Compute the negative of the Jacobian matrix. */
    const su2double a00 = V1x + parCoor[1]*V4x*oneMinTInv;
    const su2double a01 = V2x + parCoor[0]*V4x*oneMinTInv;
    const su2double a02 = V3x + parCoor[0]*parCoor[1]*V4x*oneMinTInv*oneMinTInv;

    const su2double a10 = V1y + parCoor[1]*V4y*oneMinTInv;
    const su2double a11 = V2y + parCoor[0]*V4y*oneMinTInv;
    const su2double a12 = V3y + parCoor[0]*parCoor[1]*V4y*oneMinTInv*oneMinTInv;

    const su2double a20 = V1z + parCoor[1]*V4z*oneMinTInv;
    const su2double a21 = V2z + parCoor[0]*V4z*oneMinTInv;
    const su2double a22 = V3z + parCoor[0]*parCoor[1]*V4z*oneMinTInv*oneMinTInv;

    /* Compute the update of the parametric coordinates. */
    const su2double detInv = 1.0/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22
                           +      a01*a12*a20 + a02*a10*a21 - a02*a11*a20);
    const su2double dr =  detInv*(a01*a12*f2 - a01*a22*f1 - a02*a11*f2
                       +          a02*a21*f1 + a11*a22*f0 - a12*a21*f0);
    const su2double ds = -detInv*(a00*a12*f2 - a00*a22*f1 - a02*a10*f2
                       +          a02*a20*f1 + a10*a22*f0 - a12*a20*f0);
    const su2double dt =  detInv*(a00*a11*f2 - a00*a21*f1 - a01*a10*f2
                       +          a01*a20*f1 + a10*a21*f0 - a11*a20*f0);

    /* Update the parametric coordinates. Note that the negative of the
       Jacobian is used, so the update must be added. */
    parCoor[0] += dr;
    parCoor[1] += ds;
    parCoor[2] += dt;

    /* Check for convergence. */
    if(fabs(dr) <= tolNewton && fabs(ds) <= tolNewton && fabs(dt) <= tolNewton) break;
  }

  /* Terminate if the Newton algorithm did not converge. */
  if(itCount == maxIt)
    Terminate("CADTElemClass::CoorInQuadrilateral", __FILE__, __LINE__,
              "Newton did not converge");

  /* Check if the point resides within the pyramid and compute the
     interpolation weights if it is. */
  bool coorIsInside = false;
  if((parCoor[2] >= paramLowerBound) && (parCoor[2] <= paramUpperBound)) {
    const su2double lowRSBound = 0.5*(parCoor[2]-1.0) - tolInsideElem;
    const su2double uppRSBound = -lowRSBound;

    if((parCoor[0] >= lowRSBound) && (parCoor[0] <= uppRSBound) &&
       (parCoor[1] >= lowRSBound) && (parCoor[1] <= uppRSBound)) {
      coorIsInside = true;

      su2double oneMinT = 1.0 - parCoor[2];
      if(fabs(oneMinT) < 1.e-10) {
        if(oneMinT < 0.0) oneMinT = -1.e-10;
        else              oneMinT =  1.e-10;
      }
      const su2double oneMinTInv = 1.0/oneMinT;

      const su2double omr = (1.0-parCoor[2]-2.0*parCoor[0]);
      const su2double opr = (1.0-parCoor[2]+2.0*parCoor[0]);
      const su2double oms = (1.0-parCoor[2]-2.0*parCoor[1]);
      const su2double ops = (1.0-parCoor[2]+2.0*parCoor[1]);

      weightsInterpol[0] = 0.125*oneMinTInv*omr*oms;
      weightsInterpol[1] = 0.125*oneMinTInv*opr*oms;
      weightsInterpol[2] = 0.125*oneMinTInv*opr*ops;
      weightsInterpol[3] = 0.125*oneMinTInv*omr*ops;
      weightsInterpol[4] = 0.5*(1.0+parCoor[2]);
    }
  }

  /* Return the value of coorIsInside. */
  return coorIsInside;
}

bool CADTElemClass::InitialGuessContainmentPyramid(const su2double xRelC[3],
                                                   const su2double xRel[5][3],
                                                   su2double       *parCoor) {
  /* Tetrahedron, 0-1-3-4.
     Create the coordinates of the tetrahedron and of the point. */
  su2double x1 = xRel[1][0], y1 = xRel[1][1], z1 = xRel[1][2];
  su2double x2 = xRel[3][0], y2 = xRel[3][1], z2 = xRel[3][2];
  su2double x3 = xRel[4][0], y3 = xRel[4][1], z3 = xRel[4][2];

  su2double xc = xRelC[0], yc = xRelC[1], zc = xRelC[2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  su2double detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  parCoor[0] =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  parCoor[1] = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  parCoor[2] =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* If the point is inside this tetrahedron, return true. */
  if((parCoor[0] >= paramLowerBound) && (parCoor[1] >= paramLowerBound) &&
     (parCoor[2] >= paramLowerBound) &&
     ((parCoor[0]+parCoor[1]+parCoor[2]) <= paramLowerBound)) return true;

  /* Tetrahedron, 2-3-1-4.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[3][0]-xRel[2][0]; y1 = xRel[3][1]-xRel[2][1]; z1 = xRel[3][2]-xRel[2][2];
  x2 = xRel[1][0]-xRel[2][0]; y2 = xRel[1][1]-xRel[2][1]; z2 = xRel[1][2]-xRel[2][2];
  x3 = xRel[4][0]-xRel[2][0]; y3 = xRel[4][1]-xRel[2][1]; z3 = xRel[4][2]-xRel[2][2];

  xc = xRelC[0]-xRel[2][0]; yc = xRelC[1]-xRel[2][1]; zc = xRelC[2]-xRel[2][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  parCoor[0] =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  parCoor[1] = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  parCoor[2] =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* If the point is inside this tetrahedron, adapt the parametric coordinates
     to the real pyramid and return true. */
  if((parCoor[0] >= paramLowerBound) && (parCoor[1] >= paramLowerBound) &&
     (parCoor[2] >= paramLowerBound) &&
     ((parCoor[0]+parCoor[1]+parCoor[2]) <= paramLowerBound)) {
    parCoor[0] = 1.0 - parCoor[0];
    parCoor[1] = 1.0 - parCoor[1];
    return true;
  }

  /* Tetrahedron, 1-2-0-4.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[2][0]-xRel[1][0]; y1 = xRel[2][1]-xRel[1][1]; z1 = xRel[2][2]-xRel[1][2];
  x2 = xRel[0][0]-xRel[1][0]; y2 = xRel[0][1]-xRel[1][1]; z2 = xRel[0][2]-xRel[1][2];
  x3 = xRel[4][0]-xRel[1][0]; y3 = xRel[4][1]-xRel[1][1]; z3 = xRel[4][2]-xRel[1][2];

  xc = xRelC[0]-xRel[1][0]; yc = xRelC[1]-xRel[1][1]; zc = xRelC[2]-xRel[1][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  parCoor[0] =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  parCoor[1] = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  parCoor[2] =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* If the point is inside this tetrahedron, adapt the parametric coordinates
     to the real pyramid and return true. */
  if((parCoor[0] >= paramLowerBound) && (parCoor[1] >= paramLowerBound) &&
     (parCoor[2] >= paramLowerBound) &&
     ((parCoor[0]+parCoor[1]+parCoor[2]) <= paramLowerBound)) {
    const su2double r = parCoor[0];
    parCoor[0] = 1.0 - parCoor[1];
    parCoor[1] = r;
    return true;
  }

  /* Tetrahedron, 3-0-2-4.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[0][0]-xRel[3][0]; y1 = xRel[0][1]-xRel[3][1]; z1 = xRel[0][2]-xRel[3][2];
  x2 = xRel[2][0]-xRel[3][0]; y2 = xRel[2][1]-xRel[3][1]; z2 = xRel[2][2]-xRel[3][2];
  x3 = xRel[4][0]-xRel[3][0]; y3 = xRel[4][1]-xRel[3][1]; z3 = xRel[4][2]-xRel[3][2];

  xc = xRelC[0]-xRel[3][0]; yc = xRelC[1]-xRel[3][1]; zc = xRelC[2]-xRel[3][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  parCoor[0] =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  parCoor[1] = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  parCoor[2] =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* If the point is inside this tetrahedron, adapt the parametric coordinates
     to the real pyramid and return true. */
  if((parCoor[0] >= paramLowerBound) && (parCoor[1] >= paramLowerBound) &&
     (parCoor[2] >= paramLowerBound) &&
     ((parCoor[0]+parCoor[1]+parCoor[2]) <= paramLowerBound)) {
    const su2double r = parCoor[0];
    parCoor[0] = parCoor[1];
    parCoor[1] = 1.0 - r;
    return true;
  }

  /* The coordinate is in none of the four sub-tetrahedra of the pyramid. This
     implies that the point is not inside the true pyramid either and hence
     false is returned. */
  return false;
}

bool CADTElemClass::CoorInPrism(const unsigned long elemID,
                                const su2double     *coor,
                                su2double           *parCoor,
                                su2double           *weightsInterpol) {

  /* Definition of the maximum number of iterations in the Newton solver
     and the tolerance level. */
  const unsigned short maxIt = 50;
  const su2double tolNewton  = 1.e-10;

  /* Determine the indices of the six vertices of the prism,
     multiplied by nDim (which is 3). This gives the position in the
     coordinate array where the coordinates of these points are stored. */
  unsigned long i0 = nDOFsPerElem[elemID];
  unsigned long i1 = i0 + 1, i2 = i0 + 2, i3 = i0+3, i4 = i0+4, i5 = i0+5;
  i0 = nDim*elemConns[i0]; i1 = nDim*elemConns[i1];
  i2 = nDim*elemConns[i2]; i3 = nDim*elemConns[i3];
  i4 = nDim*elemConns[i4]; i5 = nDim*elemConns[i5];

  /* Determine the coordinates relative to the vertex 0. */
  su2double xRel[6][3], xc[3];
  xc[0] = coor[0] - coorPoints[i0];
  xc[1] = coor[1] - coorPoints[i0+1];
  xc[2] = coor[2] - coorPoints[i0+2];

  xRel[0][0] = xRel[0][1] = xRel[0][2] = 0.0;

  xRel[1][0] = coorPoints[i1]   - coorPoints[i0];
  xRel[1][1] = coorPoints[i1+1] - coorPoints[i0+1];
  xRel[1][2] = coorPoints[i1+2] - coorPoints[i0+2];

  xRel[2][0] = coorPoints[i2]   - coorPoints[i0];
  xRel[2][1] = coorPoints[i2+1] - coorPoints[i0+1];
  xRel[2][2] = coorPoints[i2+2] - coorPoints[i0+2];

  xRel[3][0] = coorPoints[i3]   - coorPoints[i0];
  xRel[3][1] = coorPoints[i3+1] - coorPoints[i0+1];
  xRel[3][2] = coorPoints[i3+2] - coorPoints[i0+2];

  xRel[4][0] = coorPoints[i4]   - coorPoints[i0];
  xRel[4][1] = coorPoints[i4+1] - coorPoints[i0+1];
  xRel[4][2] = coorPoints[i4+2] - coorPoints[i0+2];

  xRel[5][0] = coorPoints[i5]   - coorPoints[i0];
  xRel[5][1] = coorPoints[i5+1] - coorPoints[i0+1];
  xRel[5][2] = coorPoints[i5+2] - coorPoints[i0+2];

  /* Obtain an initial guess of the parametric coordinates by splitting the
     prism into tetrahedra. If this approach is not successful, this means
     that the point is not inside the true prism and false can be returned. */
  if( !InitialGuessContainmentPrism(xc, xRel, parCoor) ) return false;

  /* The prism is parametrized by X-X0 = (Xi-X0)*li, where the sum runs over
     i = 0..5, although i = 0 does not give a contribution. The Lagrangian
     interpolation functions for the prism are given by:
     l0 = -(r+s)*(1-t)/4, l1 = (r+1)*(1-t)/4, l2 = (s+1)*(1-t)/4,
     l3 = -(r+s)*(1+t)/4, l4 = (r+1)*(1+t)/4, l5 = (s+1)*(1+t)/4.
     The boundaries are r,s >= -1, r+s <= 0, -1 <= t <= 1.
     As all coordinates are taken relative to vertex 0, X0 drops out and the
     nonlinear set of equations can be written as
     V0 - V1*r - V2*s - V3*t - V4*r*t - V5*s*t  = 0, where
     V0 = xc - (x1+x2+x4+x5)/4, V1 = (x1+x4-x3)/4, V2 = (x2+x5-x3)/4,
     V3 = (x4+x5-x1-x2)/4, V4 = (x4-x1-x3)/4, V5 = (x5-x2-x3)/4.
     First construct these vectors. */
  const su2double V0x = xc[0] - 0.25*(xRel[1][0]+xRel[2][0]+xRel[4][0]+xRel[5][0]);
  const su2double V0y = xc[1] - 0.25*(xRel[1][1]+xRel[2][1]+xRel[4][1]+xRel[5][1]);
  const su2double V0z = xc[2] - 0.25*(xRel[1][2]+xRel[2][2]+xRel[4][2]+xRel[5][2]);

  const su2double V1x = 0.25*(xRel[1][0]+xRel[4][0]-xRel[3][0]);
  const su2double V1y = 0.25*(xRel[1][1]+xRel[4][1]-xRel[3][1]);
  const su2double V1z = 0.25*(xRel[1][2]+xRel[4][2]-xRel[3][2]);

  const su2double V2x = 0.25*(xRel[2][0]+xRel[5][0]-xRel[3][0]);
  const su2double V2y = 0.25*(xRel[2][1]+xRel[5][1]-xRel[3][1]);
  const su2double V2z = 0.25*(xRel[2][2]+xRel[5][2]-xRel[3][2]);

  const su2double V3x = 0.25*(xRel[4][0]+xRel[5][0]-xRel[1][0]-xRel[2][0]);
  const su2double V3y = 0.25*(xRel[4][1]+xRel[5][1]-xRel[1][1]-xRel[2][1]);
  const su2double V3z = 0.25*(xRel[4][2]+xRel[5][2]-xRel[1][2]-xRel[2][2]);

  const su2double V4x = 0.25*(xRel[4][0]-xRel[1][0]-xRel[3][0]);
  const su2double V4y = 0.25*(xRel[4][1]-xRel[1][1]-xRel[3][1]);
  const su2double V4z = 0.25*(xRel[4][2]-xRel[1][2]-xRel[3][2]);

  const su2double V5x = 0.25*(xRel[5][0]-xRel[2][0]-xRel[3][0]);
  const su2double V5y = 0.25*(xRel[5][1]-xRel[2][1]-xRel[3][1]);
  const su2double V5z = 0.25*(xRel[5][2]-xRel[2][2]-xRel[3][2]);

  /* Loop over the maximum number of iterations. */
  unsigned short itCount;
  for(itCount=0; itCount<maxIt; ++itCount) {

    const su2double f0 = V0x - parCoor[0]*V1x - parCoor[1]*V2x - parCoor[2]*V3x
                       - parCoor[0]*parCoor[2]*V4x - parCoor[1]*parCoor[2]*V5x;
    const su2double f1 = V0y - parCoor[0]*V1y - parCoor[1]*V2y - parCoor[2]*V3y
                       - parCoor[0]*parCoor[2]*V4y - parCoor[1]*parCoor[2]*V5y;
    const su2double f2 = V0z - parCoor[0]*V1z - parCoor[1]*V2z - parCoor[2]*V3z
                       - parCoor[0]*parCoor[2]*V4z - parCoor[1]*parCoor[2]*V5z;

    /* Compute the negative of the Jacobian matrix. */
    const su2double a00 = V1x + parCoor[2]*V4x;
    const su2double a01 = V2x + parCoor[2]*V5x;
    const su2double a02 = V3x + parCoor[0]*V4x + parCoor[1]*V5x;

    const su2double a10 = V1y + parCoor[2]*V4y;
    const su2double a11 = V2y + parCoor[2]*V5y;
    const su2double a12 = V3y + parCoor[0]*V4y + parCoor[1]*V5y;

    const su2double a20 = V1z + parCoor[2]*V4z;
    const su2double a21 = V2z + parCoor[2]*V5z;
    const su2double a22 = V3z + parCoor[0]*V4z + parCoor[1]*V5z;

    /* Compute the update of the parametric coordinates. */
    const su2double detInv = 1.0/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22
                           +      a01*a12*a20 + a02*a10*a21 - a02*a11*a20);
    const su2double dr =  detInv*(a01*a12*f2 - a01*a22*f1 - a02*a11*f2
                       +          a02*a21*f1 + a11*a22*f0 - a12*a21*f0);
    const su2double ds = -detInv*(a00*a12*f2 - a00*a22*f1 - a02*a10*f2
                       +          a02*a20*f1 + a10*a22*f0 - a12*a20*f0);
    const su2double dt =  detInv*(a00*a11*f2 - a00*a21*f1 - a01*a10*f2
                       +          a01*a20*f1 + a10*a21*f0 - a11*a20*f0);

    /* Update the parametric coordinates. Note that the negative of the
       Jacobian is used, so the update must be added. */
    parCoor[0] += dr;
    parCoor[1] += ds;
    parCoor[2] += dt;

    /* Check for convergence. */
    if(fabs(dr) <= tolNewton && fabs(ds) <= tolNewton && fabs(dt) <= tolNewton) break;
  }

  /* Terminate if the Newton algorithm did not converge. */
  if(itCount == maxIt)
    Terminate("CADTElemClass::CoorInPrism", __FILE__, __LINE__,
              "Newton did not converge");

  /* Check if the point resides within the prism and compute the
     interpolation weights if it is. */
  bool coorIsInside = false;
  if((parCoor[0] >= paramLowerBound) && (parCoor[1] >= paramLowerBound) &&
     ((parCoor[0]+parCoor[1]) <= tolInsideElem) &&
     (parCoor[2] >= paramLowerBound) && (parCoor[2] <= paramUpperBound)) {
    coorIsInside = true;

    const su2double omt = 0.25*(1.0-parCoor[2]), opt = 0.25*(1.0+parCoor[2]);

    weightsInterpol[0] = -omt*(parCoor[0]+parCoor[1]);
    weightsInterpol[1] =  omt*(parCoor[0]+1.0);
    weightsInterpol[2] =  omt*(parCoor[1]+1.0);
    weightsInterpol[3] = -opt*(parCoor[0]+parCoor[1]);
    weightsInterpol[4] =  opt*(parCoor[0]+1.0);
    weightsInterpol[5] =  opt*(parCoor[1]+1.0);
  }

  /* Return the value of coorIsInside. */
  return coorIsInside;
}

bool CADTElemClass::InitialGuessContainmentPrism(const su2double xRelC[3],
                                                 const su2double xRel[6][3],
                                                 su2double       *parCoor) {

  /* Tetrahedron, 0-1-2-3.
     Create the coordinates of the tetrahedron and of the point. */
  su2double x1 = xRel[1][0], y1 = xRel[1][1], z1 = xRel[1][2];
  su2double x2 = xRel[2][0], y2 = xRel[2][1], z2 = xRel[2][2];
  su2double x3 = xRel[3][0], y3 = xRel[3][1], z3 = xRel[3][2];

  su2double xc = xRelC[0], yc = xRelC[1], zc = xRelC[2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  su2double detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  su2double r =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  su2double s = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  su2double t =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real prism and return true. */
  if((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) &&
     ((r+s+t) <= paramLowerBound)) {
    parCoor[0] = r; parCoor[1] = s; parCoor[2] = t;
    return true;
  }

  /* Tetrahedron, 4-1-3-2.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[1][0]-xRel[4][0]; y1 = xRel[1][1]-xRel[4][1]; z1 = xRel[1][2]-xRel[4][2];
  x2 = xRel[3][0]-xRel[4][0]; y2 = xRel[3][1]-xRel[4][1]; z2 = xRel[3][2]-xRel[4][2];
  x3 = xRel[2][0]-xRel[4][0]; y3 = xRel[2][1]-xRel[4][1]; z3 = xRel[2][2]-xRel[4][2];

  xc = xRelC[0]-xRel[4][0]; yc = xRelC[1]-xRel[4][1]; zc = xRelC[2]-xRel[4][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  r =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  s = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  t =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real prism and return true. */
  if((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) &&
     ((r+s+t) <= paramLowerBound)) {
    parCoor[0] = 1.0 - s; parCoor[1] = t; parCoor[2] = 1.0 - r;
    return true;
  }

  /* Tetrahedron, 3-5-4-2.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[5][0]-xRel[3][0]; y1 = xRel[5][1]-xRel[3][1]; z1 = xRel[5][2]-xRel[3][2];
  x2 = xRel[4][0]-xRel[3][0]; y2 = xRel[4][1]-xRel[3][1]; z2 = xRel[4][2]-xRel[3][2];
  x3 = xRel[2][0]-xRel[3][0]; y3 = xRel[2][1]-xRel[3][1]; z3 = xRel[2][2]-xRel[3][2];

  xc = xRelC[0]-xRel[3][0]; yc = xRelC[1]-xRel[3][1]; zc = xRelC[2]-xRel[3][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  r =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  s = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  t =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real prism and return true. */
  if((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) &&
     ((r+s+t) <= paramLowerBound)) {
    parCoor[0] = s; parCoor[1] = r; parCoor[2] = 1.0 - t;
    return true;
  }

  /* Tetrahedron, 3-5-4-0.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[5][0]-xRel[3][0]; y1 = xRel[5][1]-xRel[3][1]; z1 = xRel[5][2]-xRel[3][2];
  x2 = xRel[4][0]-xRel[3][0]; y2 = xRel[4][1]-xRel[3][1]; z2 = xRel[4][2]-xRel[3][2];
  x3 = xRel[0][0]-xRel[3][0]; y3 = xRel[0][1]-xRel[3][1]; z3 = xRel[0][2]-xRel[3][2];

  xc = xRelC[0]-xRel[3][0]; yc = xRelC[1]-xRel[3][1]; zc = xRelC[2]-xRel[3][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  r =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  s = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  t =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real prism and return true. */
  if((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) &&
     ((r+s+t) <= paramLowerBound)) {
    parCoor[0] = s; parCoor[1] = r; parCoor[2] = 1.0 - t;
    return true;
  }

  /* Tetrahedron, 1-0-4-5.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[0][0]-xRel[1][0]; y1 = xRel[0][1]-xRel[1][1]; z1 = xRel[0][2]-xRel[1][2];
  x2 = xRel[4][0]-xRel[1][0]; y2 = xRel[4][1]-xRel[1][1]; z2 = xRel[4][2]-xRel[1][2];
  x3 = xRel[5][0]-xRel[1][0]; y3 = xRel[5][1]-xRel[1][1]; z3 = xRel[5][2]-xRel[1][2];

  xc = xRelC[0]-xRel[1][0]; yc = xRelC[1]-xRel[1][1]; zc = xRelC[2]-xRel[1][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  r =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  s = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  t =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real prism and return true. */
  if((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) &&
     ((r+s+t) <= paramLowerBound)) {
    parCoor[0] = 1.0 - r; parCoor[1] = t; parCoor[2] = s;
    return true;
  }

  /* Tetrahedron, 0-1-2-5.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[1][0]; y1 = xRel[1][1]; z1 = xRel[1][2];
  x2 = xRel[2][0]; y2 = xRel[2][1]; z2 = xRel[2][2];
  x3 = xRel[5][0]; y3 = xRel[5][1]; z3 = xRel[5][2];

  xc = xRelC[0]; yc = xRelC[1]; zc = xRelC[2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  r =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  s = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  t =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real prism and return true. */
  if((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) &&
     ((r+s+t) <= paramLowerBound)) {
    parCoor[0] = r; parCoor[1] = s; parCoor[2] = t;
    return true;
  }

  /* The coordinate is in none of the six sub-tetrahedra of the prism. This
     implies that the point is not inside the true prism either and hence
     false is returned. */
  return false;
}

bool CADTElemClass::CoorInHexahedron(const unsigned long elemID,
                                     const su2double     *coor,
                                     su2double           *parCoor,
                                     su2double           *weightsInterpol) {

  /* Definition of the maximum number of iterations in the Newton solver
     and the tolerance level. */
  const unsigned short maxIt = 50;
  const su2double tolNewton  = 1.e-10;

  /* Determine the indices of the eight vertices of the hexahedron,
     multiplied by nDim (which is 3). This gives the position in the
     coordinate array where the coordinates of these points are stored. */
  unsigned long i0 = nDOFsPerElem[elemID];
  unsigned long i1 = i0 + 1, i2 = i0 + 2, i3 = i0+3, i4 = i0+4, i5 = i0+5, i6 = i0+6, i7 = i0+7;
  i0 = nDim*elemConns[i0]; i1 = nDim*elemConns[i1];
  i2 = nDim*elemConns[i2]; i3 = nDim*elemConns[i3];
  i4 = nDim*elemConns[i4]; i5 = nDim*elemConns[i5];
  i6 = nDim*elemConns[i6]; i7 = nDim*elemConns[i7];

  /* Determine the coordinates relative to the vertex 0. */
  su2double xRel[8][3], xc[3];
  xc[0] = coor[0] - coorPoints[i0];
  xc[1] = coor[1] - coorPoints[i0+1];
  xc[2] = coor[2] - coorPoints[i0+2];

  xRel[0][0] = xRel[0][1] = xRel[0][2] = 0.0;

  xRel[1][0] = coorPoints[i1]   - coorPoints[i0];
  xRel[1][1] = coorPoints[i1+1] - coorPoints[i0+1];
  xRel[1][2] = coorPoints[i1+2] - coorPoints[i0+2];

  xRel[2][0] = coorPoints[i2]   - coorPoints[i0];
  xRel[2][1] = coorPoints[i2+1] - coorPoints[i0+1];
  xRel[2][2] = coorPoints[i2+2] - coorPoints[i0+2];

  xRel[3][0] = coorPoints[i3]   - coorPoints[i0];
  xRel[3][1] = coorPoints[i3+1] - coorPoints[i0+1];
  xRel[3][2] = coorPoints[i3+2] - coorPoints[i0+2];

  xRel[4][0] = coorPoints[i4]   - coorPoints[i0];
  xRel[4][1] = coorPoints[i4+1] - coorPoints[i0+1];
  xRel[4][2] = coorPoints[i4+2] - coorPoints[i0+2];

  xRel[5][0] = coorPoints[i5]   - coorPoints[i0];
  xRel[5][1] = coorPoints[i5+1] - coorPoints[i0+1];
  xRel[5][2] = coorPoints[i5+2] - coorPoints[i0+2];

  xRel[6][0] = coorPoints[i6]   - coorPoints[i0];
  xRel[6][1] = coorPoints[i6+1] - coorPoints[i0+1];
  xRel[6][2] = coorPoints[i6+2] - coorPoints[i0+2];

  xRel[7][0] = coorPoints[i7]   - coorPoints[i0];
  xRel[7][1] = coorPoints[i7+1] - coorPoints[i0+1];
  xRel[7][2] = coorPoints[i7+2] - coorPoints[i0+2];

  /* Obtain an initial guess of the parametric coordinates by splitting the
     hexahedron into tetrahedra. If this approach is not successful, this means
     that the point is not inside the true hexahedron and false can be returned. */
  if( !InitialGuessContainmentHexahedron(xc, xRel, parCoor) ) return false;

  /* The hexahedron is parametrized by X-X0 = (Xi-X0)*li, where the sum runs
     over i = 0..7, although i = 0 does not give a contribution. The Lagrangian
     interpolation functions for the hexahedron are given by:
     l0 = (1-r)*(1-s)*(1-t)/8, l1 = (1+r)*(1-s)*(1-t)/8,
     l2 = (1+r)*(1+s)*(1-t)/8, l3 = (1-r)*(1+s)*(1-t)/8,
     l4 = (1-r)*(1-s)*(1+t)/8, l5 = (1+r)*(1-s)*(1+t)/8,
     l6 = (1+r)*(1+s)*(1+t)/8, l7 = (1-r)*(1+s)*(1+t)/8.
     The boundaries are -1 <= r,s,t <= 1.
     As all coordinates are taken relative to vertex 0, X0 drops out and the
     nonlinear set of equations can be written as
     V0 - V1*r - V2*s - V3*t - V4*r*s - V5*r*t - V6*s*t - V7*r*s*t = 0, where
     V0 = xc - (x1+x2+x3+x4+x5+x6+x7)/8, V1 = (x1+x2-x3-x4+x5+x6-x7)/8,
     V2 =      (x2+x3-x1-x4-x5+x6+x7)/8, V3 = (x4+x5+x6+x7-x1-x2-x3)/8,
     V4 =      (x2+x4+x6-x1-x3-x5-x7)/8, V5 = (x3+x5+x6-x1-x2-x4-x7)/8,
     V6 =      (x1+x6+x7-x2-x3-x4-x5)/8, V7 = (x1+x3+x4+x6-x2-x5-x7)/8.
     First construct these vectors. */
  const su2double V0x = xc[0] - 0.125*(xRel[1][0]+xRel[2][0]+xRel[3][0]+xRel[4][0]+xRel[5][0]+xRel[6][0]+xRel[7][0]);
  const su2double V0y = xc[1] - 0.125*(xRel[1][1]+xRel[2][1]+xRel[3][1]+xRel[4][1]+xRel[5][1]+xRel[6][1]+xRel[7][1]);
  const su2double V0z = xc[2] - 0.125*(xRel[1][2]+xRel[2][2]+xRel[3][2]+xRel[4][2]+xRel[5][2]+xRel[6][2]+xRel[7][2]);

  const su2double V1x = 0.125*(xRel[1][0]+xRel[2][0]-xRel[3][0]-xRel[4][0]+xRel[5][0]+xRel[6][0]-xRel[7][0]);
  const su2double V1y = 0.125*(xRel[1][1]+xRel[2][1]-xRel[3][1]-xRel[4][1]+xRel[5][1]+xRel[6][1]-xRel[7][1]);
  const su2double V1z = 0.125*(xRel[1][2]+xRel[2][2]-xRel[3][2]-xRel[4][2]+xRel[5][2]+xRel[6][2]-xRel[7][2]);

  const su2double V2x = 0.125*(xRel[2][0]+xRel[3][0]-xRel[1][0]-xRel[4][0]-xRel[5][0]+xRel[6][0]+xRel[7][0]);
  const su2double V2y = 0.125*(xRel[2][1]+xRel[3][1]-xRel[1][1]-xRel[4][1]-xRel[5][1]+xRel[6][1]+xRel[7][1]);
  const su2double V2z = 0.125*(xRel[2][2]+xRel[3][2]-xRel[1][2]-xRel[4][2]-xRel[5][2]+xRel[6][2]+xRel[7][2]);

  const su2double V3x = 0.125*(xRel[4][0]+xRel[5][0]+xRel[6][0]+xRel[7][0]-xRel[1][0]-xRel[2][0]-xRel[3][0]);
  const su2double V3y = 0.125*(xRel[4][1]+xRel[5][1]+xRel[6][1]+xRel[7][1]-xRel[1][1]-xRel[2][1]-xRel[3][1]);
  const su2double V3z = 0.125*(xRel[4][2]+xRel[5][2]+xRel[6][2]+xRel[7][2]-xRel[1][2]-xRel[2][2]-xRel[3][2]);

  const su2double V4x = 0.125*(xRel[2][0]+xRel[4][0]+xRel[6][0]-xRel[1][0]-xRel[3][0]-xRel[5][0]-xRel[7][0]);
  const su2double V4y = 0.125*(xRel[2][1]+xRel[4][1]+xRel[6][1]-xRel[1][1]-xRel[3][1]-xRel[5][1]-xRel[7][1]);
  const su2double V4z = 0.125*(xRel[2][2]+xRel[4][2]+xRel[6][2]-xRel[1][2]-xRel[3][2]-xRel[5][2]-xRel[7][2]);

  const su2double V5x = 0.125*(xRel[3][0]+xRel[5][0]+xRel[6][0]-xRel[1][0]-xRel[2][0]-xRel[4][0]-xRel[7][0]);
  const su2double V5y = 0.125*(xRel[3][1]+xRel[5][1]+xRel[6][1]-xRel[1][1]-xRel[2][1]-xRel[4][1]-xRel[7][1]);
  const su2double V5z = 0.125*(xRel[3][2]+xRel[5][2]+xRel[6][2]-xRel[1][2]-xRel[2][2]-xRel[4][2]-xRel[7][2]);

  const su2double V6x = 0.125*(xRel[1][0]+xRel[6][0]+xRel[7][0]-xRel[2][0]-xRel[3][0]-xRel[4][0]-xRel[5][0]);
  const su2double V6y = 0.125*(xRel[1][1]+xRel[6][1]+xRel[7][1]-xRel[2][1]-xRel[3][1]-xRel[4][1]-xRel[5][1]);
  const su2double V6z = 0.125*(xRel[1][2]+xRel[6][2]+xRel[7][2]-xRel[2][2]-xRel[3][2]-xRel[4][2]-xRel[5][2]);

  const su2double V7x = 0.125*(xRel[1][0]+xRel[3][0]+xRel[4][0]+xRel[6][0]-xRel[2][0]-xRel[5][0]-xRel[7][0]);
  const su2double V7y = 0.125*(xRel[1][1]+xRel[3][1]+xRel[4][1]+xRel[6][1]-xRel[2][1]-xRel[5][1]-xRel[7][1]);
  const su2double V7z = 0.125*(xRel[1][2]+xRel[3][2]+xRel[4][2]+xRel[6][2]-xRel[2][2]-xRel[5][2]-xRel[7][2]);

  /* Loop over the maximum number of iterations. */
  unsigned short itCount;
  for(itCount=0; itCount<maxIt; ++itCount) {

    const su2double f0 = V0x - parCoor[0]*V1x - parCoor[1]*V2x - parCoor[2]*V3x
                       - parCoor[0]*parCoor[1]*V4x - parCoor[0]*parCoor[2]*V5x
                       - parCoor[1]*parCoor[2]*V6x - parCoor[0]*parCoor[1]*parCoor[2]*V7x;
    const su2double f1 = V0y - parCoor[0]*V1y - parCoor[1]*V2y - parCoor[2]*V3y
                       - parCoor[0]*parCoor[1]*V4y - parCoor[0]*parCoor[2]*V5y
                       - parCoor[1]*parCoor[2]*V6y - parCoor[0]*parCoor[1]*parCoor[2]*V7y;
    const su2double f2 = V0z - parCoor[0]*V1z - parCoor[1]*V2z - parCoor[2]*V3z
                       - parCoor[0]*parCoor[1]*V4z - parCoor[0]*parCoor[2]*V5z
                       - parCoor[1]*parCoor[2]*V6z - parCoor[0]*parCoor[1]*parCoor[2]*V7z;

    /* Compute the negative of the Jacobian matrix. */
    const su2double a00 = V1x + parCoor[1]*V4x + parCoor[2]*V5x + parCoor[1]*parCoor[2]*V7x;
    const su2double a01 = V2x + parCoor[0]*V4x + parCoor[2]*V6x + parCoor[0]*parCoor[2]*V7x;
    const su2double a02 = V3x + parCoor[0]*V5x + parCoor[1]*V6x + parCoor[0]*parCoor[1]*V7x;

    const su2double a10 = V1y + parCoor[1]*V4y + parCoor[2]*V5y + parCoor[1]*parCoor[2]*V7y;
    const su2double a11 = V2y + parCoor[0]*V4y + parCoor[2]*V6y + parCoor[0]*parCoor[2]*V7y;
    const su2double a12 = V3y + parCoor[0]*V5y + parCoor[1]*V6y + parCoor[0]*parCoor[1]*V7y;

    const su2double a20 = V1z + parCoor[1]*V4z + parCoor[2]*V5z + parCoor[1]*parCoor[2]*V7z;
    const su2double a21 = V2z + parCoor[0]*V4z + parCoor[2]*V6z + parCoor[0]*parCoor[2]*V7z;
    const su2double a22 = V3z + parCoor[0]*V5z + parCoor[1]*V6z + parCoor[0]*parCoor[1]*V7z;

    /* Compute the update of the parametric coordinates. */
    const su2double detInv = 1.0/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22
                           +      a01*a12*a20 + a02*a10*a21 - a02*a11*a20);
    const su2double dr =  detInv*(a01*a12*f2 - a01*a22*f1 - a02*a11*f2
                       +          a02*a21*f1 + a11*a22*f0 - a12*a21*f0);
    const su2double ds = -detInv*(a00*a12*f2 - a00*a22*f1 - a02*a10*f2
                       +          a02*a20*f1 + a10*a22*f0 - a12*a20*f0);
    const su2double dt =  detInv*(a00*a11*f2 - a00*a21*f1 - a01*a10*f2
                       +          a01*a20*f1 + a10*a21*f0 - a11*a20*f0);

    /* Update the parametric coordinates. Note that the negative of the
       Jacobian is used, so the update must be added. */
    parCoor[0] += dr;
    parCoor[1] += ds;
    parCoor[2] += dt;

    /* Check for convergence. */
    if(fabs(dr) <= tolNewton && fabs(ds) <= tolNewton && fabs(dt) <= tolNewton) break;
  }

  /* Terminate if the Newton algorithm did not converge. */
  if(itCount == maxIt)
    Terminate("CADTElemClass::CoorInHexahedron", __FILE__, __LINE__,
              "Newton did not converge");

  /* Check if the point resides within the hexcahedron and compute the
     interpolation weights if it is. */
  bool coorIsInside = false;
  if((parCoor[0] >= paramLowerBound) && (parCoor[0] <= paramUpperBound) &&
     (parCoor[1] >= paramLowerBound) && (parCoor[1] <= paramUpperBound) &&
     (parCoor[2] >= paramLowerBound) && (parCoor[2] <= paramUpperBound)) {
    coorIsInside = true;

    const su2double omr = 0.5*(1.0-parCoor[0]), opr = 0.5*(1.0+parCoor[0]);
    const su2double oms = 0.5*(1.0-parCoor[1]), ops = 0.5*(1.0+parCoor[1]);
    const su2double omt = 0.5*(1.0-parCoor[2]), opt = 0.5*(1.0+parCoor[2]);

    weightsInterpol[0] = omr*oms*omt;
    weightsInterpol[1] = opr*oms*omt;
    weightsInterpol[2] = opr*ops*omt;
    weightsInterpol[3] = omr*ops*omt;
    weightsInterpol[4] = omr*oms*opt;
    weightsInterpol[5] = opr*oms*opt;
    weightsInterpol[6] = opr*ops*opt;
    weightsInterpol[7] = omr*ops*opt;
  }

  /* Return the value of coorIsInside. */
  return coorIsInside;
}

bool CADTElemClass::InitialGuessContainmentHexahedron(const su2double xRelC[3],
                                                      const su2double xRel[8][3],
                                                      su2double       *parCoor) {
  /* Tetrahedron, 0-1-2-5.
     Create the coordinates of the tetrahedron and of the point. */
  su2double x1 = xRel[1][0], y1 = xRel[1][1], z1 = xRel[1][2];
  su2double x2 = xRel[2][0], y2 = xRel[2][1], z2 = xRel[2][2];
  su2double x3 = xRel[5][0], y3 = xRel[5][1], z3 = xRel[5][2];

  su2double xc = xRelC[0], yc = xRelC[1], zc = xRelC[2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  su2double detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  su2double r =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  su2double s = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  su2double t =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real hexahedron and return true. */
  if((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) &&
     ((r+s+t) <= paramLowerBound)) {
    parCoor[0] = r; parCoor[1] = s; parCoor[2] = t;
    return true;
  }

  /* Tetrahedron, 4-7-5-0.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[7][0]-xRel[4][0]; y1 = xRel[7][1]-xRel[4][1]; z1 = xRel[7][2]-xRel[4][2];
  x2 = xRel[5][0]-xRel[4][0]; y2 = xRel[5][1]-xRel[4][1]; z2 = xRel[5][2]-xRel[4][2];
  x3 = xRel[0][0]-xRel[4][0]; y3 = xRel[0][1]-xRel[4][1]; z3 = xRel[0][2]-xRel[4][2];

  xc = xRelC[0]-xRel[4][0]; yc = xRelC[1]-xRel[4][1]; zc = xRelC[2]-xRel[4][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  r =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  s = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  t =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real hexahedron and return true. */
  if((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) &&
     ((r+s+t) <= paramLowerBound)) {
    parCoor[0] = s; parCoor[1] = r; parCoor[2] = 1.0 - t;
    return true;
  }

  /* Tetrahedron, 6-7-5-2.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[7][0]-xRel[6][0]; y1 = xRel[7][1]-xRel[6][1]; z1 = xRel[7][2]-xRel[6][2];
  x2 = xRel[5][0]-xRel[6][0]; y2 = xRel[5][1]-xRel[6][1]; z2 = xRel[5][2]-xRel[6][2];
  x3 = xRel[2][0]-xRel[6][0]; y3 = xRel[2][1]-xRel[6][1]; z3 = xRel[2][2]-xRel[6][2];

  xc = xRelC[0]-xRel[6][0]; yc = xRelC[1]-xRel[6][1]; zc = xRelC[2]-xRel[6][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  r =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  s = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  t =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real hexahedron and return true. */
  if((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) &&
     ((r+s+t) <= paramLowerBound)) {
    parCoor[0] = 1.0 - s; parCoor[1] = 1.0 - r; parCoor[2] = 1.0 - t;
    return true;
  }

  /* Tetrahedron, 3-0-2-7.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[0][0]-xRel[3][0]; y1 = xRel[0][1]-xRel[3][1]; z1 = xRel[0][2]-xRel[3][2];
  x2 = xRel[2][0]-xRel[3][0]; y2 = xRel[2][1]-xRel[3][1]; z2 = xRel[2][2]-xRel[3][2];
  x3 = xRel[7][0]-xRel[3][0]; y3 = xRel[7][1]-xRel[3][1]; z3 = xRel[7][2]-xRel[3][2];

  xc = xRelC[0]-xRel[3][0]; yc = xRelC[1]-xRel[3][1]; zc = xRelC[2]-xRel[3][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  r =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  s = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  t =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real hexahedron and return true. */
  if((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) &&
     ((r+s+t) <= paramLowerBound)) {
    parCoor[0] = 1.0 - s; parCoor[1] = r; parCoor[2] = 1.0 - t;
    return true;
  }

  /* Tetrahedron, 0-5-2-7.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[5][0]; y1 = xRel[5][1]; z1 = xRel[5][2];
  x2 = xRel[2][0]; y2 = xRel[2][1]; z2 = xRel[2][2];
  x3 = xRel[7][0]; y3 = xRel[7][1]; z3 = xRel[7][2];

  xc = xRelC[0]; yc = xRelC[1]; zc = xRelC[2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  r =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  s = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  t =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real hexahedron and return true. */
  if((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) &&
     ((r+s+t) <= paramLowerBound)) {
    parCoor[0] = 1.0+r+s; parCoor[1] = 1.0+s+t; parCoor[2] = 1.0+r+t;
    return true;
  }

  /* Tetrahedron, 0-1-3-4.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[1][0]; y1 = xRel[1][1]; z1 = xRel[1][2];
  x2 = xRel[3][0]; y2 = xRel[3][1]; z2 = xRel[3][2];
  x3 = xRel[4][0]; y3 = xRel[4][1]; z3 = xRel[4][2];

  xc = xRelC[0]; yc = xRelC[1]; zc = xRelC[2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  r =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  s = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  t =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real hexahedron and return true. */
  if((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) &&
     ((r+s+t) <= paramLowerBound)) {
    parCoor[0] = r; parCoor[1] = s; parCoor[2] = t;
    return true;
  }

  /* Tetrahedron, 7-6-4-3.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[6][0]-xRel[7][0]; y1 = xRel[6][1]-xRel[7][1]; z1 = xRel[6][2]-xRel[7][2];
  x2 = xRel[4][0]-xRel[7][0]; y2 = xRel[4][1]-xRel[7][1]; z2 = xRel[4][2]-xRel[7][2];
  x3 = xRel[3][0]-xRel[7][0]; y3 = xRel[3][1]-xRel[7][1]; z3 = xRel[3][2]-xRel[7][2];

  xc = xRelC[0]-xRel[7][0]; yc = xRelC[1]-xRel[7][1]; zc = xRelC[2]-xRel[7][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  r =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  s = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  t =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real hexahedron and return true. */
  if((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) &&
     ((r+s+t) <= paramLowerBound)) {
    parCoor[0] = r; parCoor[1] = 1.0 - s; parCoor[2] = 1.0 - t;
    return true;
  }

  /* Tetrahedron, 5-4-6-1.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[4][0]-xRel[5][0]; y1 = xRel[4][1]-xRel[5][1]; z1 = xRel[4][2]-xRel[5][2];
  x2 = xRel[6][0]-xRel[5][0]; y2 = xRel[6][1]-xRel[5][1]; z2 = xRel[6][2]-xRel[5][2];
  x3 = xRel[1][0]-xRel[5][0]; y3 = xRel[1][1]-xRel[5][1]; z3 = xRel[1][2]-xRel[5][2];

  xc = xRelC[0]-xRel[5][0]; yc = xRelC[1]-xRel[5][1]; zc = xRelC[2]-xRel[5][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  r =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  s = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  t =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real hexahedron and return true. */
  if((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) &&
     ((r+s+t) <= paramLowerBound)) {
    parCoor[0] = 1.0 - r; parCoor[1] = s; parCoor[2] = 1.0 - t;
    return true;
  }

  /* Tetrahedron, 2-3-1-6.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[3][0]-xRel[2][0]; y1 = xRel[3][1]-xRel[2][1]; z1 = xRel[3][2]-xRel[2][2];
  x2 = xRel[1][0]-xRel[2][0]; y2 = xRel[1][1]-xRel[2][1]; z2 = xRel[1][2]-xRel[2][2];
  x3 = xRel[6][0]-xRel[2][0]; y3 = xRel[6][1]-xRel[2][1]; z3 = xRel[6][2]-xRel[2][2];

  xc = xRelC[0]-xRel[2][0]; yc = xRelC[1]-xRel[2][1]; zc = xRelC[2]-xRel[2][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  r =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  s = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  t =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real hexahedron and return true. */
  if((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) &&
     ((r+s+t) <= paramLowerBound)) {
    parCoor[0] = 1.0 - r; parCoor[1] = 1.0 - s; parCoor[2] = t;
    return true;
  }

  /* Tetrahedron, 3-4-1-6.
     Create the coordinates of the tetrahedron and of the point. */
  x1 = xRel[4][0]-xRel[3][0]; y1 = xRel[4][1]-xRel[3][1]; z1 = xRel[4][2]-xRel[3][2];
  x2 = xRel[1][0]-xRel[3][0]; y2 = xRel[1][1]-xRel[3][1]; z2 = xRel[1][2]-xRel[3][2];
  x3 = xRel[6][0]-xRel[3][0]; y3 = xRel[6][1]-xRel[3][1]; z3 = xRel[6][2]-xRel[3][2];

  xc = xRelC[0]-xRel[3][0]; yc = xRelC[1]-xRel[3][1]; zc = xRelC[2]-xRel[3][2];

  /* Determine the parametric coordinates inside this tetrahedron. */
  detInv = 2.0/(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);
  r =  detInv*(x2*y3*zc - x2*yc*z3 - x3*y2*zc + x3*yc*z2 + xc*y2*z3 - xc*y3*z2) - 1.0;
  s = -detInv*(x1*y3*zc - x1*yc*z3 - x3*y1*zc + x3*yc*z1 + xc*y1*z3 - xc*y3*z1) - 1.0;
  t =  detInv*(x1*y2*zc - x1*yc*z2 - x2*y1*zc + x2*yc*z1 + xc*y1*z2 - xc*y2*z1) - 1.0;

  /* If the point is inside this tetrahedron, set the parametric coordinates for
     the real hexahedron and return true. */
  if((r >= paramLowerBound) && (s >= paramLowerBound) && (t >= paramLowerBound) &&
     ((r+s+t) <= paramLowerBound)) {
    parCoor[0] = 1.0+s+t; parCoor[1] = -1.0-r-s; parCoor[2] = 1.0+r+t;
    return true;
  }

  /* The coordinate is in none of the ten sub-tetrahedra of the hexahedron.
     This implies that the point is not inside the true hexahedron either
     and hence false is returned. */
  return false;
}

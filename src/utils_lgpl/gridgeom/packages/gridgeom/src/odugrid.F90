!----- LGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2018.                                
!                                                                               
!  This library is free software; you can redistribute it and/or                
!  modify it under the terms of the GNU Lesser General Public                   
!  License as published by the Free Software Foundation version 2.1.                 
!                                                                               
!  This library is distributed in the hope that it will be useful,              
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
!  Lesser General Public License for more details.                              
!                                                                               
!  You should have received a copy of the GNU Lesser General Public             
!  License along with this library; if not, see <http://www.gnu.org/licenses/>. 
!                                                                               
!  contact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D" and "Deltares"    
!  are registered trademarks of Stichting Deltares, and remain the property of  
!  Stichting Deltares. All rights reserved.                                     
!                                                                               
!-------------------------------------------------------------------------------

! $Id: odugrid.F90 61898 2018-09-24 15:17:11Z carniato $
! $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/trunk/src/utils_lgpl/gridgeom/packages/gridgeom/src/odugrid.F90 $

! Module for grid operations on 1d networks
   
module odugrid
!use xxx

implicit none

!user defined data types
   contains
!
! Calculates x and y coordinates from UGrid filr format
! 
function odu_get_xy_coordinates(branchids, branchoffsets, geopointsX, geopointsY, nbranchgeometrynodes, branchlengths, jsferic, meshXCoords, meshYCoords) result(ierr)

   use geometry_module, only: sphertocart3D, cart3Dtospher
   use m_missing, only : dmiss
   
   integer, intent(in)               :: branchids(:), nbranchgeometrynodes(:)
   double precision, intent(in)      :: branchoffsets(:), geopointsX(:), geopointsY(:), branchlengths(:)
   double precision, intent(inout)   :: meshXCoords(:), meshYCoords(:)
   integer, intent(in)               :: jsferic

   integer                           :: angle, i, k, ierr, ind, branchid, nsegments
   double precision, allocatable     :: branchSegmentLengths(:)
   double precision, allocatable     :: xincrement(:), yincrement(:), zincrement(:)
   double precision, allocatable     :: deltaX(:), deltaY(:), deltaZ(:)
   double precision, allocatable     :: cartMeshXCoords(:), cartMeshYCoords(:), cartMeshZCoords(:)
   double precision, allocatable     :: cartGeopointsX(:), cartGeopointsY(:), cartGeopointsZ(:)
   double precision, allocatable     :: meshZCoords(:)
   double precision, allocatable     :: geopointsZ(:)  !returned by sphertocart3D 
   double precision                  :: totalLength, afac, fractionbranchlength, maxlat, previousLength
   integer                           :: nBranchSegments, nbranches, br
   integer                           :: startGeometryNode, endGeometryNode, nGeometrySegments
   integer                           :: startMeshNode, endMeshNode 
   integer, allocatable              :: meshnodemapping(:,:)
      
   ierr = 0
   ! the number of geometry segments is always equal to number of geopoints - 1
   nBranchSegments = size(geopointsX,1) - 1
   allocate(branchSegmentLengths(nBranchSegments))
   
   allocate(deltaX(size(geopointsX,1) - 1))
   allocate(deltaY(size(geopointsX,1) - 1))
   allocate(deltaZ(size(geopointsX,1) - 1))
   
   allocate(xincrement(size(geopointsX,1) - 1))
   allocate(yincrement(size(geopointsX,1) - 1))
   allocate(zincrement(size(geopointsX,1) - 1))
   
   allocate(cartGeopointsX(size(geopointsX,1)))
   allocate(cartGeopointsY(size(geopointsX,1)))
   allocate(cartGeopointsZ(size(geopointsX,1)))
   
   allocate(cartMeshXCoords(size(meshXCoords,1)))
   allocate(cartMeshYCoords(size(meshXCoords,1)))
   allocate(cartMeshZCoords(size(meshXCoords,1)))
   
   !If the coordinates are spherical, trasform them in cartesian, so we operate in linear space
   maxlat = 0
   if (jsferic == 1) then
      do i = 1, size(geopointsX,1)
         maxlat = max(maxlat,geopointsX(i))
         call sphertocart3D(geopointsX(i),geopointsY(i),cartGeopointsX(i),cartGeopointsY(i), cartGeopointsZ(i)) 
      end do
   else
      cartGeopointsX(:) = geopointsX(:)
      cartGeopointsY(:) = geopointsY(:)
      cartGeopointsZ(:) = 0
   endif
   
   !map the mesh nodes
   nbranches = size(branchlengths,1)
   allocate(meshnodemapping(2,nbranches))
   meshnodemapping = -1
   ierr = odu_get_start_end_nodes_of_branches(branchids, meshnodemapping(1,:), meshnodemapping(2,:))
   
   ! initialization
   startGeometryNode    = 1
   do br = 1, nbranches
      ! starting and ending nodes
      startMeshNode         =  meshnodemapping(1,br)
      endMeshNode           =  meshnodemapping(2,br)
      !number of geometry segments for the current branch
      nGeometrySegments     = nbranchgeometrynodes(br) - 1
      !ending geometry point
      endGeometryNode       = startGeometryNode + nGeometrySegments
      !calculate the branch lenghts
      totalLength = 0.0d0
      do i = startGeometryNode, endGeometryNode - 1
         deltaX(i) = cartGeopointsX(i+1) - cartGeopointsX(i)
         deltaY(i) = cartGeopointsY(i+1) - cartGeopointsY(i)
         deltaZ(i) = cartGeopointsZ(i+1) - cartGeopointsZ(i)
         branchSegmentLengths(i)= sqrt(deltaX(i)**2+deltaY(i)**2+deltaZ(i)**2)
         totalLength = totalLength + branchSegmentLengths(i)
      enddo
      !correct for total segment length
      if (totalLength > 1.0d-6) then
         afac = branchlengths(br)/totalLength
         branchSegmentLengths(startGeometryNode: endGeometryNode - 1) = branchSegmentLengths(startGeometryNode: endGeometryNode - 1) * afac
      end if
   
      !calculate the increments
      do i = startGeometryNode, endGeometryNode - 1
         if (branchSegmentLengths(i) > 1.0d-6) then
            xincrement(i)  = deltaX(i)/branchSegmentLengths(i)
            yincrement(i)  = deltaY(i)/branchSegmentLengths(i)
            zincrement(i)  = deltaZ(i)/branchSegmentLengths(i)
         else
            xincrement(i)  = 0.d0
            yincrement(i)  = 0.d0   
            zincrement(i)  = 0.d0   
         endif
      enddo
      !now loop over the mesh points
      ind            = startGeometryNode
      totallength    = 0.d0
      previousLength = 0.d0 
      do i = startMeshNode, endMeshNode         
         !determine max and min lengths
         totalLength = previousLength
         do k = ind, endGeometryNode - 1
            totalLength = totalLength + branchSegmentLengths(k)
            if (totalLength > branchoffsets(i)) then
                  previousLength = totalLength - branchSegmentLengths(k)
                  ind = k
               exit
            endif
         enddo
         fractionbranchlength =  branchoffsets(i) - previousLength
         cartMeshXCoords(i) = cartGeopointsX(ind) + fractionbranchlength * xincrement(ind)
         cartMeshYCoords(i) = cartGeopointsY(ind) + fractionbranchlength * yincrement(ind)
         !TODO: this function should also return meshZCoords (it is relevant if coordinates are spheric) 
         cartMeshZCoords(i) = cartGeopointsZ(ind) + fractionbranchlength * zincrement(ind)
      enddo
      !update geometry indexes
      startGeometryNode = endGeometryNode + 1
   enddo
   
   if (jsferic == 1) then
      do i = 1, size(meshXCoords,1)
         call cart3Dtospher(cartMeshXCoords(i),cartMeshYCoords(i),cartMeshZCoords(i),meshXCoords(i),meshYCoords(i), maxlat) 
      end do
   else
      meshXCoords(:) = cartMeshXCoords(:)
      meshYCoords(:) = cartMeshYCoords(:)
      ! meshZCoords(:) = cartMeshZCoords(:)  !  meshZCoords is never used and assignment causes runtime error
   endif

end function odu_get_xy_coordinates


!Calculate the start and the end nodes of the branches
function odu_get_start_end_nodes_of_branches(branchidx, branchStartNode, branchEndNode) result(ierr)

   integer, dimension(:), intent(in)      :: branchidx

   integer, dimension(:), intent(inout)   :: branchStartNode
   integer, dimension(:), intent(inout)   :: branchEndNode
   integer                                :: ierr, i, ibran, numnode, nbranches 

   ! Get the starting and endig indexes of the grid points
   ierr  =  0
   ibran =  0
   numnode = size(branchidx)
   nbranches = size(branchStartNode)
   do i = 1, numnode
      if (branchidx(i) > ibran) then
         ibran = branchidx(i)
         branchStartNode(ibran) = i
         if (i > 2 .and. ibran >= 2) then
            branchEndNode(ibran - 1) = i - 1
         endif
      endif
   enddo
   branchEndNode(nbranches) = numnode

end function odu_get_start_end_nodes_of_branches

end module odugrid
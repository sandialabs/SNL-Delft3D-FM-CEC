!----- AGPL --------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2017-2022.                                
!                                                                               
!  This file is part of Delft3D (D-Flow Flexible Mesh component).               
!                                                                               
!  Delft3D is free software: you can redistribute it and/or modify              
!  it under the terms of the GNU Affero General Public License as               
!  published by the Free Software Foundation version 3.                         
!                                                                               
!  Delft3D  is distributed in the hope that it will be useful,                  
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU Affero General Public License for more details.                          
!                                                                               
!  You should have received a copy of the GNU Affero General Public License     
!  along with Delft3D.  If not, see <http://www.gnu.org/licenses/>.             
!                                                                               
!  contact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D",                  
!  "D-Flow Flexible Mesh" and "Deltares" are registered trademarks of Stichting 
!  Deltares, and remain the property of Stichting Deltares. All rights reserved.
!                                                                               
!-------------------------------------------------------------------------------

! $Id: m_boat.f90 140618 2022-01-12 13:12:04Z klapwijk $
! $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/tags/delft3dfm/141476/src/engines_gpl/dflowfm/packages/dflowfm_kernel/src/dflowfm_data/m_boat.f90 $

module m_boat
  implicit none
  double precision, allocatable :: XBOAT (:), YBOAT(:), ZBOAT(:)  ! INLEZEN &
  double precision, allocatable :: XBOOT (:), YBOOT(:), ZBOOT(:)  ! AFBEELDEN LIJNEN
  integer, allocatable          :: NCBOAT(:)                      !
  integer                       :: MXBOAT, MAXBOAT, NCLBOAT
  double precision              :: BLEN = 25, BHEIGHT = 4, BWIDTH = 5, BHPMAX = 1000, BHPPERC = 0.5d0
  integer                       :: KKB(10) = 0
end module m_boat

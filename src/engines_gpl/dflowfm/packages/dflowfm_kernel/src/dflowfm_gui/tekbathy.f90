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

! $Id: tekbathy.f90 140618 2022-01-12 13:12:04Z klapwijk $
! $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/tags/delft3dfm/141476/src/engines_gpl/dflowfm/packages/dflowfm_kernel/src/dflowfm_gui/tekbathy.f90 $

 subroutine tekbathy(ja)
 use unstruc_display
 use m_flowgeom
 use m_flow
 use gridoperations
 implicit none
 integer :: nodemode, nodewhat,ndraw
 integer :: k, ja, nn, ncol
 double precision :: znod, zn
 common /drawthis/ ndraw(50)
 logical inview

 if (ndraw(39) == 0) return

 nodewhat  = ndraw(28)
 ndraw(28) = 3

 do k = 1,ndxi
    if (mod(k,200) == 0) then
       call halt(ja)
       if (ja == 1) then
          ndraw(28) = nodewhat
          return
       endif
    endif

    if (inview( xz(k), yz(k) ) ) then
       zn = znod(k)
       call isocol2(zn,ncol)
       nn = size( nd(k)%x )
       call PFILLER(nd(k)%x, nd(k)%y, nn,NCOL,NCol)
    endif
 enddo

 ndraw(28) = nodewhat
 end subroutine tekbathy
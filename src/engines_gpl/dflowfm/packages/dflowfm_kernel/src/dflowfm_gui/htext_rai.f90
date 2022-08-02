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

! $Id: htext_rai.f90 140618 2022-01-12 13:12:04Z klapwijk $
! $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/tags/delft3dfm/141476/src/engines_gpl/dflowfm/packages/dflowfm_kernel/src/dflowfm_gui/htext_rai.f90 $

 subroutine htext_rai(val,x,y,xx,zz,ihv)
 use m_raaitek
 implicit none
 double precision  :: val,x,y,xx,zz
 double precision  :: fx, fy, xa, ya
 integer           :: ihv
 fx = xs2m-xs1m
 fy = ys2m-ys1m
 if (ihv == 1) then
    xa = fx*(x-xx-xw1m)/(xw2m-xw1m) + xs1m
    ya = fy*(y   -yw1m)/(yw2m-yw1m) + ys1m
    call movabs(xa,ya)
    xa = fx*(x+xx-xw1m)/(xw2m-xw1m) + xs1m
    call lnabs (xa,ya)
    xa = fx*(x-11d0*xx-xw1m)/(xw2m-xw1m) + xs1m
 else if (ihv == 2) then
    xa = fx*(x   -xw1m)/(xw2m-xw1m) + xs1m
    ya = fy*(y-zz-yw1m)/(yw2m-yw1m) + ys1m
    call movabs(xa,ya)
    ya = fy*(y+zz-yw1m)/(yw2m-yw1m) + ys1m
    call lnabs (xa,ya)
    xa = fx*(x-5d0*xx-xw1m)/(xw2m-xw1m) + xs1m
    ya = fy*(y-3d0*zz-yw1m)/(yw2m-yw1m) + ys1m
 endif
 call htext(val,xa,ya)
 end subroutine htext_rai

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

! $Id: makelongculverts_commandline.f90 141035 2022-04-07 14:34:02Z buwalda $
! $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/tags/delft3dfm/141476/src/engines_gpl/dflowfm/packages/dflowfm_kernel/src/dflowfm_kernel/prepost/makelongculverts_commandline.f90 $

!>  perform partitioning from command line
subroutine makelongculverts_commandline()
   use unstruc_model
   use m_readstructures
   use string_module, only: strsplit
   use m_longculverts
   use unstruc_netcdf, only :  unc_write_net, UNC_CONV_UGRID
   use unstruc_channel_flow, only: network
   
   character(len=1024) :: fnamesstring
   character(len=:), allocatable :: converted_fnamesstring
   character(len=:), allocatable :: converted_crsdefsstring
   character(len=:), allocatable :: tempstring_crsdef
   character(len=:), allocatable :: tempstring_fnames
   character(len=200), dimension(:), allocatable       :: fnames
  
   
    if (len_trim(md_1dfiles%structures) > 0) then
    
      fnamesstring = md_1dfiles%structures
      call strsplit(fnamesstring,1,fnames,1)
      call convertLongCulvertsAsNetwork(fnames(1), 0,md_culvertprefix,converted_fnamesstring,converted_crsdefsstring, istat)
      do ifil=2,size(fnames)
         call convertLongCulvertsAsNetwork(fnames(ifil), 1,md_culvertprefix, tempstring_fnames,tempstring_crsdef, istat)
         converted_crsdefsstring = trim(trim(converted_crsdefsstring)//', ')//tempstring_crsdef
         converted_fnamesstring  = trim(trim(converted_fnamesstring) //', ')//tempstring_fnames
      end do
      deallocate(fnames)
      call finalizeLongCulvertsInNetwork()
     
      call unc_write_net(trim(md_culvertprefix)//md_netfile, janetcell = 1, janetbnd = 0, jaidomain = 0, iconventions = UNC_CONV_UGRID)

      md_netfile = trim(md_culvertprefix)//md_netfile
      md_1dfiles%structures = converted_fnamesstring
      md_1dfiles%cross_section_definitions = converted_crsdefsstring   
      converted_fnamesstring = trim( trim(md_culvertprefix)//md_ident )//'.mdu'
      call writeMDUFile(converted_fnamesstring ,istat)
    endif
    
end subroutine makelongculverts_commandline

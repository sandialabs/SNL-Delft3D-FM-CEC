function File = qp_refresh(File)
%QP_REFRESH Refresh data for current data resource.

%----- LGPL --------------------------------------------------------------------
%                                                                               
%   Copyright (C) 2011-2022 Stichting Deltares.                                     
%                                                                               
%   This library is free software; you can redistribute it and/or                
%   modify it under the terms of the GNU Lesser General Public                   
%   License as published by the Free Software Foundation version 2.1.                         
%                                                                               
%   This library is distributed in the hope that it will be useful,              
%   but WITHOUT ANY WARRANTY; without even the implied warranty of               
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU            
%   Lesser General Public License for more details.                              
%                                                                               
%   You should have received a copy of the GNU Lesser General Public             
%   License along with this library; if not, see <http://www.gnu.org/licenses/>. 
%                                                                               
%   contact: delft3d.support@deltares.nl                                         
%   Stichting Deltares                                                           
%   P.O. Box 177                                                                 
%   2600 MH Delft, The Netherlands                                               
%                                                                               
%   All indications and logos of, and references to, "Delft3D" and "Deltares"    
%   are registered trademarks of Stichting Deltares, and remain the property of  
%   Stichting Deltares. All rights reserved.                                     
%                                                                               
%-------------------------------------------------------------------------------
%   http://www.deltaressystems.com
%   $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/tags/delft3dfm/141476/src/tools_lgpl/matlab/quickplot/progsrc/private/qp_refresh.m $
%   $Id: qp_refresh.m 140618 2022-01-12 13:12:04Z klapwijk $

if strcmp(File.FileType,'diff')
    File.Data(1) = qp_refresh(File.Data(1));
    File.Data(2) = qp_refresh(File.Data(2));
else
    qp_proxy('reload',File)
end

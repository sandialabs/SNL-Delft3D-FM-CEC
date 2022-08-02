function make_d3dmatlab(basedir,varargin)
%MAKE_D3DMATLAB Pre-compile Delft3D-MATLAB toolbox
%   Pre-compile MATLAB m-code to p-code for distribution as Delft3D-MATLAB toolbox.
%
%   MAKE_D3DMATLAB(BASEDIR)
%   Use specified directory instead of current directory as base directory

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
%   $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/tags/delft3dfm/141476/src/tools_lgpl/matlab/quickplot/make_d3dmatlab.m $
%   $Id: make_d3dmatlab.m 140743 2022-02-13 12:25:40Z jagers $

curdir = pwd;
addpath(curdir)
%
% comment lines for
% consistency with
% make_quickplot and
% make_ecoplot
%
if nargin>0
    cd(basedir);
end
err = [];
try
    localmake(varargin{:});
catch err
end
if nargin>0
    cd(curdir);
end
rmpath(curdir)
if ~isempty(err)
    rethrow(err)
end


function localmake(qpversion,T)
if ~exist('progsrc','dir')
    error('Cannot locate source folder "progsrc".')
end
sourcedir=[pwd,filesep,'progsrc'];

tdir = 'delft3d_matlab';
targetname = 'Delft3D-MATLAB interface';
targetdir = [pwd,filesep,tdir];

if ~exist(targetdir, 'dir')
    fprintf('Creating %s directory ...\n', tdir);
    mkdir(tdir);
end
cd(tdir)
% diary make_quickplot_diary % no diary to avoid clutter in the distribution folder ...

fprintf('Copying files ...\n');
exportsrc(sourcedir,targetdir)

fprintf('Including netCDF files ...\n');
if ~exist('netcdf','dir')
    mkdir('netcdf');
end
copyfile('../../../../third_party_open/netcdf/matlab/netcdfAll-4.1.jar','netcdf')
if ~exist('netcdf/mexnc','dir')
    mkdir('netcdf/mexnc');
    exportsrc('../../../../third_party_open/netcdf/matlab/mexnc', 'netcdf/mexnc')
end
if ~exist('netcdf/snctools','dir')
    mkdir('netcdf/snctools');
    exportsrc('../../../../third_party_open/netcdf/matlab/snctools', 'netcdf/snctools')
end

if nargin<2
    qpversion = read_identification(sourcedir, 'd3d_qp.m');
    T=now;
end
% strip off the 32/64 bit flag (the toolbox is platform independent)
qpversion = deblank(sscanf(qpversion,'%[^(]'));
% for the progress statement add the platform statement
qpversion_ = [qpversion, ' (all platforms)'];

TStr = datestr(T);
fprintf('\nBuilding %s version %s\n\n', targetname, qpversion_);
fprintf('Current date and time           : %s\n', TStr);

fprintf('Modifying files ...\n');
fstrrep([targetdir,filesep,'d3d_qp.m'], '<VERSION>', qpversion)
fstrrep([targetdir,filesep,'d3d_qp.m'], '<CREATIONDATE>', TStr)
fstrrep([targetdir,filesep,'Contents.m'], '<VERSION>', qpversion)
fstrrep([targetdir,filesep,'Contents.m'], '<CREATIONDATE>', TStr)

fprintf('Stripping files ...\n');
svnstripfile(targetdir)

%fprintf('Pcoding files ...\n');
%pmfile('dir',targetdir,targetdir,'-verbose')

fprintf('Cleaning up directory ...\n');
X = {'*.asv'
    '*.bak'
    '*.scc'
    'bin'
    'compileonly'};
cleanup(X)

fprintf('Removing unneeded subdirectories ...\n');
X = {'org'};
cleanup(X)
cd ..
fprintf('Finished.\n');


function exportsrc(sourcedir,targetdir)
d = dir(sourcedir);
for i = 1:length(d)
    source = [sourcedir filesep d(i).name];
    target = [targetdir filesep d(i).name];
    if d(i).isdir
        switch d(i).name
            case {'.','..','.svn'}
                % skip
            otherwise
                mkdir(target);
                exportsrc(source,target)
        end
    else
        copyfile(source,target)
    end
end
function test_nc_attput(mode)

if nargin < 1
	mode = nc_clobber_mode;
end

fprintf ('\t\tTesting NC_ATTPUT...  ' );

ncfile = 'foo.nc';
if exist(ncfile,'file')
    delete(ncfile);
end
nc_create_empty(ncfile,mode);

switch(mode)
	case nc_clobber_mode
		% netcdf-3
		run_common_tests(ncfile);

	case 'netcdf4-classic'
		run_common_tests(ncfile);
		verify_netcdf4(ncfile);

	case 'hdf4'
		run_common_tests(ncfile);
		run_hdf4_tests;
end
fprintf('OK\n');
return




%--------------------------------------------------------------------------
function run_hdf4_tests()
% HDF4 specific tests
test_hdf4_datastrs;
test_hdf4_cal;
test_hdf4_fillvalue;
return

%--------------------------------------------------------------------------
function run_common_tests(ncfile)
% write/retrieve a new double attribute
% write/retrieve a new float attribute
% write/retrieve a new int attribute
% write/retrieve a new short int attribute
% write/retrieve a new uint8 attribute
% write/retrieve a new int8 attribute
% write/retrieve a new text attribute
% write/read an empty attribute
% Verify that fill value attributes match the datatype of their dataset.
verify_fill_value_correctness(ncfile);

test_read_write_double_att ( ncfile );
test_read_write_float_att ( ncfile );
test_read_write_int_att ( ncfile );
test_read_write_short_att ( ncfile );
test_read_write_uint8_att ( ncfile );
test_read_write_int8_att ( ncfile );
test_read_write_char_att ( ncfile );
test_read_write_empty_att(ncfile);

return



%--------------------------------------------------------------------------
function verify_fill_value_correctness(ncfile)
% It was possible to write fill values with the wrong datatype on netcdf-3
% and not know what you did wrong until trying to read or write something.
% This test verifies that the _FillValue datatype is now correct.

nc_adddim(ncfile,'xx',5);
v.Name = 'yy';
v.Attribute.Name = '_FillValue';
v.Dimension = {'xx'};
v.Attribute.Value = single(-99);
nc_addvar(ncfile,v);

info = nc_getvarinfo(ncfile,'yy');
if ~strcmp(info.Attribute.Datatype,'double')
    error('failed');
end

%--------------------------------------------------------------------------
function test_read_write_empty_att(ncfile )
% BACKGROUND:  in R2008b, the TMW mex-file incorrectly disallowed empty 
% attributes, which are most definitely allowed. 
%
% REFERENCE:  http://www.mathworks.com/support/bugreports/609383


info = nc_info(ncfile);
if strcmp(info.Format,'HDF4')
    return
end

warning('off','SNCTOOLS:NCATTPUT:emptyAttributeBug');
nc_attput ( ncfile, nc_global, 'emptyAtt', '' );
x = nc_attget ( ncfile, nc_global, 'emptyAtt' );

if ( ~ischar(x) )
    error('class of retrieved attribute was not char.' );
end



v = version('-release');
mv = mexnc('inq_libvers');
switch(v)
case { '2008b','2009a','2009b'}
    if mv(1) == '4'
        % netcdf-4 capable mex-file
        if numel(x) ~= 0
            error('failed');
        end
    else
        if( numel(x) ~= 1 )
            error ( 'retrieved attribute was not one char in length' );
        end
    end
    %%% If you have applied the fix for bug #609383, then
    %%% comment out the code above and uncomment the
	%%% code below.
    %%% if ( numel(x) ~= 0 )
    %%%     error ( 'retrieved attribute was not empty' );
    %%% end
otherwise
    if ( numel(x) ~= 0 )
        error ( 'retrieved attribute was not empty' );
    end
end

warning('on','SNCTOOLS:NCATTPUT:emptyAttributeBug');
return



%--------------------------------------------------------------------------
function test_read_write_double_att ( ncfile )
% Verify that we can read/write double precision attributes.

nc_attput ( ncfile, nc_global, 'new_att', 0 );
x = nc_attget ( ncfile, nc_global, 'new_att' );

if ( ~strcmp(class(x), 'double' ) )
	error('class of retrieved attribute was not double.' );
end

if ( double(x) ~= 0 )
	error ( 'retrieved attribute was not same as written value' );
end

return




%--------------------------------------------------------------------------
function test_read_write_float_att ( ncfile )

nc_attput ( ncfile, nc_global, 'new_att2', single(0) );
x = nc_attget ( ncfile, nc_global, 'new_att2' );

if ( ~strcmp(class(x), 'single' ) )
	error('%class of retrieved attribute was not single.');
end
if ( double(x) ~= 0 )
	error ( 'retrieved attribute was not same as written value' );
end


%--------------------------------------------------------------------------
function test_read_write_int_att ( ncfile )

nc_attput ( ncfile, nc_global, 'new_att3', int32(0) );
x = nc_attget ( ncfile, nc_global, 'new_att3' );

if ( ~strcmp(class(x), 'int32' ) )
	error('class of retrieved attribute was not int32.');
end
if ( double(x) ~= 0 )
	error ( 'retrieved attribute was not same as written value' );
end


%--------------------------------------------------------------------------
function test_read_write_short_att ( ncfile )

nc_attput ( ncfile, nc_global, 'new_att4', int16(0) );
x = nc_attget ( ncfile, nc_global, 'new_att4' );

if ( ~strcmp(class(x), 'int16' ) )
	error('class of retrieved attribute was not int16.');
end
if ( double(x) ~= 0 )
	error ( 'retrieved attribute was not same as written value' );
end


%--------------------------------------------------------------------------
function test_read_write_uint8_att ( ncfile )

nc_attput ( ncfile, nc_global, 'new_att5', uint8(0) );
x = nc_attget ( ncfile, nc_global, 'new_att5' );

        
info = nc_info(ncfile);
if strcmp(info.Format,'HDF4')
    if ~strcmp(class(x), 'uint8' )
        error('class of retrieved attribute was not uint8.' );
    end
elseif strfind(info.Format,'NetCDF-4')
    if ~isa(x,'uint8')
        error('class of retrieved attribute was not uint8.' );
    end
elseif  ~strcmp(class(x), 'int8' )
    error('class of retrieved attribute was not int8.' );
end
if ( double(x) ~= 0 )
	error ( 'retrieved attribute was not same as written value' );
end


%--------------------------------------------------------------------------
function test_read_write_int8_att ( ncfile )

nc_attput ( ncfile, nc_global, 'new_att6', int8(0));
x = nc_attget ( ncfile, nc_global, 'new_att6' );

if  ~strcmp(class(x), 'int8' )
    error('class of retrieved attribute was not int8.' );
end

if ( double(x) ~= 0 )
	error ( 'retrieved attribute was not same as written value' );
end


%--------------------------------------------------------------------------
function test_read_write_char_att ( ncfile )

nc_attput ( ncfile, nc_global, 'new_att7', '0' );
x = nc_attget ( ncfile, nc_global, 'new_att7' );

if ( ~ischar(x ) )
	error('class of retrieved attribute was not char.');
end
if (x ~= '0' )
	error ( 'retrieved attribute was not same as written value' );
end


return





%--------------------------------------------------------------------------
function test_hdf4_datastrs (  )

nc_create_empty('foo.hdf','hdf4');
nc_adddim('foo.hdf','x', 4);
nc_attput('foo.hdf','x','coordsys','dud');
x = nc_attget('foo.hdf','x','coordsys');

if ( ~ischar(x ) )
	error('class of retrieved attribute was not char.');
end
if ~strcmp(x,'dud')
	error ( 'retrieved attribute was not same as written value' );
end


return

%--------------------------------------------------------------------------
function test_hdf4_cal (  )

nc_create_empty('foo.hdf','hdf4');
nc_adddim('foo.hdf','x', 4);
nc_attput('foo.hdf','x','scale_factor',1);
x = nc_attget('foo.hdf','x','scale_factor');
y = nc_attget('foo.hdf','x','add_offset');

if (x ~= 1) || (y ~= 0)
	error ( 'retrieved attribute was not same as written value' );
end


return

%--------------------------------------------------------------------------
function test_hdf4_fillvalue (  )

nc_create_empty('foo.hdf','hdf4');
nc_adddim('foo.hdf','x',4); 
nc_attput('foo.hdf','x','_FillValue',99);
x = nc_attget('foo.hdf','x','_FillValue');


if (x ~= 99)
	error ( 'retrieved attribute was not same as written value' );
end


return


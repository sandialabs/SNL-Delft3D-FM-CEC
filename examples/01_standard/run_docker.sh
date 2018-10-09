#!/bin/bash
    #
    # This script is an example for running Delft3D-FLOW in Docker (Linux container)
    # Adapt and use it for your own purpose
    #
    # 
    #
    # This script starts a single-domain Delft3D-FLOW computation on Linux
    #

    #
    # Set the config file here
    # 
argfile=config_d_hydro.xml




    #
    # Set the directories containing the binaries here
    #
export ARCH=lnx64
export D3D_HOME=/opt/delft3d_latest
flowexedir=$D3D_HOME/$ARCH/flow2d3d/bin
    #
    # No adaptions needed below
    #

    # Set some (environment) parameters
export LD_LIBRARY_PATH=$flowexedir:$LD_LIBRARY_PATH 

    # Run
$flowexedir/d_hydro.exe $argfile

    # Wait until all child processes are finished
wait

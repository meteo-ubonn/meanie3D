__author__ = "juergen.simon@uni-bonn.de"

import glob
import os
import os.path
import shutil
import time
from subprocess import call
from ..meanie3D import external
from ..meanie3D import utils

## Checks the configuration for mandatory items etc.
# \param The configuration to check.
# \throws RuntimeError in case something is amiss.
def checkConfiguration(conf):
    if not conf['data']:
        raise RuntimeError("Missing 'data' section")

    data = conf['data']
    if not data['variables']:
        raise RuntimeError("Missing 'data.variables'")

    # TODO: finish this
    return

## Runs a batch of files through the clustering and tracking.
#
# \param configuration with the following keys:
#
#    description : a description of the configuration
#    netcdf_dir : directory containing files to process
#    output_dir : directory to write results to
#    m3d_home   : ${MEANIE3D_HOME} environment variable
#    resume     : if True then pick off where you left off, if False
#                 all previous results are deleted and processing
#                 starts from scratch
#    meanie3D-detect : parameters for meanie3D-detect
#    meanie3D-track: parameters for meanie3D-track
#    scale : scale parameter (optional)
#    use_previous: True: using previous results to enhance tracking
#                  Absent or False: no use of previous results.
#    use_ci_score : Use the CI-Score algorithm
#
# \param time index the index in time dimension to read and process.
#
def run(config,time_index):

    print "---------------------------------------------------"
    print "Configuration: " + config['description']
    print "---------------------------------------------------"

    checkConfiguration(config)
    data = config['data']

    # Piece together the command line params for detection and tracking
    tracking_params = ""
    detect_params = ""

    if config['detection']:
        detect_params = "-v " + ','.join(data['variables'])
        detect_params = "%s -d %s" % (detect_params, ','.join(data['dimensions']))
        if data['vtk_dimensions']:
            detect_params = "%s --vtk-dimensions %s" % (detect_params, ','.join(data['vtk_dimensions']))
        if data['lowerThreshold']:
            detect_params = "%s --lower-thresholds %s" % (detect_params, ','.join(data['lowerThreshold']))
        if data['upperThreshold']:
            detect_params = "%s --upper-thresholds %s" % (detect_params, ','.join(data['upperThreshold']))
        if data['replacementValues']:
            detect_params = "%s --replacement-values %s" % (detect_params, ','.join(data['replacementValues']))

    if data['detection']:
        detection = data['detection']
        if detection['meanie3D-detect']:
            detect_params = "%s %s" % (detect_params, detection['meanie3D-detect'])

    if config['tracking']:
        tracking = config['tracking']
        if tracking['meanie3D-track']:
            tracking_params = "%s %s" % (tracking_params, tracking['meanie3D-track'])
        if tracking['histogramVariable']:
            tracking_params = "%s -t %s" % (tracking_params, tracking['histogramVariable'])
        if tracking['histogramWeight']:
            tracking_params = "%s --wt %s" % (tracking_params, tracking['histogramWeight'])
        if tracking['positionWeight']:
            tracking_params = "%s --wr %s" % (tracking_params, tracking['positionWeight'])
        if tracking['sizeWeight']:
            tracking_params = "%s --ws %s" % (tracking_params, tracking['sizeWeight'])
        if tracking['maxSpeed']:
            tracking_params = "%s --max-speed %f" % (tracking_params, tracking['maxSpeed'])
        if tracking['maxTime']:
            tracking_params = "%s --max-time %d" % (tracking_params, tracking['maxTime'])
        if tracking['useDisplacementVectors'] == True:
            tracking_params = "%s -v"

    resume_at_index = 0;

    # In case a scale T parameter is given, the output dir
    # is scale<T>. Otherwise it's 'clustering'
    output_dir = config['output_dir']
    if (config.get('scale') == "None"):
        output_dir = output_dir + "/clustering"
    else:
        output_dir = output_dir + "/scale"+str(config['scale'])

    print "Writing output to " + output_dir

    # Resume?
    if config['resume'] == False:
        # consider time index. Even if not resuming, the
        # output directories should only be created at
        # the first time step
        if time_index <= 0:
            print "Removing results from previous runs"
            utils.create_ouput_directories(output_dir)
    else:
        resume_at_index = utils.number_of_netcdf_files(output_dir+"/netcdf")
        print "Resuming at index " + str(resume_at_index)

    # Get a list of the files we need to process

    netcdf_pattern = config['source_directory'] + "/*.nc"
    netcdf_list=sorted(glob.glob(netcdf_pattern))
    last_cluster_file=""
    run_count = 0

    # Process the files one by one

    for netcdf_file in netcdf_list:

        basename = os.path.basename(netcdf_file)
        cluster_file= ""

        if time_index < 0:
            cluster_file= output_dir+"/netcdf/" + os.path.splitext(basename)[0] + "-clusters.nc"
        else:
            cluster_file= output_dir+"/netcdf/" + os.path.splitext(basename)[0] +"-clusters_" +str(time_index) + ".nc"
            last_cluster_file = output_dir+"/netcdf/" + os.path.splitext(basename)[0] +"-clusters_" +str(time_index-1) + ".nc"

        # if there is a resume counter, keep skipping
        # until the count is right
        if (resume_at_index > 0) and (run_count < resume_at_index):
            last_cluster_file = cluster_file
            run_count = run_count + 1
            continue

        print "-------------------------------------------------------------------------------------"
        print "Processing " + netcdf_file
        print "-------------------------------------------------------------------------------------"

        # ----------------------------------------------
        # Clustering
        # ----------------------------------------------

        if config['detection']:
            print "-- Clustering --"

            # build the clustering command
            params = detect_params

            # amend for time index and make proper logfile name
            if time_index < 0:
                logfile = output_dir+"/log/clustering_" + str(run_count)+".log"
            else:
                params = params + " -t " + str(time_index)
                logfile = output_dir+"/log/clustering_" + str(time_index)+".log"

            # scale?
            if not config.get('scale') == "None":
                scale_param = " -s " + config['scale']
                params += scale_param

            # use previous result to enhance current?
            if ((run_count > 0) or (time_index > 0)) and detection['use_previous']:
                params += " -p " + last_cluster_file

            # add ci-comparison-file if applicable
            if run_count >= 3 and config['use_ci_score']:
                params += " --ci-comparison-file " + netcdf_list[run_count-3]
                proto_file = os.path.splitext(os.path.basename(netcdf_list[run_count-3]))[0] + "-protoclusters.nc"
                params += " --ci-comparison-protocluster-file " + proto_file;

            # complete command with directing output to logfile
            params = params + " > " + logfile

            # execute
            start_time = time.time()
            print "meanie3D-detect " + params
            print external.execute_command("meanie3D-detect",params,True)
            print "Finished. (%.2f seconds)" % (time.time()-start_time)

        # ----------------------------------------------
        # Tracking
        # ----------------------------------------------

        if config['tracking']:
            # if we have a previous scan, run the tracking command
            if (run_count > 0) or (time_index > 0):

                if time_index < 0:
                    logfile = output_dir+"/log/tracking_" + str(run_count)+".log"
                else:
                    logfile = output_dir+"/log/tracking_" + str(time_index)+".log"

                print "-- Tracking --"
                params = tracking_params + " -p "+last_cluster_file+" -c "+cluster_file
                params = params + " > " + logfile

                # execute
                start_time = time.time()
                print "meanie3D-track" + params
                print external.execute_command("meanie3D-track",params,True)
                print "Finished. (%.2f seconds)" % (time.time()-start_time)

        # keep track
        last_cluster_file = cluster_file
        # don't forget to increment run counter
        run_count = (run_count + 1)
        # Clean out the trash if there is any
        print "Cleaning up *.vt*"
        call("rm -f *.vt*", shell=True)

    print "Done."
    return
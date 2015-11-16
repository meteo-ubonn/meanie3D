import glob
import os
import os.path
import time
from subprocess import call

import external
import utils

## Checks the configuration for mandatory items etc.
# \param The configuration to check.
# \throws RuntimeError in case something is amiss.
def checkConfiguration(conf):
    if not conf['data']:
        raise RuntimeError("Missing 'data' section")

    data = conf['data']
    if not 'variables' in data:
        raise RuntimeError("Missing 'data.variables'")
    if not 'dimensions' in data:
        raise RuntimeError("Missing 'data.dimensions'")


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
    data = utils.getSafe(config,'data')
    detection = utils.getSafe(config,'detection')
    tracking = utils.getSafe(config,'tracking')

    # Piece together the command line params for detection and tracking
    tracking_params = ""
    detect_params = ""

    if data and detection:

        variables = list()
        lowerThresholds = list()
        upperThresholds = list()
        replacementValues = list()

        detect_params = "%s -d %s" % (detect_params, ','.join(data['dimensions']))

        for variable in data['variables']:
            variables.append(variable['name'])
            if (utils.getSafe(variable,'lowerThreshold')):
                lowerThresholds.append("%s=%s" % (variable['name'], str(variable['lowerThreshold'])))

            if (utils.getSafe(variable,'upperThreshold')):
                upperThresholds.append("%s=%s" % (variable['name'], str(variable['upperThreshold'])))

            if (utils.getSafe(variable,'replacementValue')):
                replacementValues.append("%s=%s" % (variable['name'], str(variable['replacementValue'])))

        detect_params += " -v %s" % ','.join(variables)

        if lowerThresholds:
            detect_params += " --lower-thresholds %s" % ','.join(lowerThresholds)

        if upperThresholds:
            detect_params += " --upper-thresholds %s" % ','.join(upperThresholds)

        if replacementValues:
            detect_params += " --replacement-values %s" % ','.join(lowerThresholds)

        if utils.getSafe(data,'vtkDimensions'):
            detect_params += " --vtk-dimensions %s" % ','.join(data['vtkDimensions'])

    if detection:
        minClusterSize = utils.getSafe(detection,'minClusterSize')
        if minClusterSize:
            detect_params = "%s -m %d" % (detect_params, minClusterSize)

        additionalParams = utils.getSafe(detection,'meanie3D-detect')
        if detect_params:
            detect_params = "%s %s" % (detect_params, additionalParams)

    if tracking:
        
        if utils.getSafe(tracking,'meanie3D-track'):
            tracking_params = "%s %s" % (tracking_params, tracking['meanie3D-track'])
            
        if utils.getSafe(tracking,'histogramVariable'):
            tracking_params = "%s -t %s" % (tracking_params, tracking['histogramVariable'])
            
        if utils.getSafe(tracking,'histogramWeight'):
            tracking_params = "%s --wt %s" % (tracking_params, tracking['histogramWeight'])
        else:
            tracking_params = "%s --wt 0" % tracking_params
            
        if utils.getSafe(tracking,'positionWeight'):
            tracking_params = "%s --wr %s" % (tracking_params, tracking['positionWeight'])
        else:
            tracking_params = "%s --wr 0" % tracking_params
            
        if utils.getSafe(tracking,'sizeWeight'):
            tracking_params = "%s --ws %s" % (tracking_params, tracking['sizeWeight'])
        else:
            tracking_params = "%s --ws 0" % tracking_params
            
        if utils.getSafe(tracking,'maxSpeed'):
            tracking_params = "%s --max-speed %f" % (tracking_params, tracking['maxSpeed'])
            
        if utils.getSafe(tracking,'maxTime'):
            tracking_params = "%s --max-time %d" % (tracking_params, tracking['maxTime'])

        if utils.getSafe(tracking,'maxSizeDeviation'):
            tracking_params = "%s --max-size-deviation %d" % (tracking_params, tracking['maxSizeDeviation'])

        if utils.getSafe(tracking,'useDisplacementVectors'):
            tracking_params = "%s -v"  % tracking_params

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

        print "-------------------------------------------------------------------------------------"
        print "Processing " + netcdf_file
        print "-------------------------------------------------------------------------------------"

        # ----------------------------------------------
        # Clustering
        # ----------------------------------------------

        if detection:
            # if there is a resume counter, keep skipping until the count is right
            # Note: this only applies to clustering, which is expensive. Tracking
            # will be re-run
            if (resume_at_index > 0) and (run_count < resume_at_index):
                last_cluster_file = cluster_file
                run_count = run_count + 1
                continue

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
            if ((run_count > 0) or (time_index > 0)) and detection['usePrevious']:
                params += " -p " + last_cluster_file

            # add ci-comparison-file if applicable
            if run_count >= 3 and utils.getSafe(detection,'useCIScore'):
                params += " --ci-comparison-file " + netcdf_list[run_count-3]
                proto_file = os.path.splitext(os.path.basename(netcdf_list[run_count-3]))[0] + "-protoclusters.nc"
                params += " --ci-comparison-protocluster-file " + proto_file;

            # Input file
            params += " -f %s -o %s" % (netcdf_file,cluster_file)

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

        if tracking:
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
        if config['cleanup_vtk']:
            print "Cleaning up *.vt*"
            call("rm -f *.vt*", shell=True)

    print "Done."
    return
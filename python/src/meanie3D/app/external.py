__author__ = 'simon'

''' This module contains code used to call external programs '''

import os
import subprocess
import platform
from subprocess import Popen

# ---------------------------------------------
# Locating commands
# ---------------------------------------------

# Constant denoting failure to locate a command
COMMAND_NOT_FOUND="NOT_FOUND"
COMMAND_MAP={}

## Get DYLD_LIBRARY_PATH depending on operating system. OSX needs
# -------------------------------------------------------------------
# special work because of the homebrew gfx libraries, which get
# in the way of the system libraries.
# \return DYLD_LIBRARY_PATH
def get_dyld_library_path():
    path = "/usr/local/lib"
    if platform.system() == 'Darwin':
        path = "/System/Library/Frameworks/ImageIO.framework/Versions/A/Resources/:"+path
    return path


# -------------------------------------------------------------------
## Get the complete shell command with added DYLD_LIBRARY_PATH etc.
# for the given executable.
# \param executable name (like meanie3D-detect, meanie3D-track etc.)
# \return shell command to run the binary
def get_executable_command(executable):
    bin_prefix = "export DYLD_LIBRARY_PATH="+get_dyld_library_path()+";"
    command = bin_prefix + executable
    return command


# Returns the standard search path for the binaries
# @return search path
def command_search_paths():
    paths=['/usr/local/bin','/usr/bin/','/bin','/sbin','/Applications/VisIt.app/Contents/Resources/bin']
    return paths

# Checks if the given command can be found in the given path.
# If the recurse flag is set, all subdirectories are searched recursively.
#
# @param command
# @param path
# @param recurse
# @return fully qualified command path or "NOT_FOUND"
# @throws IOError if a command is not executable
def find_command_in_path(command, path, recurse):
    result = COMMAND_NOT_FOUND
    if os.path.exists(path) and os.path.isdir(path):
        files = os.listdir(path)
        for file_ in files:
            if file_ == command:
                # Check if the entry is a file and if it is executable
                full_path = path + os.sep + command
                if os.path.isfile(full_path):
                    if os.access(full_path, os.X_OK):
                        return full_path
                    else:
                        raise IOError('Command ' + full_path + ' is not executable.')

        # If it gets here, the command was not located.
        if recurse:
            for file_ in files:
                if os.path.isdir(file_):
                    full_path = path + os.pathsep + file_
                    result = find_command_in_path(command, full_path, recurse)
                    if result != COMMAND_NOT_FOUND:
                        return result

    return result


## Attempts to locate the given executables in the given filesystem paths.
#
# \param command_list - list of strings containing the commands
# \param path_list - list of strings containing the search paths
# \param recurse - boolean. If <true> the method recurses into each path.
# \return a dictionary containing the command names mapping to the command paths
# \throws IOError if a command can not be located or is not executable
def find_ext_cmds_in_paths(command_list, path_list, recurse):
    for command in command_list:
        result = COMMAND_NOT_FOUND
        for path in path_list:
            result = find_command_in_path(command, path, recurse)
            if result != COMMAND_NOT_FOUND:
                COMMAND_MAP[command] = result
                break
        if result == COMMAND_NOT_FOUND:
            raise IOError('Could not locate command '+command)
    return COMMAND_MAP


# Attempts to locate the given executables in the standard filesystem
# paths (/usr/local/bin /var/opt/bin /usr/bin)
#
# @param list of strings containing the commands
# @return a dictionary containing the command names mapping to the command paths
# @throws IOError if a command can not be located or is not executable
def find_ext_cmds(command_list, recurse=False):
    return find_ext_cmds_in_paths(command_list, command_search_paths(), recurse)


## Executes the given command and returns the code
#
# \param command
# \param parameters
# \return process returncode
# \throws IOError if a command is not executable
# \throws Exception if a command was not searched for before attempting to run
def execute_command(command, parameters, with_stdout=False):
    # Make sure we searched for the command first
    if not COMMAND_MAP.has_key(command):
        raise Exception(
            "Requested command: %s was not found in the map. All commands must be searched before they can be used")

    cmd = get_executable_command(COMMAND_MAP.get(command)) + ' ' + parameters
    # print("Going to execute, with_stdout: %s, cmd: %s" % (with_stdout, cmd))

    p = Popen(cmd, shell=True, stdout=subprocess.PIPE)
    if with_stdout:
        out, err = p.communicate()
        return out
    else:
        p.communicate() # We don't care about the output
        # print("returncode: %s" % p.returncode)
        return p.returncode


## Executes the given command and returns the command's output
# \param command
# \param parameters
# \return process returncode
# \throws IOError if a command is not executable
# \throws Exception if a command was not searched for before attempting to run
def execute_command_with_stdout(command, parameters):
    return execute_command(command, parameters, True)


## Executes the given command and returns the command's output. Throws
# exception if the return code of the command process is non zero.
#
# \param command The command to execute.
# \param parameters Parameters to pass to the command.
# \return Stdout of the command.
# \throws IOError: If a command is not executable.
# \throws subprocess.CalledProcessError: If return code is non zero.
#
def execute_command_with_check(command, parameters):
    if not COMMAND_MAP.has_key(command):
        # User did not locate this command. Try to recover using standard
        # search paths
        print("WARNING: command {} was not verified. Attempting to locate it now.".format(command))
        find_ext_cmds([command], True)

    cmd = "{} {}".format(get_executable_command(COMMAND_MAP.get(command)), parameters)
    # print(cmd)

    p = Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

    if err:
        # Using warning logging level here as subprocess output to STDERR is
        # not necessarily indicate an overall error.
        print("ERROR: STDERR in subprocess {}: {}".format(command, err.strip()))

    if p.returncode:
        raise subprocess.CalledProcessError(p.returncode, cmd, out)

    return out


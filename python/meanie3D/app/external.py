'''
The MIT License (MIT)

(c) Juergen Simon 2014 (juergen.simon@uni-bonn.de)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

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
def __get_dyld_library_path():
    path = "/usr/local/lib"
    if platform.system() == 'Darwin':
        path = "/System/Library/Frameworks/ImageIO.framework/Versions/A/Resources/:"+path
    return path

# Returns the standard search path for the binaries
# @return search path
def __command_search_paths():
    paths=['/usr/local/bin',
           '/usr/bin/',
           '/bin',
           '/sbin',
           '/usr/local/visit/bin',
           '/Applications/VisIt.app/Contents/Resources/bin']
    return paths

# -------------------------------------------------------------------
## Get the complete shell command with added DYLD_LIBRARY_PATH etc.
# for the given executable.
# \param executable name (like meanie3D-detect, meanie3D-track etc.)
# \return shell command to run the binary
def getCommand(executable):
    bin_prefix = "export DYLD_LIBRARY_PATH="+__get_dyld_library_path()+";"
    command = bin_prefix + executable
    return command

# Checks if the given command can be found in the given path.
# If the recurse flag is set, all subdirectories are searched recursively.
#
# @param command
# @param path
# @param recurse
# @return fully qualified command path or "NOT_FOUND"
# @throws IOError if a command is not executable
def locateCommandInPath(command, path, recurse):
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
                    result = locateCommandInPath(command, full_path, recurse)
                    if result != COMMAND_NOT_FOUND:
                        return result

    return result


def hasCommand(command):
    """
    Tests if the command is available
    @param command name
    @returns True or False
    """
    try:
        locateCommands([command])
    except IOError:
        return False
        

def locateCommandsInPaths(command_list, path_list, recurse):
    """
    Attempts to locate the given executables in the given filesystem paths.
    @param command_list - list of strings containing the commands
    @param path_list - list of strings containing the search paths
    @param recurse - boolean. If <true> the method recurses into each path.
    @returns a dictionary containing the command names mapping to the command paths
    @throws IOError if a command can not be located or is not executable
    """
    for command in command_list:
        result = COMMAND_NOT_FOUND
        for path in path_list:
            result = locateCommandInPath(command, path, recurse)
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
def locateCommands(command_list, recurse=True):
    return locateCommandsInPaths(command_list, __command_search_paths(), recurse)

## Executes the given command and returns the code
#
# \param command
# \param parameters
# \return process returncode
# \throws IOError if a command is not executable
# \throws Exception if a command was not searched for before attempting to run
def execute_command(command, parameters, withStdOut=False, silent=True):
    # Make sure we searched for the command first
    if not COMMAND_MAP.has_key(command):
        raise Exception(
            "Requested command: %s was not found. All commands must be located before they can be used")

    cmd = getCommand(COMMAND_MAP.get(command)) + ' ' + parameters
    if withStdOut:
        p = Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        return p.returncode, out
    else:
        if silent:
            p = Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            # TODO: this does seem to be broken
            p = Popen(cmd, shell=True, stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
        p.communicate()
        return p.returncode

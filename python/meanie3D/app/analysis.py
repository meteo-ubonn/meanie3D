'''
Analyse and visualise tracking runs
'''

import getopt
import os
import sys
import webbrowser
import json

import meanie3D

def loadTrackingDictionary(path):
    '''
    Loads a tracking dictionary as output by meanie3D-trackstats -t
    :param path:
    :return: parsed JSON dictionary
    '''
    return

def constructTreeJson(trackingDictionary):
    '''
    Parses a tracking dictionary and creates data for visualising
    the tracks in d3.js
    :param trackingDictionary:
    :return:
    '''

def showTrackTree(tree):
    '''
    Displays the given tree JSON in d3.js
    :param tree:
    :return:
    '''
    return


def usage():
    '''
    Prints help and exits
    :return:
    '''
    print "meanie3D-trackgraph -f <tracking dictionary file>"
    print "Analyses a track dictionary and shows a track graph with splits, merges etc."
    print "-f : tracking dictionary file"
    print "--help, -h  : print this message and exits."
    print "--version   : prints the version information and exits"
    sys.exit(1)
    return


def run():
    '''
    Parses the command line and performs the analysis of the track dictionary.
    :return:
    '''

    # Parse command line
    try:
        argv = sys.argv[1:]
        opts, args = getopt.getopt(argv, "c:f:s:o:r:h", ["json-example","resume","help","version","start=","end="])
    except getopt.GetoptError as detail:
        print detail
        sys.exit(2)

    num_params = 0
    dictionaryPath = None

    for o, a in opts:

        if o in ['--file','-f']:
            dictionaryPath = a
            num_params = num_params + 1

        elif o in ["--help"]:
            usage()
            sys.exit()

        elif o in ["--version"]:
            meanie3D.getVersion()

        else:
            usage()

    if num_params < 2:
        usage()

    uses_time = False
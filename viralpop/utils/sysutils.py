#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import gzip
import shutil
import tempfile
import argparse
from subprocess import check_output, check_call, Popen, STDOUT, PIPE
from subprocess import CalledProcessError

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

class PipelineStepError(Exception):
    pass

"""
def run_command(cmd, out_fn):
    with open(out_fn, 'w') as outh:
        p = Popen(cmd, stderr=PIPE, stdout=outh)
        ret = p.communicate()
    if ret[1]: print >>sys.stderr, ret[1]
    if p.returncode != 0:
        sys.exit('Command %s failed with exit code %d' % (cmd[0], p.returncode))
    return
"""

def check_dependency(prog):
    """ Check whether shell command can be called
    """
    try:
        _ = check_output('which %s' % prog, stderr=STDOUT, shell=True)
    except CalledProcessError:
        raise PipelineStepError('Dependency "%s" not found.' % prog)

def command_runner(cmds, stage=None, debug=False):
    """ Run a list of commands
    """
    # Formatted print of each command
    for i,args in enumerate(cmds):
        print >>sys.stderr, '\n[--- %s command %d ---]' % (stage, (i+1))
        s = '%s' % args[0]
        prev = 'init'
        for a in args[1:]:
            if a.startswith('-'):
                s += '\n    %s' % a
                prev = 'opt'
            elif a in ['>', '>>', '2>', '&>', '|', ]:
                s += '\n        %s' % a
                prev = 'redir'
            else:
                if prev == 'opt':
                    s += ' %s' % a
                elif prev == 'redir':
                    s += ' %s' % a
                else:
                    s += '\n    %s' % a
                prev = 'val'
        print >>sys.stderr, s
    # Join each command with whitespace and join commands with "&&"
    cmdstr = ' && '.join(' '.join(c) for c in cmds)
    
    if debug:
        # Print the joined command
        print >>sys.stderr, '\n[--- %s commands ---]' % stage
        print >>sys.stderr, cmdstr
        print >>sys.stderr, '\n[--- %s ---]' % stage        
    else:
        # Call using subprocess
        p = check_call(cmdstr, shell=True)

"""
Helpers for parsing command-line arguments
"""
def existing_file(f):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.isfile(f):
        raise argparse.ArgumentTypeError("{0} does not exist".format(f))
    return f

def existing_dir(f):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.isdir(f):
        raise argparse.ArgumentTypeError("{0} does not exist".format(f))
    return f

def args_params(args):
    """ Returns a dictionary from argparse namespace
        Excludes "func" argument
    """
    d = {k:v for k,v in vars(args).iteritems() if v is not None}
    if 'func' in d: d.pop('func')
    return d

"""
def command_runner_stdout(cmd, out_fn):
    with open(out_fn, 'w') as outh:
        p = Popen(cmd, stderr=PIPE, stdout=outh)
        ret = p.communicate()
    if ret[1]: print >>sys.stderr, ret[1]
    if p.returncode != 0:
        sys.exit('Command %s failed with exit code %d' % (cmd[0], p.returncode))
    return
"""

def create_tempdir(step='pipeline_step', basedir=None, quiet=False):
    """ Creates temporary directory
    """
    import tempfile
    checkdirs = ['/tmp', '/scratch', '/Temp']
    # Temporary directory
    if basedir is None:
        if 'TMPDIR' in os.environ:
            basedir = os.environ['TMPDIR']
        else:
            for d in checkdirs:
                if os.path.isdir(d):
                    basedir = d
                    break
    
    if not basedir or not os.path.isdir(basedir):
        raise PipelineStepError("Could not identify temporary directory")
    
    curdir = tempfile.mkdtemp(prefix='tmp_%s' % step, dir=basedir)
    if not quiet:
        print >>sys.stderr, '\n[--- %s ---] Using temporary directory %s' % (step, curdir)
    return curdir

def remove_tempdir(d, step='pipeline_step', quiet=False):
    """ Removes temporary directory
    """
    if os.path.isdir(d):
        if not quiet:
            print >>sys.stderr, '\n[--- %s ---] Removing temporary directory %s' % (step, d)
        shutil.rmtree(d)


def get_filehandle(fh):
    """ Resolve string or filehandle to filehandle
    """
    if isinstance(fh, basestring):
        if os.path.splitext(fh)[-1] == '.gz':
            return gzip.open(fh, 'rb')
        else:
            return open(fh, 'rU')
    else:
        return fh

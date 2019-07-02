# -*- coding: utf-8 -*-

from __future__ import print_function
from builtins import str
from past.builtins import basestring
import sys
import os
import gzip
import shutil
import tempfile
import argparse
import subprocess


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"


class PipelineStepError(Exception):
    pass


class MissingRequiredArgument(Exception):
    pass


class ArgumentDefaultsHelpFormatterSkipNone(argparse.HelpFormatter):
    """Help message formatter which adds default values to argument help.

    Only the name of this class is considered a public API. All the methods
    provided by the class are considered an implementation detail.

    Modified from argparse.HelpFormatter to skip default when it is "None"

    """
    def _get_help_string(self, action):
        help = action.help
        if '%(default)' not in action.help:
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    if action.default is not None:
                        help += ' (default: %(default)s)'
        return help


def check_dependency(prog):
    """ Check whether shell command can be called
    """
    try:
        _ = subprocess.check_output('which %s' % prog, stderr=subprocess.STDOUT, shell=True)
    except subprocess.CalledProcessError:
        raise PipelineStepError('Dependency "%s" not found.' % prog)


def determine_dependency_path(choices):
    messages = ['Dependencies not found: ']
    for prog in choices:
        try:
            check_dependency(prog)
            return prog
        except PipelineStepError as e:
            messages.append(str(e))
    raise PipelineStepError('\n'.join(messages))


def log_message(msg, quiet, logfile):
    if not quiet:
        try:
            sys.stderr.write(msg.encode('utf-8')) # python2
        except TypeError:
            sys.stderr.write(msg)  # python3
    if logfile is not None:
        try:
            logfile.write(msg.encode('utf-8')) # python2
        except TypeError:
            logfile.write(msg) # python3


def pretty_print_commands(cmds, stage, out_fh=sys.stderr):
    # Formatted print of each command
    for i,args in enumerate(cmds):
            print('\n[--- %s command %d ---]' % (stage, (i+1)), file=out_fh)
            s = '%s' % args[0]
            prev = 'init'
            for a in args[1:]:
                if a.startswith('-'):
                    s += ' \\\n    %s' % a
                    prev = 'opt'
                elif a in ['>', '>>', '2>', '&>', '|', ]:
                    s += ' \\\n        %s' % a
                    prev = 'redir'
                else:
                    if prev == 'opt':
                        s += ' %s' % a
                    elif prev == 'redir':
                        s += ' %s' % a
                    else:
                        s += ' \\\n    %s' % a
                    prev = 'val'
            print(s, file=out_fh)


def command_runner(cmds, stage=None, quiet=False, logfile=None, debug=False):
    """ Run a list of commands
    """
    # Join each command with whitespace and join commands with "&&"
    cmdstr = ' && '.join(' '.join(c) for c in cmds)

    if debug:
        # Print the joined command
        print('\n[--- %s commands ---]' % stage, file=sys.stderr)
        print(cmdstr, file=sys.stderr)
        print('\n[--- %s ---]' % stage, file=sys.stderr)
        return

    if not quiet:
        pretty_print_commands(cmds, stage, sys.stderr)

    if logfile is not None:
        pretty_print_commands(cmds, stage, logfile)

    p = subprocess.Popen(
        cmdstr, shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    for line in p.stdout:
        log_message(line.decode("utf-8"), quiet, logfile)
    p.wait()
    return


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
    d = {k:v for k,v in list(vars(args).items()) if v is not None}
    if 'func' in d: d.pop('func')
    return d


def create_tempdir(step='HPstep', basedir=None, quiet=False, logfile=None):
    """ Creates temporary directory
    """
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
    
    curdir = tempfile.mkdtemp(prefix='tmpHP_%s' % step, dir=basedir)
    msg = '\n[--- %s ---] Using temporary directory %s\n' % (step, curdir)
    log_message(msg, quiet, logfile)
    return curdir


def remove_tempdir(d, step='HPstep', quiet=False, logfile=None):
    """ Removes temporary directory
    """
    if os.path.isdir(d):
        msg = '\n[--- %s ---] Removing temporary directory %s\n' % (step, d)
        log_message(msg, quiet, logfile)
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


def get_java_heap_size():
    """ Determine a reasonable JVM heap size for this system

     Have not yet implemented the logic here. Just returns 32

    Returns:
        heap_size (int): Heap size in GB

    """
    return 32

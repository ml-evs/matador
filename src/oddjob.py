#!/usr/bin/python
# coding :utf-8
""" Submits the same job across multiple nodes
on a cluster with a shared filesystem, and no queuing,
handling killed jobs and crashes.
"""

import argparse
import multiprocessing as mp
import subprocess as sp
from os import getcwd
from utils.print_utils import print_notify, print_warning
from traceback import print_exc
from sys import stdout


class Hive:
    """ Implements the spawning of workers. """
    def __init__(self, *args, **kwargs):
        self.args = kwargs
        self.command = self.args.get('command')
        self.nodes = self.args.get('nodes')
        if len(set(self.nodes)) < len(list(self.nodes)):
            print_warning('Skipping duplicate nodes...')
        self.nodes = set(map(int, self.nodes))

    def spawn(self):
        """ Spawn processes to perform calculations (borrowed from run3). """
        procs = []
        for node in self.nodes:
            procs.append(mp.Process(target=self.run_command, args=[node]))
        try:
            for proc in procs:
                proc.start()
        except(KeyboardInterrupt, SystemExit, RuntimeError):
            print_exc()
            for proc in procs:
                proc.terminate()
            exit('Killing running jobs and exiting...')
        except:
            print_exc()
            pass

    def run_command(self, node):
        """ Parse command and run on particular node. """
        cwd = getcwd()
        compute_command = self.command.replace('$ALL_CORES', str(mp.cpu_count()))
        print_notify('Executing {} on node{}...'.format(compute_command, node))
        process = sp.Popen(['ssh', 'node{}'.format(node),
                            'cd', '{};'.format(cwd),
                            compute_command], stdout=stdout, shell=False)
        return process


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='oddjob',
        description='Run a single script on multiple nodes, numbered node<n>.',
        epilog='Written by Matthew Evans (2017).')
    parser.add_argument('command', type=str,
                        help='command to run, with arguments - must be apostrophized. \
                        Use the $ALL_CORES macro to use, unsurprisingly, all the cores \
                        on the node in question. e.g. oddjob \'pyairss -c $ALL_CORES -ha \
                        <seed>\' -n 1 16')
    parser.add_argument('-n', '--nodes', type=str, nargs='+',
                        help='list node numbers to run job on with space delimiters, e.g. \
                        -n 3 14 15')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='debug output')
    args = parser.parse_args()
    runner = Hive(command=args.command, nodes=args.nodes, debug=args.debug)
    try:
        runner.spawn()
    except(KeyboardInterrupt, SystemExit):
        exit('Exiting top-level...')

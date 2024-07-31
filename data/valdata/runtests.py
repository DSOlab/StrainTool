#! /usr/bin/python
#-*- coding: utf-8 -*-

import sys
import os
import argparse
import subprocess
import atexit

def fileql(f1, f2, start_at_line=None):
    cmp_started = 0
    with open(f1, 'r') as fin1:
        with open(f2, 'r') as fin2:
            for l1,l2 in zip(fin1.readlines(), fin2.readlines()):
                if start_at_line:
                    if cmp_started:
                        if l1.strip() != l2.strip():
                            return False
                    else:
                        if l1.strip() == start_at_line:
                            cmp_started = 1
                else:
                    if l1.strip() != l2.strip():
                        return False
    if start_at_line and (cmp_started == 0):
        return False
    return True

def cleanup(fstr):
    if os.path.isfile(fstr): os.remove(fstr)

class myFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawTextHelpFormatter):
    pass

parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description='Automatic testing of StrainTool.',
    epilog=('''National Technical University of Athens,
    Dionysos Satellite Observatory\n
Send bug reports to:
  Xanthos Papanikolaou, xanthos@mail.ntua.gr
  Dimitris Anastasiou,danastasiou@mail.ntua.gr
July, 2024'''))

parser.add_argument('-p', '--root-dir',
    default=argparse.SUPPRESS,
    metavar='ROOT_DIR',
    dest='root_dir',
    required=True,
    help='Path to the project\'s root directory (i.e. where the pystrain and bin folders are located)')

data_files = ['midas001.vel', 'midas002.vel', 'midas003.vel', 'midas004.vel']
test_step = [0.5, 1]
test_wt = [6, 12, 24, 48]
test_dmax = [300, 500, 700]

if __name__ == '__main__':
    args  = parser.parse_args()

    if not os.path.isdir(args.root_dir):
        print('ERROR. Failed locating directory: {:}'.format(args.root_dir))
        sys.exit(1)

    data_dir = os.path.join(args.root_dir, 'data/valdata')

    exe = os.path.join(args.root_dir, 'bin/StrainTensor.py')
    if not os.path.isfile(exe):
        print('ERROR. Failed locating main executable : {:}'.format(exe))
        sys.exit(1)

    pid = os.getpid()
    count = 0

    for tf in data_files:
        fn = os.path.join(data_dir, tf)
        if not os.path.isfile(fn):
            print('ERROR. Failed locating test data file : {:}'.format(fn))
            sys.exit(1)

        for step in test_step:
            for wt in test_wt:
                for dmax in test_dmax:

                    reffn = '0' + tf.split('.')[0][-1] + str(int(step)) + '00' + '{0:02d}'.format(wt) + '01' + str(dmax) + '0210'

                    cmd = ['-i', fn, '--x-grid-step', '{:.2f}'.format(step), '--y-grid-step', '{:.2f}'.format(step), '--Wt', '{:d}'.format(wt), '--dmax', '{:d}'.format(dmax), '-g']
                    # print([exe] + cmd)
                    with open('.{:}'.format(pid), 'w') as fout:
                        result = subprocess.run([exe] + cmd, stdout=fout, stderr=subprocess.DEVNULL, check=False)

                    atexit.register(cleanup, '.{:}'.format(pid))

                    if result.returncode != 0:
                        print('ERROR Failed executing command: {:}'.format(' '.join([exe]+cmd)), file=sys.stderr)
                        sys.exit(5)

                    for ref_type, start_at in zip(['_station_info.dat', '_strain_info.dat', '_strain_stats.dat'], [None,None,"Longtitude  Latitude  # stations D (optimal)  CutOff dis.     Sigma"]):
                        #rfn = os.path.join(data_dir, 'refdata', reffn + ref_type)
                        rfn = os.path.join(data_dir, reffn + ref_type)
                        if not os.path.isfile(rfn):
                            print('ERROR. Failed locating reference result file {:}'.format(rfn), file=sys.stderr);
                            sys.exit(6)
                        if not fileql('{:}'.format(ref_type[1:]), rfn, start_at):
                            print('{:140s} Failed@{:}'.format(' '.join([exe]+cmd), ref_type[1:]))
                        else:
                            print('{:140s} Passed (@{:})'.format(' '.join([exe]+cmd), ref_type[1:]))

#!/bin/python

import sys;
import getopt;
import os;
from subprocess import call;
import signal;

def sig_handler(signum, frame):
    print "segfault"

signal.signal(signal.SIGSEGV, sig_handler)

readFile = ''
trueFile = ''
baseDist = '100'
minOverlap = '0.5'
try:
	opts, args = getopt.getopt(sys.argv[1:],"hr:t:d:m:",["read=", "true=", "dist=", "overlap="])
except getopt.GetoptError:
	print 'USAGE: run_evaluate.py -r reads.fa -t true.m5 [-d baseDist -m minOverlap] label1,mapping1.m5 [label2,mapping2.m5 ...]'
	sys.exit(2)

# print opts;
# print args;

for opt, arg in opts:
	if opt == '-h':
		print 'USAGE: run_evaluate.py -r reads.fa -t true.m5 [-d baseDist -m minOverlap] label1,mapping1.m5 [label2,mapping2.m5 ...]'
		sys.exit()
	elif opt in ("-r", "--read"):
		readFile = arg;
	elif opt in ("-t", "--true"):
		trueFile = arg;
	elif opt in ("-d", "--dist"):
		baseDist = arg;
	elif opt in ("-m", "--overlap"):
		minOverlap = arg;

if readFile == '':
	print 'Please enter the read file!'
	print 'USAGE: run_evaluate.py -r reads.fa -t true.m5 [-d baseDist -m minOverlap] label1,mapping1.m5 [label2,mapping2.m5 ...]'
	sys.exit()
if trueFile == '':
	print 'Please enter the true mapping file!'
	print 'USAGE: run_evaluate.py -r reads.fa -t true.m5 [-d baseDist -m minOverlap] label1,mapping1.m5 [label2,mapping2.m5 ...]'
	sys.exit()
if len(args) == 0:
	print 'Please enter the reported mapping file'
	print 'USAGE: run_evaluate.py -r reads.fa -t true.m5 [-d baseDist -m minOverlap] label1,mapping1.m5 [label2,mapping2.m5 ...]'
	sys.exit()

# print 'Input read file: ', readFile
# print 'Input true mapping file: ', trueFile

evaluatePath = os.path.dirname(os.path.realpath(__file__)) + '/evaluate';

sys.stdout.write('          & # r     & # rMap  & # rCor  & # b        & # bCor     & # bIncor   & # bUnmap   & % sens & % prec ');
sys.stdout.write('\n');
sys.stdout.flush();

for m in args:
	[l, p] = m.split(',');
	sys.stdout.write(l.ljust(9) + ' & ');
	sys.stdout.flush();
	exitCode = call([evaluatePath, readFile, trueFile, p, baseDist, minOverlap]);
	if exitCode != 0:
		print '...';
	# print evaluatePath;
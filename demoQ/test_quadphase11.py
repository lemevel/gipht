#!/usr/bin/env python3


import sys
import os
import glob
sname = sys.argv[0]
verbose = "v"

#---------------------------------------------------


def check_file(file):
	cwd = os.getcwd()
	if not os.path.exists(file):
		print("cannot find file %s in %s", file, cwd)
		sys.exit(1)

	return(0)
#--------------------------------------------------
def link_file(file, path, sname, verbose):
	check_file(path)
	if not (os.path.exists(file)):
		run_cmd(('ln -s ' + path + ' .'), sname, verbose)

	return(0)
#--------------------------------------------------
def copy_file(file, path, sname, verbose):
	check_file(path)
	if not (os.path.exists(file)):
		run_cmd(('cp ' + path + ' .'), sname, verbose)

	return(0)
#--------------------------------------------------
def create_dir(dir, sname, verbose):
	if not os.path.exists(dir):
		run_cmd(('mkdir ' + dir), sname, verbose)

	return 0
#--------------------------------------------------
def find_max_min(value, min, max):

	if (value < min):
		min = value

	if (value > max):
		max = value

	return min, max
#--------------------------------------------------
def print_help(name):

	if os.path.exists((name + ".txt")):
		f = open((name + ".txt"), "r")
		for line in f.readlines():
			print(line)
		f.close()


	sys.exit(1)

	return 0
#------------------------------------------------
def run_cmd(cmd,sname,verbose):

	(dir,name) = os.path.split(sname)
	if (verbose == "v"):
		print (sname, ":", cmd)
	os.system(cmd)
	cwd = os.getcwd()

# now log it to current working directory
	f = open((cwd + "/" + name + ".log"),"a")
	f.write((cmd + "\n"))
	f.close()

	return(0)

#--------------------------------------------------
if (len(sys.argv) == 2):
	phafile = sys.argv[1]
	print("reading file",phafile, "\n")
	ierr = check_file(phafile)
	print("ierr is now ",ierr,"\n")

	# srcname = '../src/pha2qls.c';
	# exeext  = mexext;
	# exename = strrep(srcname,'.c',sprintf('.%s',exeext(4:end)))
	exename = "../src/pha2qls.maci64"
	# number of columns in each interferogram
	ncols =  1420
	# number of lines in each interferogram
	nrows = 1230
 
 	# east component of gradient
	grxfile = phafile.replace('.pha','.i2')
	grxfile = grxfile.replace('pha','grx')

	# north component of gradient
	gryfile = phafile.replace('.pha','.i2')
	gryfile = gryfile.replace('pha','gry')

	#cmd1 = (exename + ' pha_11176_21540_ort.pha  1420 1230 -P qha_11176_21540_ort.pha -X grx_11176_21540_ort.i2 -Y gry_11176_21540_ort.i2 -L 16 -M 8 -N 4')
	#cmd1 = "%s %s %d %d %s" % (exename, 'pha_11176_21540_ort.pha',ncols, nrows,' -P qha_11176_21540_ort.pha -X grx_11176_21540_ort.i2 -Y gry_11176_21540_ort.i2 -L 16 -M 8 -N 4')
	cmd1 = "%s %s %d %d -P %s -X %s -Y %s %s" % (exename, phafile, ncols, nrows\
		,phafile.replace('pha_','qha_'), grxfile, gryfile \
		,'-L 16 -M 8 -N 4')
	print(cmd1)
	run_cmd(cmd1, sname, verbose)
else:
	print("Usage: test_quadphase11\n")

sys.exit(1)

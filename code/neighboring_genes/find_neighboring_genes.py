#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import modules
import os
import sys
import argparse
import pybedtools as pybed

def main():

	parser = argparse.ArgumentParser(description = "Find neighboring gene pairs within a GFF file")
	parser.add_argument("-g", "--gff", dest = "gff", help = "GFF3 formatted annotation file", required = True)
	parser.add_argument("-p", "--out-prefix", dest = "prefix", help = "Output prefix to use for each output file", required = True)
	parser.add_argument("-i", "--ignore", dest = "ignore", help = "List of genes to ignore", required = True)

	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	args = parser.parse_args()

	# read in ignore list
	genes_to_ignore = dict()
	with open(args.ignore) as ignore_list:
		for line in ignore_list:
			line = line.rstrip()
			genes_to_ignore[line] = 1

	# create BedTool object of GFF3 file
	anno = pybed.BedTool(args.gff)
	# create dictionaries to keep track of every pairing
	divergent  = dict()
	convergent = dict()
	tandem     = dict()

	# go through each interval one at a time
	length = len(anno)
	for i in range(length):
		print("Working through interval",i,file=sys.stdout)
		interval = anno[i]
		# remove werid PlasmoDB gene IDs
		interval_name = interval.name.replace(".1","").replace(".2","").replace(".3","").replace(".4.","")
		left = "left"
		right = "right"
		if i > 0:
			left = anno[i - 1]
		if i < length - 1:
			right = anno[i + 1]
		if interval_name in genes_to_ignore.keys():
			continue
		if isinstance(left, pybed.cbedtools.Interval):
			# remove weird PlasmoDB gene IDs
			left_name = left.name.replace(".1","").replace(".2","").replace(".3","").replace(".4","")
			if left_name in genes_to_ignore.keys():
				continue
			if interval.chrom == left.chrom:
				if interval.strand == "+":
					if left.strand == "+":
						dist = interval.start - left.stop
						if left_name + "-" + interval_name not in tandem.keys() and interval_name + "-" + left_name not in tandem.keys():
							tandem[left_name + "-" + interval_name] = {"left" : left_name, "right" : interval_name, "dist" : dist}
					if left.strand == "-":
						dist = interval.start - left.stop
						if left_name + "-" + interval_name not in divergent.keys() and interval_name + "-" + left_name not in divergent.keys():
							divergent[left_name + "-" + interval_name] = {"left" : left_name, "right" : interval_name, "dist" : dist}
				if interval.strand == "-":
					if left.strand == "+":
						dist = interval.start - left.stop
						if left_name + "-" + interval_name not in convergent.keys() and interval_name + "-" + left_name not in convergent.keys():
							convergent[left_name + "-" + interval_name] = {"left" : left_name, "right" : interval_name, "dist" : dist}
					if left.strand == "-":
						dist = interval.start - left.stop
						if left_name + "-" + interval_name not in tandem.keys() and interval_name + "-" + left_name not in tandem.keys():
							tandem[left_name + "-" + interval_name] = {"left" : left_name, "right" : interval_name, "dist" : dist}
		if isinstance(right, pybed.cbedtools.Interval):
			# remove PlasmoDB faulty gene IDs
			right_name = right.name.replace(".1","").replace(".2","").replace(".3","").replace(".4","")
			if right_name in genes_to_ignore.keys():
				continue
			if interval.chrom == right.chrom:
				if interval.strand == "+":
					if right.strand == "+":
						dist = right.start - interval.stop
						if interval_name + "-" + right_name not in tandem.keys() and right_name + "-" + interval_name not in tandem.keys():
							tandem[interval_name + "-" + right_name] = {"left" : interval_name, "right" : right_name, "dist" : dist}
					if right.strand == "-":
						dist = right.start - interval.stop
						if interval_name + "-" + right_name not in convergent.keys() and right_name + "-" + interval_name not in convergent.keys():
							convergent[interval_name + "-" + right_name] = {"left" : interval_name, "right" : right_name, "dist" : dist}
				if interval.strand == "-":
					if right.strand == "+":
						dist = right.start - interval.stop
						if right_name + "-" + interval_name not in divergent.keys() and right_name + "-" + interval_name not in divergent.keys():
							divergent[interval_name + "-" + right_name] = {"left" : interval_name, "right" : right_name, "dist" : dist}
					if right.strand == "-":
						dist = right.start - interval.stop
						if right_name + "-" + interval_name not in tandem.keys() and right_name + "-" + interval_name not in tandem.keys():
							tandem[interval_name + "-" + right_name] = {"left" : interval_name, "right" : right_name, "dist" : dist}

	# output results to file
	f = open(args.prefix + "_divergent.tsv", "w")
	print("left_gene\tright_gene\tdist",file=f)
	for key in divergent.keys():
		print("{0}\t{1}\t{2}".format(divergent[key]["left"],divergent[key]["right"],str(divergent[key]["dist"])),file=f)
	f.close()

	f = open(args.prefix + "_convergent.tsv", "w")
	print("left_gene\tright_gene\tdist",file=f)
	for key in convergent.keys():
		print("{0}\t{1}\t{2}".format(convergent[key]["left"], convergent[key]["right"], str(convergent[key]["dist"])),file=f)
	f.close()

	f = open(args.prefix + "_tandem.tsv", "w")
	print("left_gene\tright_gene\tdist",file=f)
	for key in tandem.keys():
		print("{0}\t{1}\t{2}".format(tandem[key]["left"], tandem[key]["right"], str(tandem[key]["dist"])),file=f)
	f.close()

if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("User interrupt! Ciao!\n")
		sys.exit(0)

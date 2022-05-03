#!/usr/bin/env python
# Written by mhchoi

import sys, argparse, math, pysam, re, pickle, os
import numpy as np
import queue
import threading
from multiprocessing import Process, Queue

# Cigar information
CIGAR_MAP = {"M":0, "I":1, "D":2, "N":3, "S":4, "H":5, "P":6, "=":7, "X":8, "B":9}

def vcf_to_Dic(input_VCF, target_type):
	varDic = dict()
	with open(input_VCF, 'r') as f:
		for lines in f:
			# The SV record line
			if not lines.startswith('#'):
				items = lines.strip().split('\t')
				svType = [a.split('=')[-1] for a in items[7].split(';') if 'SVTYPE=' in a][0]
				chr2_n = [a.split('=')[-1] for a in items[7].split(';') if 'CHR2=' in a][0]
				chr2_p = [a.split('=')[-1] for a in items[7].split(';') if 'END=' in a][0]

				if target_type == svType:
					# Save format: {ChrA: set{ChrA\tPosA\tChrB\tPosB}}
					if not items[0] in varDic:
						varDic.setdefault(items[0], set()) \
							  .add('{0}\t{1}\t{2}\t{3}'.format(items[0], items[1], chr2_n, chr2_p))
					else: varDic[items[0]].add('{0}\t{1}\t{2}\t{3}'.format(items[0], items[1], chr2_n, chr2_p))

				del items; del svType; del chr2_n; del chr2_p
				
			# The header line
			else: pass

	return varDic

class genome_interval:
	def __init__(self, chrm, start, end):
		self.chrm = chrm
		self.start = start
		self.end = end

	def __str__(self):
		return "(" + self.chrm + "," + str(self.start) + "," + str(self.end) + ")"

	def __repr__(self):
		return str(self)

	def __eq__(self, gi2):
		return self.chrm == gi2.chrm and self.start == gi2.start and self.end == gi2.end
		""" return -1 if before, 0 if in, 1 if after """

	def intersect(self, gi):
		if gi.chrm.replace("chr", "") < self.chrm.replace("chr", "") or gi.end < self.start:
			return -1
		elif gi.chrm.replace("chr", "") > self.chrm.replace("chr", "") or gi.start > self.end:
			return 1
		else:
			return 0

##Long Read methods
class Alignment:
	def __init__(self, chrm, start, end, strand, query_position):
		self.pos = genome_interval(chrm, start, end)
		self.strand = strand
		self.query_position = query_position

	def __str__(self):
		return ",".join([str(x)for x in [self.pos.chrm, self.pos.start, self.pos.end, self.strand, self.query_position,]])

	def __repr__(self):
		return "Alignment(%s,%s,%s,%s,%s)" % (self.pos.chrm, self.pos.start, self.pos.end, self.strand, self.query_position,)

class LongRead:
    def __init__(self, alignments):
        self.alignments = alignments

    def __str__(self):
        return ",".join([str(x) for x in self.alignments])

    def __repr__(self):
        return "LongRead(" + str(self) + ")"

class plan_step:
	step_events = ["Align", "ANNOTATION"]

	def __init__(self, start_pos, end_pos, event, info=None):
		self.start_pos = start_pos
		self.end_pos = end_pos
		self.event = event
		self.info = info

	def __str__(self):
		if self.info:
			return ("Step(" + str(self.start_pos) + ", " + str(self.end_pos) + ", " + self.event + ", " + str(self.info) + ")")
		else:
			return ("Step(" + str(self.start_pos) + ", " + str(self.end_pos) + ", " + self.event + ")")
			
	def __repr__(self):
		return str(self)

def estimate_fragment_len(long_bam, reference_seq):
	try:
		if not reference_seq:
			bam_file = pysam.AlignmentFile(long_bam, "rb")
		else:
			bam_file = pysam.AlignmentFile(long_bam, "rc", reference_filename=reference_seq)
	except Exception as err:
		print("Error:", err, file=sys.stderr)
		sys.exit(1)
		
	frag_lens = []
	for i, read in enumerate(bam_file):
		if i >= 10000:
			break
		frag_lens.append(abs(read.tlen))

	if len(frag_lens) >= 5000:
		return np.median(frag_lens)
	else:
		print("Insufficient reads for fragment length estimate.\nContinuing with unmodified window size", file=sys.stderr,)
		return 0

def get_alignments_from_cigar(chrm, curr_pos, strand, cigartuples, reverse=False):
	alignments = []
	q_pos = 0
	if reverse: cigartuples = cigartuples[::-1]

	for op, length in cigartuples:
		if op in [CIGAR_MAP["M"], CIGAR_MAP["="], CIGAR_MAP["X"]]:
			alignments.append(Alignment(chrm, curr_pos, curr_pos + length, strand, q_pos))
			curr_pos += length
			q_pos += length
		elif op == CIGAR_MAP["I"]: q_pos += length
		elif op == CIGAR_MAP["D"]: curr_pos += length
		elif op == CIGAR_MAP["N"]: curr_pos += length
		elif op == CIGAR_MAP["S"]: q_pos += length
	
	return alignments

def get_cigartuples_from_string(cigarstring):
	cigartuples = []
	
	for match in re.findall(r"(\d+)([A-Z]{1})", cigarstring):
		length = int(match[0])
		op = match[1]
		cigartuples.append((CIGAR_MAP[op], length))
		
	return cigartuples

def merge_alignments(min_gap, alignments):
	merged_alignments = []

	for alignment in alignments:
		if len(merged_alignments) == 0:
			merged_alignments.append(alignment)
		else:
			if (alignment.pos.chrm == merged_alignments[-1].pos.chrm \
				and alignment.pos.start < merged_alignments[-1].pos.end + min_gap):
				merged_alignments[-1].pos.end = alignment.pos.end
			else: merged_alignments.append(alignment)

	return merged_alignments

def get_range_hit(ranges, chrm, point):
	for j in range(len(ranges)):
		r = ranges[j]
		if (r.chrm.replace("chr", "") == chrm.replace("chr", "") and r.start <= point and r.end >= point):
			return j
	return None

def add_align_step(alignment, steps, ranges):
	# alignment can span ranges
	start_range_hit_i = get_range_hit(ranges, alignment.pos.chrm, alignment.pos.start)
	end_range_hit_i = get_range_hit(ranges, alignment.pos.chrm, alignment.pos.end)
	
	# neither end is in range, add nothing
	if start_range_hit_i == None and end_range_hit_i == None:
		return

    # start is not in range, use end hit
	if start_range_hit_i == None:
		start = genome_interval(alignment.pos.chrm, max(alignment.pos.start, ranges[end_range_hit_i].start), \
								max(alignment.pos.start, ranges[end_range_hit_i].start),)
		end = genome_interval(alignment.pos.chrm, min(alignment.pos.end, ranges[end_range_hit_i].end), \
							  min(alignment.pos.end, ranges[end_range_hit_i].end),)
		steps.append(plan_step(start, end, "LONGREAD", info={"TYPE": "Align"}))

	# end is not in range, use start hit
	elif end_range_hit_i == None:
		start = genome_interval(alignment.pos.chrm, max(alignment.pos.start, ranges[start_range_hit_i].start), \
								max(alignment.pos.start, ranges[start_range_hit_i].start),)
		end = genome_interval(alignment.pos.chrm, min(alignment.pos.end, ranges[start_range_hit_i].end), \
							  min(alignment.pos.end, ranges[start_range_hit_i].end),)
		steps.append(plan_step(start, end, "LONGREAD", info={"TYPE": "Align"}))
		
	# both are in the same range
	elif start_range_hit_i == end_range_hit_i:
		start = genome_interval(alignment.pos.chrm, max(alignment.pos.start, ranges[start_range_hit_i].start), \
								max(alignment.pos.start, ranges[start_range_hit_i].start),)
		end = genome_interval(alignment.pos.chrm, min(alignment.pos.end, ranges[end_range_hit_i].end), \
							  min(alignment.pos.end, ranges[end_range_hit_i].end),)
		steps.append(plan_step(start, end, "LONGREAD", info={"TYPE": "Align"}))

	# in different ranges
	else:
		start_1 = genome_interval(alignment.pos.chrm, max(alignment.pos.start, ranges[start_range_hit_i].start), \
								  max(alignment.pos.start, ranges[start_range_hit_i].start),)
		end_1 = genome_interval(alignment.pos.chrm, ranges[start_range_hit_i].end, ranges[start_range_hit_i].end,)
		steps.append(plan_step(start_1, end_1, "LONGREAD", info={"TYPE": "Align"}))
		
		start_2 = genome_interval(alignment.pos.chrm, ranges[end_range_hit_i].start, ranges[end_range_hit_i].start,)
		end_2 = genome_interval(alignment.pos.chrm, min(alignment.pos.end, ranges[end_range_hit_i].end), \
								min(alignment.pos.end, ranges[end_range_hit_i].end),)
		steps.append(plan_step(start_2, end_2, "LONGREAD", info={"TYPE": "Align"}))

	return steps

def SV_judgement_runner(bam_file, ranges, svType, target_range, sv_list, q):
	# Fetch the target region
	long_reads_Dic = {}
	for r in ranges:
		## Extend target ragnes: +- 1000bp 
		try:
			bam_iter = bam_file.fetch(r.chrm, max(0, r.start - 1000), r.end + 1000)
		except ValueError:
			chrm = r.chrm
			if chrm[:3] == "chr" or chrm[:3] == "Chr":
				chrm = chrm[3:]
			else:
				chrm = "chr" + chrm
			bam_iter = bam_file.fetch(chrm, max(0, r.start - 1000), r.end + 1000)

		## Get the alignment information in the target region
		for read in bam_iter:
			## Pass the QC fail OR Unmapped OR Duplicate OR less than Min.MQ 0
			if read.is_qcfail or read.is_unmapped or \
				read.is_duplicate or int(read.mapping_quality) == 0:
				continue

			alignments = get_alignments_from_cigar(
									bam_file.get_reference_name(read.reference_id), read.pos, 
									not read.is_reverse, read.cigartuples,)

			merged_alignments = merge_alignments(50, alignments)

			read_strand = not read.is_reverse

			## Check supplymentary alignment information
			if read.has_tag("SA"):
				for sa in read.get_tag("SA").split(";"):
					if len(sa) == 0: continue
					
					rname, pos, strand, cigar, mapq, nm = sa.split(",")
					sa_pos = int(pos)
					sa_strand = strand == "+"
					strand_match = read_strand != sa_strand
					sa_cigartuples = get_cigartuples_from_string(cigar)
					sa_alignments = get_alignments_from_cigar(rname, sa_pos, sa_strand, sa_cigartuples, reverse=strand_match)
					sa_merged_alignments = merge_alignments(50, sa_alignments)

					if len(sa_merged_alignments) > 0:
						merged_alignments += sa_merged_alignments

			## Save alignment results of long-read bam
			if read.query_name not in long_reads_Dic:
				long_reads_Dic.setdefault(read.query_name, []).append(LongRead(merged_alignments))
			else:					
				long_reads_Dic[read.query_name].append(LongRead(merged_alignments))
			
			del alignments; del merged_alignments; del read_strand

		del bam_iter

	val_Set = validate_SVs_in_short_read(long_reads_Dic, ranges, svType, target_range, sv_list)

	q.put(val_Set)

def get_long_read_plan(read_name, long_reads_Dic, ranges):
	alignments = []
	
	# only keep alignments that intersect a range
	seen = {}
	
	if read_name not in long_reads_Dic:
		sys.stderr.write("ERROR: Read name " + read_name + " not in list of long reads")
		sys.exit(1)
		
	for long_read in long_reads_Dic[read_name]:
		for alignment in long_read.alignments:
			if alignment.query_position in seen:
				continue
			seen[alignment.query_position] = 1

			in_range = False
			for r in ranges:
				if r.intersect(alignment.pos) == 0:
					in_range = True
			if in_range:
				alignments.append(alignment)
				
	if len(alignments) <= 0:
		return None

	alignments.sort(key=lambda x: x.query_position)

	longest_alignment = 0
	longest_alignment_i = -1
	for i in range(len(alignments)):
		l = alignments[i].pos.end - alignments[i].pos.start
		if longest_alignment < l:
			longest_alignment = l
			longest_alignment_i = i
	primary_strand = alignments[longest_alignment_i].strand

	steps = []
	curr = alignments[0]
	
	add_align_step(curr, steps, ranges)
	
	for i in range(1, len(alignments)):
		last = alignments[i - 1]
		curr = alignments[i]

		##############################
		# figure out what the event is
		## Event: INTER CHROM
		if curr.pos.chrm != last.pos.chrm:
			if curr.strand != last.strand:
				start = genome_interval(last.pos.chrm, last.pos.end, last.pos.end)
				end = genome_interval(curr.pos.chrm, curr.pos.end, curr.pos.end)
				info = {"TYPE": "TRA"}
				steps.append(plan_step(start, end, "LONGREAD", info=info))
			else:
				start = genome_interval(last.pos.chrm, last.pos.end, last.pos.end)
				end = genome_interval(curr.pos.chrm, curr.pos.start, curr.pos.start)
				info = {"TYPE": "TRA"}
				steps.append(plan_step(start, end, "LONGREAD", info=info))
			
			add_align_step(curr, steps, ranges)

		## Inversion
		elif curr.strand != last.strand:
			# it is possible that we have a complex even that
			# is an inverted Duplication
			if curr.pos.start < last.pos.end:
				start = genome_interval(last.pos.chrm, last.pos.end, last.pos.end)
				end = genome_interval(curr.pos.chrm, curr.pos.start, curr.pos.start)
				info = {"TYPE": "DEL"}
				steps.append(plan_step(start, end, "LONGREAD", info=info))
			
			if curr.strand != primary_strand:
				start = genome_interval(last.pos.chrm, last.pos.end, last.pos.end)
				end = genome_interval(curr.pos.chrm, curr.pos.end, curr.pos.end)
				info = {"TYPE": "INV"}
				steps.append(plan_step(start, end, "LONGREAD", info=info))
			else:

				if curr.pos.start < last.pos.end:
					start = genome_interval(last.pos.chrm, last.pos.end, last.pos.end)
					end = genome_interval(curr.pos.chrm, curr.pos.start, curr.pos.start)
					info = {"TYPE": "DUP"}
					steps.append(plan_step(start, end, "LONGREAD", info=info))

				start = genome_interval(last.pos.chrm, last.pos.start, last.pos.start)
				end = genome_interval(curr.pos.chrm, curr.pos.start, curr.pos.start)
				info = {"TYPE": "INV"}
				steps.append(plan_step(start, end, "LONGREAD", info=info))

			add_align_step(curr, steps, ranges)

		## Duplication
		elif curr.pos.start < last.pos.end:
			start = genome_interval(last.pos.chrm, last.pos.end, last.pos.end)
			end = genome_interval(curr.pos.chrm, curr.pos.start, curr.pos.start)
			info = {"TYPE": "DUP"}
			steps.append(plan_step(start, end, "LONGREAD", info=info))
			add_align_step(curr, steps, ranges)

		## Deletion
		else:
			start = genome_interval(last.pos.chrm, last.pos.end, last.pos.end)
			end = genome_interval(curr.pos.chrm, curr.pos.start, curr.pos.start)
			info = {"TYPE": "DEL"}
			steps.append(plan_step(start, end, 'LONGREAD', info=info))
			add_align_step(curr, steps, ranges)

		# if either end is in a range, then add its gap to the list

	max_gap = None

	chrms = set([s.start_pos.chrm for s in steps] + [s.end_pos.chrm for s in steps])

	# set interchrm dist to 5000
	if len(chrms) > 1:
		max_gap = 5000
	else:
		step_sizes = [abs(step.end_pos.end - step.start_pos.start) \
					  for step in steps if step.info["TYPE"] != "Align" \
					  and get_range_hit(ranges, step.start_pos.chrm, step.start_pos.start) != None \
					  and get_range_hit(ranges, step.end_pos.chrm, step.end_pos.end) != None]

		max_gap = max(step_sizes) if len(step_sizes) > 0 else 0
		
	plan = [max_gap, steps]
	
	return plan

def validate_SVs_in_short_read(long_reads_Dic, ranges, svType, target_range, sv_list):
	# Set up target SV position
	# DEL, DUP, and INV
	if sv_list[0] == sv_list[2]:
		if int(sv_list[1]) < int(sv_list[3]):
			sv_c1 = sv_list[0]
			sv_c2 = sv_list[2]
			sv_s = int(sv_list[1])
			sv_e = int(sv_list[3])
			sv_ls = (sv_e - sv_s + 1) - math.ceil((sv_e - sv_s + 1) * 0.25)
			sv_lr = (sv_e - sv_s + 1) + math.ceil((sv_e - sv_s + 1) * 0.25)
		else:
			sv_c1 = sv_list[2]
			sv_c2 = sv_list[0]
			sv_s = int(sv_list[3])
			sv_e = int(sv_list[1])
			sv_ls = (sv_e - sv_s + 1) - math.ceil((sv_e - sv_s + 1) * 0.25)
			sv_lr = (sv_e - sv_s + 1) + math.ceil((sv_e - sv_s + 1) * 0.25)
	# TRA
	else:
		sv_c1 = sv_list[0]
		sv_c2 = sv_list[2]
		sv_s = int(sv_list[1])
		sv_e = int(sv_list[3])

	# Check the target SV position & length (only DEL & INV)
	val_Set = set()

	for read_name in long_reads_Dic.keys():
		long_read_plan = get_long_read_plan(read_name, long_reads_Dic, ranges)

		if long_read_plan is None:
			continue
		else:
			max_gap = long_read_plan[0]
			steps = long_read_plan[1]

			if len(steps) > 0:
				for step in steps:
					val_string = ''

					# DEL, DUP, and INV
					if svType != 'TRA':
						# Check Start vs BP-1
						if str(step.start_pos.chrm) == sv_c1:
							# Tier-1: +-5bp near BP-1
							if (max(0, sv_s-5) <= int(step.start_pos.start) and sv_s+5 >= int(step.start_pos.start)) or \
								(max(0, sv_s-5) <= int(step.start_pos.end) and sv_s+5 >= int(step.start_pos.end)):
								val_string += 'T1'
							# Tier-2: +-10bp near BP-1
							elif (max(0, sv_s-10) <= int(step.start_pos.start) and sv_s+10 >= int(step.start_pos.start)) or \
								(max(0, sv_s-10) <= int(step.start_pos.end) and sv_s+10 >= int(step.start_pos.end)):
								val_string += 'T2'
							# Tier-3: +-25bp near BP-1
							elif (max(0, sv_s-25) <= int(step.start_pos.start) and sv_s+25 >= int(step.start_pos.start)) or \
								(max(0, sv_s-25) <= int(step.start_pos.end) and sv_s+25 >= int(step.start_pos.end)):
								val_string += 'T3'
							# Tier-4: +-50bp near BP-1
							elif (max(0, sv_s-50) <= int(step.start_pos.start) and sv_s+50 >= int(step.start_pos.start)) or \
								(max(0, sv_s-50) <= int(step.start_pos.end) and sv_s+50 >= int(step.start_pos.end)):
								val_string += 'T4'
							# Tier-5: specific range near BP-1
							elif (max(0, sv_s-target_range) <= int(step.start_pos.start) and sv_s+target_range >= int(step.start_pos.start)) or \
								(max(0, sv_s-target_range) <= int(step.start_pos.end) and sv_s+target_range >= int(step.start_pos.end)):
								val_string += 'T5'
							else:
								val_string += 'F'
						else:
							val_string += 'F'

						# Check End vs BP-2
						if str(step.end_pos.chrm) == sv_c2:
							# Tier-1: +-5bp near BP-2
							if (max(0, sv_e-5) <= int(step.end_pos.start) and sv_e+5 >= int(step.end_pos.start)) or \
								(max(0, sv_e-5) <= int(step.end_pos.end) and sv_e+5 >= int(step.end_pos.end)):
								val_string += '_T1'
							# Tier-2: +-10bp near BP-2
							elif (max(0, sv_e-10) <= int(step.end_pos.start) and sv_e+10 >= int(step.end_pos.start)) or \
								(max(0, sv_e-10) <= int(step.end_pos.end) and sv_e+10 >= int(step.end_pos.end)):
								val_string += '_T2'
							# Tier-3: +-25bp near BP-2
							elif (max(0, sv_e-25) <= int(step.end_pos.start) and sv_e+25 >= int(step.end_pos.start)) or \
								(max(0, sv_e-25) <= int(step.end_pos.end) and sv_e+25 >= int(step.end_pos.end)):
								val_string += '_T3'
							# Tier-4: +-50bp near BP-2
							elif (max(0, sv_e-50) <= int(step.end_pos.start) and sv_e+50 >= int(step.end_pos.start)) or \
								(max(0, sv_e-50) <= int(step.end_pos.end) and sv_e+50 >= int(step.end_pos.end)):
								val_string += '_T4'
							# Tier-5: specific ragne near BP-2
							elif (max(0, sv_e-target_range) <= int(step.end_pos.start) and sv_e+target_range >= int(step.end_pos.start)) or \
								(max(0, sv_e-target_range) <= int(step.end_pos.end) and sv_e+target_range >= int(step.end_pos.end)):
								val_string += '_T5'
							else:
								val_string += '_F'
						else:
							val_string += '_F'

						## Check the SV length
						### Allow 25%
						### 75bp <= 100bp <= 125bp
						if (svType == 'DEL') or (svType == 'INV'):
							if 'T' in val_string:
								if sv_ls <= int(step.end_pos.start)-int(step.start_pos.end)+1 and \
									int(step.end_pos.start)-int(step.start_pos.end)+1 <= sv_lr:
									pass
								else:
									val_string = 'F_F'
							else:
								pass
						else: pass
						
						## Check the SV type
						if svType == step.info['TYPE']:
							val_string += '_T'
						### False SV type
						else:
							val_string += '_{0}'.format(step.info['TYPE'])

						val_Set.add(val_string)
					
					# TRA type
					else:
						# Check Start vs BP-1
						if str(step.start_pos.chrm) == sv_c1:
							# Tier-1: +-5bp near BP-1
							if (max(0, sv_s-5) <= int(step.start_pos.start) and sv_s+5 >= int(step.start_pos.start)) or \
								(max(0, sv_s-5) <= int(step.start_pos.end) and sv_s+5 >= int(step.start_pos.end)):
								val_string += 'T1'
							# Tier-2: +-10bp near BP-1
							elif (max(0, sv_s-10) <= int(step.start_pos.start) and sv_s+10 >= int(step.start_pos.start)) or \
								(max(0, sv_s-10) <= int(step.start_pos.end) and sv_s+10 >= int(step.start_pos.end)):
								val_string += 'T2'
							# Tier-3: +-25bp near BP-1
							elif (max(0, sv_s-25) <= int(step.start_pos.start) and sv_s+25 >= int(step.start_pos.start)) or \
								(max(0, sv_s-25) <= int(step.start_pos.end) and sv_s+25 >= int(step.start_pos.end)):
								val_string += 'T3'
							# Tier-4: +-50bp near BP-1
							elif (max(0, sv_s-50) <= int(step.start_pos.start) and sv_s+50 >= int(step.start_pos.start)) or \
								(max(0, sv_s-50) <= int(step.start_pos.end) and sv_s+50 >= int(step.start_pos.end)):
								val_string += 'T4'
							# Tier-5: specific range near BP-1
							elif (max(0, sv_s-target_range) <= int(step.start_pos.start) and sv_s+target_range >= int(step.start_pos.start)) or \
								(max(0, sv_s-target_range) <= int(step.start_pos.end) and sv_s+target_range >= int(step.start_pos.end)):
								val_string += 'T5'
							else:
								val_string += 'F'
						else:
							val_string += 'F'

						# Check End vs BP-2
						if str(step.end_pos.chrm) == sv_c2:
							# Tier-1: +-5bp near BP-2
							if (max(0, sv_e-5) <= int(step.end_pos.start) and sv_e+5 >= int(step.end_pos.start)) or \
								(max(0, sv_e-5) <= int(step.end_pos.end) and sv_e+5 >= int(step.end_pos.end)):
								val_string += '_T1'
							# Tier-2: +-10bp near BP-2
							elif (max(0, sv_e-10) <= int(step.end_pos.start) and sv_e+10 >= int(step.end_pos.start)) or \
								(max(0, sv_e-10) <= int(step.end_pos.end) and sv_e+10 >= int(step.end_pos.end)):
								val_string += '_T2'
							# Tier-3: +-25bp near BP-2
							elif (max(0, sv_e-25) <= int(step.end_pos.start) and sv_e+25 >= int(step.end_pos.start)) or \
								(max(0, sv_e-25) <= int(step.end_pos.end) and sv_e+25 >= int(step.end_pos.end)):
								val_string += '_T3'
							# Tier-4: +-50bp near BP-2
							elif (max(0, sv_e-50) <= int(step.end_pos.start) and sv_e+50 >= int(step.end_pos.start)) or \
								(max(0, sv_e-50) <= int(step.end_pos.end) and sv_e+50 >= int(step.end_pos.end)):
								val_string += '_T4'
							# Tier-5: specific ragne near BP-2
							elif (max(0, sv_e-target_range) <= int(step.end_pos.start) and sv_e+target_range >= int(step.end_pos.start)) or \
								(max(0, sv_e-target_range) <= int(step.end_pos.end) and sv_e+target_range >= int(step.end_pos.end)):
								val_string += '_T5'
							else:
								val_string += '_F'
						else:
							val_string += '_F'

						# Reverse checking
						if  not 'T' in val_string:
							val_string = ''
							
							# Check End vs BP-1
							if str(step.end_pos.chrm) == sv_c1:
								# Tier-1: +-5bp near BP-1
								if (max(0, sv_s-5) <= int(step.end_pos.start) and sv_s+5 >= int(step.end_pos.start)) or \
									(max(0, sv_s-5) <= int(step.end_pos.end) and sv_s+5 >= int(step.end_pos.end)):
									val_string += 'T1'
								# Tier-2: +-10bp near BP-1
								elif (max(0, sv_s-10) <= int(step.end_pos.start) and sv_s+10 >= int(step.end_pos.start)) or \
									(max(0, sv_s-10) <= int(step.end_pos.end) and sv_s+10 >= int(step.end_pos.end)):
									val_string += 'T2'
								# Tier-3: +-25bp near BP-1
								elif (max(0, sv_s-25) <= int(step.end_pos.start) and sv_s+25 >= int(step.end_pos.start)) or \
									(max(0, sv_s-25) <= int(step.end_pos.end) and sv_s+25 >= int(step.end_pos.end)):
									val_string += 'T3'
								# Tier-4: +-50bp near BP-1
								elif (max(0, sv_s-50) <= int(step.end_pos.start) and sv_s+50 >= int(step.end_pos.start)) or \
									(max(0, sv_s-50) <= int(step.end_pos.end) and sv_s+50 >= int(step.end_pos.end)):
									val_string += 'T4'
								# Tier-5: specific range near BP-1
								elif (max(0, sv_s-target_range) <= int(step.end_pos.start) and sv_s+target_range >= int(step.end_pos.start)) or \
									(max(0, sv_s-target_range) <= int(step.end_pos.end) and sv_s+target_range >= int(step.end_pos.end)):
									val_string += 'T5'
								else:
									val_string += 'F'
							else:
								val_string += 'F'

							# Check Start vs BP-2
							if str(step.start_pos.chrm) == sv_c2:
								# Tier-1: +-5bp near BP-2
								if (max(0, sv_e-5) <= int(step.start_pos.start) and sv_e+5 >= int(step.start_pos.start)) or \
									(max(0, sv_e-5) <= int(step.start_pos.end) and sv_e+5 >= int(step.start_pos.end)):
									val_string += '_T1'
								# Tier-2: +-10bp near BP-2
								elif (max(0, sv_e-10) <= int(step.start_pos.start) and sv_e+10 >= int(step.start_pos.start)) or \
									(max(0, sv_e-10) <= int(step.start_pos.end) and sv_e+10 >= int(step.start_pos.end)):
									val_string += '_T2'
								# Tier-3: +-25bp near BP-2
								elif (max(0, sv_e-25) <= int(step.start_pos.start) and sv_e+25 >= int(step.start_pos.start)) or \
									(max(0, sv_e-25) <= int(step.start_pos.end) and sv_e+25 >= int(step.start_pos.end)):
									val_string += '_T3'
								# Tier-4: +-50bp near BP-2
								elif (max(0, sv_e-50) <= int(step.start_pos.start) and sv_e+50 >= int(step.start_pos.start)) or \
									(max(0, sv_e-50) <= int(step.start_pos.end) and sv_e+50 >= int(step.start_pos.end)):
									val_string += '_T4'
								# Tier-5: specific range near BP-2
								elif (max(0, sv_e-target_range) <= int(step.start_pos.start) and sv_e+target_range >= int(step.start_pos.start)) or \
									(max(0, sv_e-target_range) <= int(step.start_pos.end) and sv_e+target_range >= int(step.start_pos.end)):
									val_string += '_T5'
								else:
									val_string += '_F'
							else:
								val_string += '_F'
						else:
							pass

						## Check the SV type
						if svType == step.info['TYPE']:
							val_string += '_T'
						### False SV type
						else:
							val_string += '_{0}'.format(step.info['TYPE'])
			
						val_Set.add(val_string)

					del val_string

			else: pass

		del long_read_plan

	return val_Set

def pattern_checker(val_info_set, tag_info_set):
	for val_string in val_info_set:

		# Tier-1: Tier-1 and Tier-1
		if 'T1_T1' in val_string:
			if val_string == 'T1_T1_T':
				tag_info_set.add('T1')
			else:
				val_type = val_string.split('_')[-1]
				if val_type == 'Align':
					tag_info_set.add('T1_Align')
				else:
					tag_info_set.add('D1_{0}'.format(val_string.split('_')[-1]))
				del val_type

		# Tier-2: Tier-1 and Tier-2 ||  Tier-2 and Tier-1
		elif 'T1_T2' in val_string or \
		     'T2_T1' in val_string:
			if val_string == 'T1_T2_T' or \
			   val_string == 'T2_T1_T':
				tag_info_set.add('T2')
			else:
				val_type = val_string.split('_')[-1]
				if val_type == 'Align':
					tag_info_set.add('T2_Align')
				else:
					tag_info_set.add('D2_{0}'.format(val_string.split('_')[-1]))
				del val_type

		# Tier-2: Tier-2 and Tier-2
		elif 'T2_T2' in val_string:
			if val_string == 'T2_T2_T':
				tag_info_set.add('T2')
			else:
				val_type = val_string.split('_')[-1]
				if val_type == 'Align':
					tag_info_set.add('T2_Align')
				else:
					tag_info_set.add('D2_{0}'.format(val_string.split('_')[-1]))
				del val_type

		# Tier-3: Tier-1 and Tier-3 ||  Tier-3 and Tier-1
		elif 'T1_T3' in val_string or \
		     'T3_T1' in val_string:
			if val_string == 'T1_T3_T' or \
				val_string == 'T3_T1_T':
				tag_info_set.add('T3')
			else:
				val_type = val_string.split('_')[-1]
				if val_type == 'Align':
					tag_info_set.add('T3_Align')
				else:
					tag_info_set.add('D3_{0}'.format(val_string.split('_')[-1]))
				del val_type

		# Tier-3: Tier-2 and Tier-3 ||  Tier-3 and Tier-2
		elif 'T2_T3' in val_string or \
		     'T3_T2' in val_string:
			if val_string == 'T2_T3_T' or \
				val_string == 'T3_T2_T':
				tag_info_set.add('T3')
			else:
				val_type = val_string.split('_')[-1]
				if val_type == 'Align':
					tag_info_set.add('T3_Align')
				else:
					tag_info_set.add('D3_{0}'.format(val_string.split('_')[-1]))
				del val_type

		# Tier-3: Tier-3 and Tier-3
		elif 'T3_T3' in val_string:
			if val_string == 'T3_T3_T':
				tag_info_set.add('T3')
			else:
				val_type = val_string.split('_')[-1]
				if val_type == 'Align':
					tag_info_set.add('T3_Align')
				else:
					tag_info_set.add('D3_{0}'.format(val_string.split('_')[-1]))
				del val_type

		# Tier-4: Tier-1 and Tier-4 ||  Tier-4 and Tier-1
		elif 'T1_T4' in val_string or \
		     'T4_T1' in val_string:
			if val_string == 'T1_T4_T' or \
				val_string == 'T4_T1_T':
				tag_info_set.add('T4')
			else:
				val_type = val_string.split('_')[-1]
				if val_type == 'Align':
					tag_info_set.add('T4_Align')
				else:
					tag_info_set.add('D4_{0}'.format(val_string.split('_')[-1]))
				del val_type

		# Tier-4: Tier-2 and Tier-4 ||  Tier-4 and Tier-2
		elif 'T2_T4' in val_string or \
		     'T4_T2' in val_string:
			if val_string == 'T2_T4_T' or \
				val_string == 'T4_T2_T':
				tag_info_set.add('T4')
			else:
				val_type = val_string.split('_')[-1]
				if val_type == 'Align':
					tag_info_set.add('T4_Align')
				else:
					tag_info_set.add('D4_{0}'.format(val_string.split('_')[-1]))
				del val_type

		# Tier-4: Tier-3 and Tier-4 ||  Tier-4 and Tier-3
		elif 'T3_T4' in val_string or \
		     'T4_T3' in val_string:
			if val_string == 'T3_T4_T' or \
				val_string == 'T4_T3_T':
				tag_info_set.add('T4')
			else:
				val_type = val_string.split('_')[-1]
				if val_type == 'Align':
					tag_info_set.add('T4_Align')
				else:
					tag_info_set.add('D4_{0}'.format(val_string.split('_')[-1]))
				del val_type

		# Tier-4: Tier-4 and Tier-4
		elif 'T4_T4' in val_string:
			if val_string == 'T4_T4_T':
				tag_info_set.add('T4')
			else:
				val_type = val_string.split('_')[-1]
				if val_type == 'Align':
					tag_info_set.add('T4_Align')
				else:
					tag_info_set.add('D4_{0}'.format(val_string.split('_')[-1]))
				del val_type

		# Tier-5: Tier-1 and Tier-5 ||  Tier-5 and Tier-1
		elif 'T1_T5' in val_string or \
		     'T5_T1' in val_string:
			if val_string == 'T1_T5_T' or \
				val_string == 'T5_T1_T':
				tag_info_set.add('T5')
			else:
				val_type = val_string.split('_')[-1]
				if val_type == 'Align':
					tag_info_set.add('T5_Align')
				else:
					tag_info_set.add('D5_{0}'.format(val_string.split('_')[-1]))
				del val_type

		# Tier-5: Tier-2 and Tier-5 ||  Tier-5 and Tier-2
		elif 'T2_T5' in val_string or \
		     'T5_T2' in val_string:
			if val_string == 'T2_T5_T' or \
				val_string == 'T5_T2_T':
				tag_info_set.add('T5')
			else:
				val_type = val_string.split('_')[-1]
				if val_type == 'Align':
					tag_info_set.add('T5_Align')
				else:
					tag_info_set.add('D5_{0}'.format(val_string.split('_')[-1]))
				del val_type

		# Tier-5: Tier-3 and Tier-5 ||  Tier-5 and Tier-3
		elif 'T3_T5' in val_string or \
		     'T5_T3' in val_string:
			if val_string == 'T3_T5_T' or \
				val_string == 'T5_T3_T':
				tag_info_set.add('T5')
			else:
				val_type = val_string.split('_')[-1]
				if val_type == 'Align':
					tag_info_set.add('T5_Align')
				else:
					tag_info_set.add('D5_{0}'.format(val_string.split('_')[-1]))
				del val_type

		# Tier-5: Tier-4 and Tier-5 ||  Tier-5 and Tier-4
		elif 'T4_T5' in val_string or \
		     'T5_T4' in val_string:
			if val_string == 'T4_T5_T' or \
				val_string == 'T5_T4_T':
				tag_info_set.add('T5')
			else:
				val_type = val_string.split('_')[-1]
				if val_type == 'Align':
					tag_info_set.add('T5_Align')
				else:
					tag_info_set.add('D5_{0}'.format(val_string.split('_')[-1]))
				del val_type

		# Tier-5: Tier-5 and Tier-5
		elif 'T5_T5' in val_string:
			if val_string == 'T5_T5_T':
				tag_info_set.add('T5')
			else:
				val_type = val_string.split('_')[-1]
				if val_type == 'Align':
					tag_info_set.add('T5_Align')
				else:
					tag_info_set.add('D5_{0}'.format(val_string.split('_')[-1]))
				del val_type

		# False case
		else:
			tag_info_set.add('F')

	return tag_info_set

def pattern_checker_final(tag_info_set):
	# Tier-1
	if len([a for a in tag_info_set if 'T1' in a]) > 0:
		return [a for a in tag_info_set if 'T1' in a]

	# Tier-2
	elif len([a for a in tag_info_set if 'T1' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T2' in a]) > 0:
		return [a for a in tag_info_set if 'T2' in a]

	# Tier-3
	elif len([a for a in tag_info_set if 'T1' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T2' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T3' in a]) > 0:
		return [a for a in tag_info_set if 'T3' in a]

	# Tier-4
	elif len([a for a in tag_info_set if 'T1' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T2' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T3' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T4' in a]) > 0:
		return [a for a in tag_info_set if 'T4' in a]

	# Tier-5
	elif len([a for a in tag_info_set if 'T1' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T2' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T3' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T4' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T5' in a]) > 0:
		return [a for a in tag_info_set if 'T5' in a]

	# D1
	elif len([a for a in tag_info_set if 'T1' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T2' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T3' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T3' in a]) == 0 and \
		 len([a for a in tag_info_set if 'D1' in a]) > 0:
		return [a for a in tag_info_set if 'D1' in a]

	# D2
	elif len([a for a in tag_info_set if 'T1' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T2' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T3' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T3' in a]) == 0 and \
		 len([a for a in tag_info_set if 'D1' in a]) == 0 and \
		 len([a for a in tag_info_set if 'D2' in a]) > 0:
		return [a for a in tag_info_set if 'D2' in a]

	# D3
	elif len([a for a in tag_info_set if 'T1' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T2' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T3' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T3' in a]) == 0 and \
		 len([a for a in tag_info_set if 'D1' in a]) == 0 and \
		 len([a for a in tag_info_set if 'D2' in a]) == 0 and \
		 len([a for a in tag_info_set if 'D3' in a]) > 0:
		return [a for a in tag_info_set if 'D3' in a]

	# D4
	elif len([a for a in tag_info_set if 'T1' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T2' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T3' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T3' in a]) == 0 and \
		 len([a for a in tag_info_set if 'D1' in a]) == 0 and \
		 len([a for a in tag_info_set if 'D2' in a]) == 0 and \
		 len([a for a in tag_info_set if 'D3' in a]) == 0 and \
		 len([a for a in tag_info_set if 'D4' in a]) > 0:
		return [a for a in tag_info_set if 'D4' in a]

	# D5
	elif len([a for a in tag_info_set if 'T1' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T2' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T3' in a]) == 0 and \
		 len([a for a in tag_info_set if 'T3' in a]) == 0 and \
		 len([a for a in tag_info_set if 'D1' in a]) == 0 and \
		 len([a for a in tag_info_set if 'D2' in a]) == 0 and \
		 len([a for a in tag_info_set if 'D3' in a]) == 0 and \
		 len([a for a in tag_info_set if 'D4' in a]) == 0 and \
		 len([a for a in tag_info_set if 'D5' in a]) > 0:
		return [a for a in tag_info_set if 'D5' in a]

	# False
	else:
		return ['F']

def SV_judgement_step(val_SV_Dic, judge_SVs_Dic):
	for sv_line, val_info_Dic in val_SV_Dic.items():
		t_tag_info_set = set(); n_tag_info_set = set()

		# Check patterns in tumor bam
		if len(val_info_Dic['T']) > 0:
			t_tag_info_set = pattern_checker(val_info_Dic['T'], t_tag_info_set)

			if len(t_tag_info_set) > 0:
				t_final_list = pattern_checker_final(t_tag_info_set)

				# Check patterns in normal bam
				if len(val_info_Dic['N']) > 0:
					n_tag_info_set = pattern_checker(val_info_Dic['N'], n_tag_info_set)

					if len(n_tag_info_set) == 0:
						if not sv_line in judge_SVs_Dic:
							judge_SVs_Dic.setdefault(sv_line, t_final_list)
				
						return judge_SVs_Dic
					else:
						n_final_list = pattern_checker_final(n_tag_info_set)
						n_tier_list = [a for a in n_final_list if re.match('T\d+$', a)]

						# Normal: not Tier-1~4
						if len(n_tier_list) == 0:
							if not sv_line in judge_SVs_Dic:
								judge_SVs_Dic.setdefault(sv_line, t_final_list)
							
							return judge_SVs_Dic

						else:
							t_final_list.append('G')

							if not sv_line in judge_SVs_Dic:
								judge_SVs_Dic.setdefault(sv_line, t_final_list)
					
							return judge_SVs_Dic
				else:
					if not sv_line in judge_SVs_Dic:
						judge_SVs_Dic.setdefault(sv_line, t_final_list)
						
					return judge_SVs_Dic
			else:
				if not sv_line in judge_SVs_Dic:
					judge_SVs_Dic.setdefault(sv_line, ['F'])
					
				return judge_SVs_Dic
		else:
			if not sv_line in judge_SVs_Dic:
				judge_SVs_Dic.setdefault(sv_line, ['F'])
				
			return judge_SVs_Dic

def readInfo_in_LongBams(svType, target_range, svInfo_list, long_bam_t, long_bam_n, \
						 reference_seq, frag_len_t, frag_len_n, o):

	# Load the BAM file
	try:
		if not reference_seq:
			bam_file_t = pysam.AlignmentFile(long_bam_t, "rb")
			bam_file_n = pysam.AlignmentFile(long_bam_n, "rb")
		else:
			bam_file_t = pysam.AlignmentFile(long_bam_t, "rc", reference_filename = reference_seq)
			bam_file_n = pysam.AlignmentFile(long_bam_n, "rc", reference_filename = reference_seq)
	except Exception as err:
		print("Error:", err, 
				"This can be caused by issues with the alignment file. "+
				"Please make sure that it is sorted and indexed before trying again",
				file = sys.stderr)
		sys.exit(0)

	# Take the features in long-read bams
	judge_SVs_Dic = dict()

	for sv_line in svInfo_list:
		sv_list = sv_line.strip().split('\t')

		# Set up zoom size about SV types
		## SV type: DEL, DUP, and INV
		if svType != 'TRA':

			# Set up breakpoint information of target SV
			sv_c1 = ''; sv_p1 = 0; sv_c2 = ''; sv_p2 = 0

			if int(sv_list[1]) < int(sv_list[3]):
				sv_c1 = sv_list[0]
				sv_p1 = int(sv_list[1])
				sv_c2 = sv_list[2]
				sv_p2 = int(sv_list[3])
			else:
				sv_c1 = sv_list[2]
				sv_p1 = int(sv_list[3])
				sv_c2 = sv_list[0]
				sv_p2 = int(sv_list[1])

			windows_t = int((sv_p2 - sv_p1) / 2)
			if (0 < frag_len_t) and (windows_t < 1.5 * frag_len_t):
				old_windows_t = windows_t
				windows_t = int(1.5 * frag_len_t)
				print("Window size is under 1.5x the estimated fragment length " \
						+ "and will be resized to {}. Rerun with -w {} to override".format(windows_t, old_windows_t), file=sys.stderr,)
			# if region is larger than zoom, set window to zoom and set two ranges
			if windows_t >= 500000:
				windows_t = 500000
				ranges_t = [genome_interval(sv_c1, max(0, sv_p1 - windows_t), sv_p1 + windows_t), \
						    genome_interval(sv_c2, max(0, sv_p2 - windows_t), sv_p2 + windows_t),]
			else:
				ranges_t = [genome_interval(sv_c1, max(0, sv_p1 - windows_t), sv_p2 + windows_t)]

			windows_n = int((sv_p2 - sv_p1) / 2)
			if (0 < frag_len_n) and (windows_n < 1.5 * frag_len_n):
				old_windows_n = windows_n
				windows_n = int(1.5 * frag_len_n)
				print("Window size is under 1.5x the estimated fragment length " \
						+ "and will be resized to {}. Rerun with -w {} to override".format(windows_n, old_windows_n), file=sys.stderr,)
			# if region is larger than zoom, set window to zoom and set two ranges
			if windows_n >= 500000:
				windows_n = 500000
				ranges_n = [genome_interval(sv_c1, max(0, sv_p1 - windows_n), sv_p1 + windows_n), \
						    genome_interval(sv_c2, max(0, sv_p2 - windows_n), sv_p2 + windows_n),]
			else:
				ranges_n = [genome_interval(sv_c1, max(0, sv_p1 - windows_n), sv_p2 + windows_n)]

		## SV type: TRA
		else:
			windows_t = 1000
			ranges_t = [genome_interval(sv_list[0], max(0, int(sv_list[1]) - windows_t), int(sv_list[1]) + windows_t), \
					    genome_interval(sv_list[2], max(0, int(sv_list[3]) - windows_t), int(sv_list[3]) + windows_t),]

			windows_n = 1000
			ranges_n = [genome_interval(sv_list[0], max(0, int(sv_list[1]) - windows_n), int(sv_list[1]) + windows_n), \
					    genome_interval(sv_list[2], max(0, int(sv_list[3]) - windows_n), int(sv_list[3]) + windows_n),]

		###########################################
		# Check the input bam file: [tumor, normal]
		val_SV_Dic = dict()

		q_1 = queue.Queue()
		t_1 = threading.Thread(target=SV_judgement_runner, args=(bam_file_t, ranges_t, svType, target_range, sv_list, q_1,))
		t_1.start()

		q_2 = queue.Queue()
		t_2 = threading.Thread(target=SV_judgement_runner, args=(bam_file_n, ranges_n, svType, target_range, sv_list, q_2,))
		t_2.start()

		## Tumor long-read bam
		if not sv_line in val_SV_Dic:
			val_SV_Dic.setdefault(sv_line, dict()) \
						.setdefault('T', q_1.get())
		else:
			val_SV_Dic[sv_line].setdefault('T', q_1.get())

		## Tumor long-read bam
		if not sv_line in val_SV_Dic:
			val_SV_Dic.setdefault(sv_line, dict()) \
						.setdefault('N', q_2.get())
		else:
			val_SV_Dic[sv_line].setdefault('N', q_2.get())

		t_1.join(); t_2.join()

		del t_1; del q_1; del t_2; del q_2

		# Check SVs by short-read on long-read Bam file
		judge_SVs_Dic = SV_judgement_step(val_SV_Dic, judge_SVs_Dic)

	o.put(judge_SVs_Dic)

def get_read_info_in_LongBams(varDic, cpus, long_bam_t, long_bam_n, \
							  reference_seq, outname, target_type, \
							  target_range, frag_len_t, frag_len_n):

	t_judge_SVs_Dic = dict()

	for svInfo_set in varDic.values():

		# Using: 1 cpu
		if len(svInfo_set) < cpus:
			p_list = list(); o_list = list()
			tmp_list = [a for a in svInfo_set]

			for x in range(0, len(tmp_list)):
				o = Queue(); o_list.append(o)
				p = Process(target = readInfo_in_LongBams, \
							args = (target_type, target_range, [tmp_list[x]], long_bam_t, long_bam_n, reference_seq, 
									frag_len_t, frag_len_n, o))
				p_list.append(p)
			
			for e_p in p_list: e_p.start()
			for e_o in o_list:
				r = e_o.get()

				# Save format: {SV line: [SV validantion results]}
				for sv_line, sv_judge_list in r.items():
					t_judge_SVs_Dic.setdefault(sv_line, sv_judge_list)

				e_o.close()
				del r
			for e_p in p_list: e_p.join()

			del p_list; del o_list; del tmp_list

		# Using: multiprocessing
		else:
			p_list = list(); o_list = list()
			tmp_list = [a for a in svInfo_set]
			split_flag = math.floor(len(tmp_list)/float(cpus)); step_flag = 0

			for x in range(0, cpus):
				if x == 0:
					o = Queue(); o_list.append(o)
					p = Process(target = readInfo_in_LongBams, \
								args = (target_type, target_range, tmp_list[0:split_flag], long_bam_t, long_bam_n, reference_seq, \
										frag_len_t, frag_len_n, o))
					p_list.append(p); step_flag += split_flag
				elif x < cpus-1:
					o = Queue(); o_list.append(o)
					p = Process(target = readInfo_in_LongBams, \
								args = (target_type, target_range, tmp_list[step_flag:step_flag+split_flag], \
										long_bam_t, long_bam_n, reference_seq, frag_len_t, frag_len_n, o))
					p_list.append(p); step_flag += split_flag
				else:
					o = Queue(); o_list.append(o)
					p = Process(target = readInfo_in_LongBams, \
								args = (target_type, target_range, tmp_list[step_flag:], long_bam_t, long_bam_n, \
										reference_seq, frag_len_t, frag_len_n, o))
					p_list.append(p)

			for e_p in p_list: e_p.start()
			for e_o in o_list:
				r = e_o.get()

				# Save format: {SV line: [SV validantion results]}
				for sv_line, sv_judge_list in r.items():
					t_judge_SVs_Dic.setdefault(sv_line, sv_judge_list)
				
				e_o.close()
				del r
			for e_p in p_list: e_p.join()

			del p_list; del o_list; del tmp_list; del split_flag; del step_flag

	# Output the validation results to outfile
	outfile = open('{0}.{1}.info.list'.format(outname, target_type), 'w')
	outfile.write('SV_type\tChr_A\tPos_A\tChr_B\tPos_B\tJudgement_by_longRead\n')

	for sv_line, sv_judge_list in t_judge_SVs_Dic.items():
		sv_judge_list.sort()
		sv_judge_line = ';'.join(sv_judge_list)

		outfile.write(('{0}\t{1}\t{2}\n').format(target_type, sv_line, sv_judge_line))

		del sv_judge_line

	outfile.close()

def SV_judgement_step_without_normal(val_SV_Dic, judge_SVs_Dic):
	for sv_line, val_info_Dic in val_SV_Dic.items():
		t_tag_info_set = set(); n_tag_info_set = set()

		# Check patterns in tumor bam
		if len(val_info_Dic['T']) > 0:
			t_tag_info_set = pattern_checker(val_info_Dic['T'], t_tag_info_set)

			if len(t_tag_info_set) > 0:
				t_final_list = pattern_checker_final(t_tag_info_set)

				if not sv_line in judge_SVs_Dic:
					judge_SVs_Dic.setdefault(sv_line, t_final_list)
					
				return judge_SVs_Dic

			else:
				if not sv_line in judge_SVs_Dic:
					judge_SVs_Dic.setdefault(sv_line, ['F'])
					
				return judge_SVs_Dic
		else:
			if not sv_line in judge_SVs_Dic:
				judge_SVs_Dic.setdefault(sv_line, ['F'])
				
			return judge_SVs_Dic

def readInfo_in_LongBams_without_normal(svType, target_range, svInfo_list, long_bam_t, \
						 				reference_seq, frag_len_t, o):
	# Load the BAM file
	try:
		if not reference_seq:
			bam_file_t = pysam.AlignmentFile(long_bam_t, "rb")
		else:
			bam_file_t = pysam.AlignmentFile(long_bam_t, "rc", reference_filename = reference_seq)
	except Exception as err:
		print("Error:", err, 
				"This can be caused by issues with the alignment file. "+
				"Please make sure that it is sorted and indexed before trying again",
				file = sys.stderr)
		sys.exit(0)

	# Take the features in long-read bams
	judge_SVs_Dic = dict()

	for sv_line in svInfo_list:
		sv_list = sv_line.strip().split('\t')

		# Set up zoom size about SV types
		## SV type: DEL, DUP, and INV
		if svType != 'TRA':

			# Set up breakpoint information of target SV
			sv_c1 = ''; sv_p1 = 0; sv_c2 = ''; sv_p2 = 0

			if int(sv_list[1]) < int(sv_list[3]):
				sv_c1 = sv_list[0]
				sv_p1 = int(sv_list[1])
				sv_c2 = sv_list[2]
				sv_p2 = int(sv_list[3])
			else:
				sv_c1 = sv_list[2]
				sv_p1 = int(sv_list[3])
				sv_c2 = sv_list[0]
				sv_p2 = int(sv_list[1])

			windows_t = int((sv_p2 - sv_p1) / 2)
			if (0 < frag_len_t) and (windows_t < 1.5 * frag_len_t):
				old_windows_t = windows_t
				windows_t = int(1.5 * frag_len_t)
				print("Window size is under 1.5x the estimated fragment length " \
						+ "and will be resized to {}. Rerun with -w {} to override".format(windows_t, old_windows_t), file=sys.stderr,)
			# if region is larger than zoom, set window to zoom and set two ranges
			if windows_t >= 500000:
				windows_t = 500000
				ranges_t = [genome_interval(sv_c1, max(0, sv_p1 - windows_t), sv_p1 + windows_t), \
						    genome_interval(sv_c2, max(0, sv_p2 - windows_t), sv_p2 + windows_t),]
			else:
				ranges_t = [genome_interval(sv_c1, max(0, sv_p1 - windows_t), sv_p2 + windows_t)]

		## SV type: TRA
		else:
			windows_t = 1000
			ranges_t = [genome_interval(sv_list[0], max(0, int(sv_list[1]) - windows_t), int(sv_list[1]) + windows_t), \
					    genome_interval(sv_list[2], max(0, int(sv_list[3]) - windows_t), int(sv_list[3]) + windows_t),]

		###########################################
		# Check the input bam file: [tumor, normal]
		val_SV_Dic = dict()

		q_1 = queue.Queue()
		t_1 = threading.Thread(target=SV_judgement_runner, args=(bam_file_t, ranges_t, svType, target_range, sv_list, q_1,))
		t_1.start()

		## Tumor long-read bam
		if not sv_line in val_SV_Dic:
			val_SV_Dic.setdefault(sv_line, dict()) \
						.setdefault('T', q_1.get())
		else:
			val_SV_Dic[sv_line].setdefault('T', q_1.get())

		t_1.join()

		del t_1; del q_1

		# Check SVs by short-read on long-read Bam file
		judge_SVs_Dic = SV_judgement_step_without_normal(val_SV_Dic, judge_SVs_Dic)

	o.put(judge_SVs_Dic)

def get_read_info_in_LongBams_without_normal(varDic, cpus, long_bam_t, reference_seq, \
											 outname, target_type, \
							  				 target_range, frag_len_t):

	t_judge_SVs_Dic = dict()

	for svInfo_set in varDic.values():

		# Using: 1 cpu
		if len(svInfo_set) < cpus:
			p_list = list(); o_list = list()
			tmp_list = [a for a in svInfo_set]

			for x in range(0, len(tmp_list)):
				o = Queue(); o_list.append(o)
				p = Process(target = readInfo_in_LongBams_without_normal, \
							args = (target_type, target_range, [tmp_list[x]], long_bam_t, reference_seq, 
									frag_len_t, o))
				p_list.append(p)
			
			for e_p in p_list: e_p.start()
			for e_o in o_list:
				r = e_o.get()

				# Save format: {SV line: [SV validantion results]}
				for sv_line, sv_judge_list in r.items():
					t_judge_SVs_Dic.setdefault(sv_line, sv_judge_list)

				e_o.close()
				del r
			for e_p in p_list: e_p.join()

			del p_list; del o_list; del tmp_list

		# Using: multiprocessing
		else:
			p_list = list(); o_list = list()
			tmp_list = [a for a in svInfo_set]
			split_flag = math.floor(len(tmp_list)/float(cpus)); step_flag = 0

			for x in range(0, cpus):
				if x == 0:
					o = Queue(); o_list.append(o)
					p = Process(target = readInfo_in_LongBams_without_normal, \
								args = (target_type, target_range, tmp_list[0:split_flag], long_bam_t, reference_seq, \
										frag_len_t, o))
					p_list.append(p); step_flag += split_flag
				elif x < cpus-1:
					o = Queue(); o_list.append(o)
					p = Process(target = readInfo_in_LongBams_without_normal, \
								args = (target_type, target_range, tmp_list[step_flag:step_flag+split_flag], \
										long_bam_t, reference_seq, frag_len_t, o))
					p_list.append(p); step_flag += split_flag
				else:
					o = Queue(); o_list.append(o)
					p = Process(target = readInfo_in_LongBams_without_normal, \
								args = (target_type, target_range, tmp_list[step_flag:], long_bam_t, \
										reference_seq, frag_len_t, o))
					p_list.append(p)

			for e_p in p_list: e_p.start()
			for e_o in o_list:
				r = e_o.get()

				# Save format: {SV line: [SV validantion results]}
				for sv_line, sv_judge_list in r.items():
					t_judge_SVs_Dic.setdefault(sv_line, sv_judge_list)
				
				e_o.close()
				del r
			for e_p in p_list: e_p.join()

			del p_list; del o_list; del tmp_list; del split_flag; del step_flag

	# Output the validation results to outfile
	outfile = open('{0}.{1}.info.list'.format(outname, target_type), 'w')
	outfile.write('SV_type\tChr_A\tPos_A\tChr_B\tPos_B\tJudgement_by_longRead\n')

	for sv_line, sv_judge_list in t_judge_SVs_Dic.items():
		sv_judge_list.sort()
		sv_judge_line = ';'.join(sv_judge_list)

		outfile.write(('{0}\t{1}\t{2}\n').format(target_type, sv_line, sv_judge_line))

		del sv_judge_line

	outfile.close()

def main(args):
	# Save the SVs to dictionary
	varDic = vcf_to_Dic(args.input_VCF, args.target_type)

	# Estimate fragment size
	frag_len_t = estimate_fragment_len(args.long_bam_t, args.reference_seq)

	if len(args.long_bam_n) > 0:
		frag_len_n = estimate_fragment_len(args.long_bam_n, args.reference_seq)

		# Validate the SVs on short-reads using long-read bam files
		get_read_info_in_LongBams(varDic, args.cpus, args.long_bam_t, args.long_bam_n, \
								  args.reference_seq, args.outname, args.target_type, \
								  args.target_range, frag_len_t, frag_len_n)
	else:
		frag_len_n = 0
		# Validate the SVs on short-reads using long-read bam files
		get_read_info_in_LongBams_without_normal(varDic, args.cpus, args.long_bam_t, \
								  				 args.reference_seq, args.outname, args.target_type, \
								  				 args.target_range, frag_len_t)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'LR_Filter is SV filtering tool using long-read BAM file')
	parser.add_argument('-i', '--input_VCF', type = str, help = 'The input VCF file by short-read SV caller, etc ETCHING, Lumpy, Delly ... ')
	parser.add_argument('-t', '--target_type', type = str, help = 'The target SV type, etc DEL, DUP, INV, or TRA')
	parser.add_argument('-lbt', '--long_bam_t', type = str, default = '', help = 'A tumor long-read mapping sorted BAM file by Minimap2 or NGMLR')
	parser.add_argument('-lbn', '--long_bam_n', type = str, default = '', help = ''A normal long-read mapping sorted BAM file by Minimap2 or NGMLR'')
	parser.add_argument('-rf', '--reference_seq', type = str, help = 'The reference fasta file, etc hg19 or hg38')
	parser.add_argument('-c', '--cpus', type = int, help = 'Set the number of CPUs, this option is for multi-processing')
	parser.add_argument('-tr', '--target_range', type = int, default = 500, help = 'Set the range of SV BP for verifying by long-reads')
	parser.add_argument('-o', '--outname', type = str, help = 'Output name')
	args = parser.parse_args()
	main(args)

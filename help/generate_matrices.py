# Input: 
#       -d: sample.vcf in a patient folder
# Author: Jingyi Wang
# Created date: 2017_09_25
# Modified date: 2017_10_02

##################
##### IMPORT #####
##################

import sys      
import os       
import argparse 
import vcf      
import numpy as np
import operator
import random

# custom imports
import file_manager as fm

#####################
##### FUNCTIONS #####
#####################

#  input: in_dir (str) full path to input directory containing .vcf file(s)
def get_mats(in_dir):
	sampleList = fm._fnames_with_extension(in_dir, '.vcf')

	m = len(sampleList)

	bp_id_to_mate_id, bp_id_to_tuple = {}, {}

	BP_sample_dict, CN_sample_dict, CN_sample_rec_dict = dict(), dict(), dict()
	for i, sample in enumerate(sampleList):
		input_vcf_file = in_dir + '/' + sample
		reader = vcf.Reader(open(input_vcf_file, 'r'))
		BP_sample_dict[sample], CN_sample_dict[sample], CN_sample_rec_dict[sample], mateIDs, toTuple = get_sample_dict(reader)
		
		# prepend sample index to each breakpoint ID
		for k, v in mateIDs.iteritems():
			bp_id_to_mate_id[str(i+1) + k] = str(i+1) + v # add all entries from mateIDs (dict) to bp_id_to_mate_id (dict)
			bp_id_to_tuple[str(i+1) + k] = toTuple[k]     # add all entries from toTuple (dict) to bp_id_to_tuple (dict)

	BP_idx_dict, l = get_BP_idx_dict(BP_sample_dict)
	CN_startPos_dict, CN_endPos_dict, r = get_CN_indices_dict(CN_sample_dict)

	F, Q, A, H = make_matrices(m, l, r, sampleList, BP_sample_dict, BP_idx_dict, CN_sample_rec_dict, CN_startPos_dict, CN_endPos_dict)
	
	F = np.array(F).astype(float)
	Q = np.array(Q)
	A = np.array(A)
	H = np.array(H)

	G = make_G(BP_idx_dict, bp_id_to_mate_id, bp_id_to_tuple)

	return F, Q, G, A, H


# init a r*c 2d list filling with 0
def make_2d_list(r,c):
    result = list()
    for i in range(r):
        result.append([0] * c)
    return result


# input dictionaries, output four matrices
def make_matrices(m, l, r, sampleList, BP_sample_dict, BP_idx_dict, CN_sample_rec_dict, CN_startPos_dict, CN_endPos_dict):
	F, Q, A, H = make_2d_list(m, l + r), make_2d_list(l, r), make_2d_list(m, l), make_2d_list(m, l)

	for sample_idx in range(len(sampleList)):
		sample = sampleList[sample_idx]
		for chrom in BP_sample_dict[sample]:
			for pos in BP_sample_dict[sample][chrom]:
				temp_bp_info_dict = BP_sample_dict[sample][chrom][pos] # dictionary
				cn, direction, bdp, dp = temp_bp_info_dict['cn'], temp_bp_info_dict['dir'], temp_bp_info_dict['bdp'], temp_bp_info_dict['dp']
				bp_idx = BP_idx_dict[(chrom, pos, direction)]
				F[sample_idx][bp_idx] = cn
				A[sample_idx][bp_idx] = bdp
				H[sample_idx][bp_idx] = dp

				if direction == False and (chrom, pos) in CN_endPos_dict:
					cn_idx = CN_endPos_dict[(chrom, pos)]
					Q[bp_idx][cn_idx] = 1
				elif direction == True and (chrom, pos) in CN_startPos_dict:
					cn_idx = CN_startPos_dict[(chrom, pos)]
					Q[bp_idx][cn_idx] = 1

		for chrom in CN_sample_rec_dict[sample]:
			for (s,e) in CN_sample_rec_dict[sample][chrom]:
				cn_idx_list = get_CN_indices(CN_startPos_dict, CN_endPos_dict, chrom, s, e)
				cn = CN_sample_rec_dict[sample][chrom][(s,e)]
				for cn_idx in cn_idx_list:
					F[sample_idx][cn_idx + l] = cn

	return F, Q, A, H


# key: (chrom, pos, dir)
# val: idx (idx starts from 0)
def get_BP_idx_dict(BP_sample_dict):
	chrom_pos_dir_dict = dict() # key: chrom, val: set of (pos, dir) tuples
	for sample in BP_sample_dict:
		for chrom in BP_sample_dict[sample]:
			if chrom not in chrom_pos_dir_dict:
				chrom_pos_dir_dict[chrom] = set()
			for pos in BP_sample_dict[sample][chrom]:
				direction = BP_sample_dict[sample][chrom][pos]['dir']
				chrom_pos_dir_dict[chrom].add((pos, direction))

	chrom_pos_dict = dict() # key: chrom, val: sorted pos list (ascending)
	for chrom in chrom_pos_dir_dict:
		sorted_pos_list = sorted(map(operator.itemgetter(0), chrom_pos_dir_dict[chrom]))
		chrom_pos_dict[chrom] = sorted_pos_list

	BP_patient_dict = dict()
	for chrom in chrom_pos_dict:
		BP_patient_dict[chrom] = list()
		for pos in chrom_pos_dict[chrom]:
			if (pos, False) in chrom_pos_dir_dict[chrom]:
				BP_patient_dict[chrom].append((pos, False))
			if (pos, True) in chrom_pos_dir_dict[chrom]:
				BP_patient_dict[chrom].append((pos, True))

	BP_idx_dict = dict()
	idx = 0
	sorted_chrom = sorted(BP_patient_dict.keys())
	for chrom in sorted_chrom:
		for (pos, direction) in BP_patient_dict[chrom]:
			BP_idx_dict[(chrom, pos, direction)] = idx
			idx += 1
	l = idx
	return BP_idx_dict, l


# return three dictionaries: BP_sample_dict, CN_sample_dict, and CN_sample_rec_dict
# 1. BP_sample_dict: 
#    key: sample
#    val: {chr1:{ pos1:{dir: T/F, cn: , a: , h: }, pos2: {dir: , cn: , a:, h: }, chr2: {}, ...}
#    mate_dir = rec.ALT[0].remoteOrientation. if mate_dir == False: ], if mate_dir == True: [
# 2. CN_sample_dict:
#    key: sample
#    value: {chr1: {pos1: s/e, pos2: s/e, ...}, chr2: {}, ... }
# 3. CN_sample_rec_dict: 
#    key: sample
#    value: {chr1: {(s1, e1): cn1, (s2, e2): cn2, ...}, chr2: {...}...}
# 4. bp_id_to_mate_id (dict) key (str) is ID of breakpoint. val (str) is ID of mate
# 5. bp_id_to_tuple   (dict) key (str) is ID of breakpoint. val (tuple) is (chrm_num, pos, direction)
def get_sample_dict(reader):
	BP_sample_dict, CN_sample_dict, CN_sample_rec_dict = dict(), dict(), dict()
	bp_id_to_mate_id = {} # key is id (str). val is mate id (str)
	bp_id_to_tuple = {}   # key is (chrm_num, pos, direction). key is id (str)
	for rec in reader:
		if is_sv_record(rec) == True:
			if rec.CHROM not in BP_sample_dict:
				BP_sample_dict[rec.CHROM] = dict()
			if rec.POS not in BP_sample_dict[rec.CHROM]:
				BP_sample_dict[rec.CHROM][rec.POS] = dict()
			BP_sample_dict[rec.CHROM][rec.POS]['id'] = rec.ID
			BP_sample_dict[rec.CHROM][rec.POS]['cn'] = rec.samples[0].data.CNADJ
			BP_sample_dict[rec.CHROM][rec.POS]['mate_dir'] = rec.ALT[0].remoteOrientation
			# BP_sample_dict[rec.CHROM][rec.POS]['mate_id'] = rec.INFO['MATEID']
			BP_sample_dict[rec.CHROM][rec.POS]['mate_pos'] = rec.ALT[0].pos

			bp_id_to_mate_id[rec.ID] = rec.INFO['MATEID'][0]
			
			BP_sample_dict[rec.CHROM][rec.POS]['bdp'] = rec.samples[0].data.BDP
			BP_sample_dict[rec.CHROM][rec.POS]['dp'] = rec.samples[0].data.DP
		elif is_cnv_record(rec) == True:
			if rec.CHROM not in CN_sample_dict:
				CN_sample_dict[rec.CHROM] = dict()
				CN_sample_rec_dict[rec.CHROM] = dict()
			CN_sample_dict[rec.CHROM][rec.POS] = ['s']
			CN_sample_dict[rec.CHROM][rec.INFO['END']] = ['e']
			CN_sample_rec_dict[rec.CHROM][(rec.POS, rec.INFO['END'])] = sum(rec.samples[0].data.CN)

	for chrom in BP_sample_dict:
		for pos in BP_sample_dict[chrom]:
			BP_sample_dict[chrom][pos]['dir'] = BP_sample_dict[chrom][BP_sample_dict[chrom][pos]['mate_pos']]['mate_dir']

			bp_id = BP_sample_dict[chrom][pos]['id']
			bp_id_to_tuple[bp_id] = (chrom, pos, BP_sample_dict[chrom][pos]['dir'])

	return BP_sample_dict, CN_sample_dict, CN_sample_rec_dict, bp_id_to_mate_id, bp_id_to_tuple


# CN_startPos_dict: key: (chrom, startPos), val: idx
# CN_endPos_dict: key: (chrom, endPos), val: idx
def get_CN_indices_dict(CN_sample_dict):

	chrom_dict = dict() # key: chrom, value: list of samples contain this chrom
	for sample in CN_sample_dict:
		for chrom in CN_sample_dict[sample]:
			if chrom not in chrom_dict:
				chrom_dict[chrom] = list()
			chrom_dict[chrom].append(sample)

	CN_patient_dict = dict()
	for chrom in chrom_dict:
		posSet = set()
		pos_dir_dict = dict() # given pos, output ['e'] or ['s'] or ['s', 'e']
		for sample in chrom_dict[chrom]:
			for pos in CN_sample_dict[sample][chrom]:
				posSet.add(pos)
				if pos not in pos_dir_dict:
					pos_dir_dict[pos] = CN_sample_dict[sample][chrom][pos]
				else:
					if CN_sample_dict[sample][chrom][pos] != pos_dir_dict[pos]:
						pos_dir_dict[pos] += CN_sample_dict[sample][chrom][pos]
		posList = sorted(list(posSet))

		CN_patient_dict[chrom] = list()
		tempS = posList[0]
		idx = 0
		while idx < len(posList) - 1:
			if 's' in pos_dir_dict[posList[idx + 1]] and 'e' not in pos_dir_dict[posList[idx + 1]]:
				tempE = posList[idx + 1] - 1
				CN_patient_dict[chrom].append((tempS, tempE))
				tempS = posList[idx + 1]
				idx += 1
			elif 'e' in pos_dir_dict[posList[idx + 1]] and 's' not in pos_dir_dict[posList[idx + 1]]:
				tempE = tempE = posList[idx + 1]
				CN_patient_dict[chrom].append((tempS, tempE))
				if idx + 1 < len(posList) - 1:
					tempS = posList[idx + 2]
				idx += 2
			elif 's' in pos_dir_dict[posList[idx + 1]] and 'e' in pos_dir_dict[posList[idx + 1]]:
				tempE = posList[idx + 1] - 1
				CN_patient_dict[chrom].append((tempS, tempE))
				CN_patient_dict[chrom].append((posList[idx + 1], posList[idx + 1]))
				tempS = posList[idx + 1] + 1
				idx += 1

	CN_startPos_dict = dict()
	CN_endPos_dict = dict()
	idx = 0
	sorted_chrom = sorted(CN_patient_dict.keys())
	for chrom in sorted_chrom:
		for (s,e) in CN_patient_dict[chrom]:
			CN_startPos_dict[(chrom, s)] = idx
			CN_endPos_dict[(chrom, e)] = idx
			idx += 1
	r = idx

	return CN_startPos_dict, CN_endPos_dict, r


#  input: bp_tuple_to_idx (dict) key is bp tuple (chrm, pos, direction). val is index of output G
#         bp_id_to_mate_id (dict) key (str) is ID of breakpoint. val (str) is ID of mate
#         bp_id_to_tuple   (dict) key (str) is ID of breakpoint. val (tuple) is (chrm_num, pos, direction)
# output: G (np.array of 0 or 1) [l, l] g_s,t == 1 if breakpoints s and t are mates. 0 otherwise
def make_G(bp_tuple_to_idx, bp_id_to_mate_id, bp_id_to_tuple):
	l = len(bp_tuple_to_idx.keys())
	G = np.zeros((l, l))

	for i in xrange(0, l):
		G[i, i] = 1          # breakpoint being its own mate is a requirement for the solver

	bp_idx_to_tuple = inv_dict(bp_tuple_to_idx)
	bp_tuple_to_mate_tuple = get_bp_tuple_to_mate_tuple(bp_id_to_mate_id, bp_id_to_tuple)

	for i in sorted(bp_idx_to_tuple.iterkeys()):
		cur_tuple = bp_idx_to_tuple[i]
		mate_tup = bp_tuple_to_mate_tuple[cur_tuple]
		j = bp_tuple_to_idx[mate_tup]
		G[i, j] = 1

	return G

#  input: bp_id_to_mate_id (dict) key (str) is ID of breakpoint. val (str) is ID of mate
#         bp_id_to_tuple   (dict) key (str) is ID of breakpoint. val (tuple) is (chrm_num, pos, direction)
# output: bp_tuple_to_mate_tuple (dict) key is (tuple) is (chrm_num, pos, direction). val is mate tuple
def get_bp_tuple_to_mate_tuple(bp_id_to_mate_id, bp_id_to_tuple):
	out_dic = {}
	bp_tuple_to_id = inv_dict(bp_id_to_tuple)
	for tup, cur_id in bp_tuple_to_id.iteritems():
		mate_id = bp_id_to_mate_id[cur_id]
		mate_tup = bp_id_to_tuple[mate_id]
		out_dic[tup] = mate_tup
	return out_dic

# inverts dictionary. keys become values. values become keys. must be able to be inverted so
#   input keys must be static
def inv_dict(dic):
	idic = {}
	for k, v in dic.iteritems():
		idic[v] = k
	return idic


# given start and end position, output list of segment indices (continuous)
def get_CN_indices(CN_startPos_dict, CN_endPos_dict, chrom, s, e):
	result = list()
	firstIdx = CN_startPos_dict[(chrom,s)]
	endIdx = CN_endPos_dict[(chrom, e)]
	for i in range(firstIdx, endIdx + 1):
		result.append(i)
	return result


def is_cnv_record(rec):
	return rec.ID[0:3] == 'cnv'


def is_sv_record(rec):
	return rec.ID[0:2] == 'sv'

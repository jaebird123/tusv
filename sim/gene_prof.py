#     file: gene_prof.py
#   author: Jesse Eaton
#  created: 9/27/2017
# modified: 9/27/2017
#  purpose: GeneProf class. GeneProf is a gene profile and contains a mutated genome with
#             references to the mutations' original position. Keeps track of copy numbers.

# imports
import sys
import copy

# helpers
def printnow(s):
	sys.stdout.write(s)
	sys.stdout.flush()
def printerr(s):
	sys.stderr.write(s)
	sys.stderr.flush()

class ChrmProf:
	# chrm (str) entire chromosome string. ex: "ATCCGA"
	def __init__(self, chrm):
		self.chrm = chrm
		self.org = _OrgNode(0, len(chrm) - 1)
		self.mut = _MutNode(0, len(chrm) - 1)
		self.org.children.append(self.mut)
		self.mut.parent = self.org

	# output: bgns (list of int) [n] beginning positions for each segment
	#         ends (list of int) [n] ending positions for each segment
	#         cps  (list of int) [n] copy number for each segment
	def get_copy_nums(self):
		cur = self.org
		bgns, ends, cps = [], [], []
		while cur != None:
			bgns.append(cur.bgn)
			ends.append(cur.end)
			cps.append(len(cur.children))
			cur = cur.r
		return bgns, ends, cps

	# duplicate region from bgn to end
	def dup(self, bgn, end):
		n = len(self.chrm)
		if bgn < 0 or bgn > end or end >= n:
			printerr('chromosome has length ' + str(n) + '. cannot dup [' + str(bgn) + ', ' + str(end) + ']')
			return
		if bgn > 0:
			self.split(bgn)     # do not split if bgn is 0 (start of chrm)
		if end + 1 < n:
			self.split(end + 1) # do not split if end is n-1 (end of chrm)
		
		insR, head, tail = _copy_from_to(self.mut, bgn, end) # copy list from bgn to end
		insL = insR.r # node to go after tail
		insR.r = head
		head.l = insR
		if insL != None:
			tail.r = insL
			insL.l = tail

		# increment bgn and end vales
		seg_len = end - bgn + 1
		while head != None:
			head.bgn += seg_len
			head.end += seg_len
			head = head.r

		# duplicate string
		s1 = self.chrm[:bgn]
		s2 = self.chrm[bgn:end+1]
		s3 = self.chrm[end+1:]
		self.chrm = s1 + s2 + s2 + s3

	# splits node at position k in mut and org. k is bgn of right offspring
	def split(self, k):

		# return if k out of bounds
		if k == 0 or k >= len(self.chrm):
			printerr('cannot split at k = ' + str(k) + '\n')
			return

		# find orgNode corresponding to the mutNode where split will occur
		splitMut = self.mut
		while splitMut != None and not (splitMut.bgn <= k and k <= splitMut.end):
			splitMut = splitMut.r
		orgNode1 = splitMut.parent

		k = k - splitMut.bgn # make k the new length of the old node (node that will be split)
		                     # is is now the number of nucletides in from a node where it should be split

		# split orgNode1 and all its children
		orgNode2 = orgNode1.split(k)
		for mutNode1 in orgNode1.children:
			mutNode2 = mutNode1.split(k)
			mutNode2.parent = orgNode2
			orgNode2.children.append(mutNode2)

	def pprint(self):
		printnow('sequence: ' + str(self.chrm) + '\n')
		printnow('lists:\n')
		for lst_name, cur in {'org': self.org, 'mut': self.mut}.iteritems():
			printnow(lst_name + ': ')
			while cur != None:
				cur.pprint()
				cur = cur.r
			printnow('\n')
		printnow('relations:\n')
		cur = self.org
		while cur != None:
			kid_pos_strs = [ kid.get_pos_str() for kid in cur.children ]
			printnow('[' + str(cur.bgn) + ',' + str(cur.end) + '] -> ' + ','.join(kid_pos_strs) + '\n')
			cur = cur.r

class _Node:
	def pprint(self):
		s = self.get_pos_str()
		if self.r != None:
			s += '->'
		else:
			s += '-v'
		printnow(s)
	def get_pos_str(self):
		return '[' + str(self.bgn) + ',' + str(self.end) + ']'

class _OrgNode(_Node):
	def __init__(self, bgn, end):
		self.children = [] # no mutated sections
		self.l = None      # no left or right pointers
		self.r = None
		self.bgn = bgn
		self.end = end

	# returns pointer to new sibling on right. k (int) means k + self.begin is bgn of new sibling
	def split(self, k):
		r = self.r
		self.r = _OrgNode(self.bgn + k, self.end)
		self.r.r = r    # set right of new node to the old node's old right
		self.r.l = self # set left of new node to old node (self)
		self.end = self.bgn + k - 1
		return self.r

class _MutNode(_Node):
	def __init__(self, bgn, end):
		self.parent = None
		self.l = None
		self.r = None
		self.bgn = bgn
		self.end = end

	# returns pointer to new sibling on right. k (int) means k + self.begin is bgn of new sibling
	def split(self, k):
		r = self.r
		self.r = _MutNode(self.bgn + k, self.end)
		self.r.r = r    # set right of new node to the old node's old right
		self.r.l = self # set left of new node to old node (self)
		self.end = self.bgn + k - 1
		return self.r

# helpers

# fm (int) is bgn index of one of the nodes. to (int) is end of one of the nodes
def _copy_from_to(head, fm, to):
	while head != None and head.bgn != fm: # make head the beggining of where to copy
		head = head.r

	curA = head
	i = 0
	prevB = None
	while curA != None:
		curB = copy.copy(curA)
		curB.parent.children.append(curB) # update parent's children pointers
		if i == 0:
			headB = curB
			i += 1
		curB.l = prevB
		if prevB != None:
			prevB.r = curB
		curB.r = None
		prevB = curB
		if curA.end == to:
			return curA, headB, curB
		curA = curA.r
	printerr('should not get here')

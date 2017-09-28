import gene_prof as gp

cp = gp.ChrmProf("AAABBBCCCTTT")
cp.pprint()

print '\n- - - - -\n'

cp.dup(3, 8)
cp.pprint()

print '\n- - - - -\n'

cp.dup(4, 15)
cp.pprint()

print '\n\n'
bgns, ends, cps = cp.get_copy_nums()
print bgns
print ends
print cps

print '\n- - - - -\n'

cp.dup(1, 28)
cp.pprint()

print '\n\n'
bgns, ends, cps = cp.get_copy_nums()
print bgns
print ends
print cps

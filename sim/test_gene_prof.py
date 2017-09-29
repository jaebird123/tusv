import gene_prof as gp

cp = gp.ChrmProf("AAABBBCCCDDDEEEFFF")
print ''

print cp.amp(4,10)
cp.pprint()
print ''

print cp.amp(11, 13)
cp.pprint()
print ''

print cp.amp(2, 19)
cp.pprint()
print ''

print cp.rem(0, 6)
cp.pprint()
print ''

print cp.rem(2, 8)
cp.pprint()
print ''





# cp.pprint()

# print '\n- - - - -\n'

# cp.amp(3, 8)
# cp.pprint()

# print '\n- - - - -\n'

# cp.amp(4, 15)
# cp.pprint()

# print '\n\n'
# bgns, ends, cps = cp.get_copy_nums()
# print bgns
# print ends
# print cps

# print '\n- - - - -\n'

# cp.amp(1, 28)
# cp.pprint()

# print '\n\n'
# bgns, ends, cps = cp.get_copy_nums()
# print bgns
# print ends
# print cps

import chrm_prof as chpr

cp = chpr.ChrmProf("AAABBBCCCDDDEEEFFF")
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

print cp.amp(9, 25)
cp.pprint()
print ''

cp2 = chpr.ChrmProf("RRRRRSSSSSTTTTTUUUUUVVVVV")
print ''

print cp2.amp(5, 14)
cp2.pprint()
print ''

print cp2.inv(3, 20)
cp2.pprint()
print ''

print cp2.inv(15, 30)
cp2.pprint()
print ''

print cp2.rem(10, 16)
cp2.pprint()
print ''

cp3 = chpr.ChrmProf("XXXXYYYYZZZZ")
print ''

print cp3.inv(0, 6)
cp3.pprint()
print ''

print cp3.inv(4, 11)
cp3.pprint()
print ''

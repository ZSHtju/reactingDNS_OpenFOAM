import numpy as np
from decimal import Decimal
import sys
import os

scheme = 'transport fitting'

with open("./tran.out") as o:
    content = o.readlines()
# break

for i in range(len(content)):
    content[i] = content[i].rstrip('\n')
    content[i] = content[i].rstrip('\r')
    if "COEFFICIENTS FOR SPECIES CONDUCTIVITIES" in content[i]:
        ini_lambda = i+3
        continue
    if "COEFFICIENTS FOR SPECIES VISCOSITIES" in content[i]:
        fin_lambda = i-3
        ini_mu = i+3
        continue
    if "COEFFICIENTS FOR SPECIES DIFFUSION COEFFICIENTS" in content[i]:
        fin_mu = i-3
        ini_bindiff = i+3
        continue
    if "COEFFICIENTS FOR THERMAL DIFFUSION RATIOS" in content[i]:
        fin_bindiff = i-3
        ini_thermdiff = i+3
        continue
species_lambda = content[ini_lambda:fin_lambda]
species_lambda = [name for name in species_lambda if name.strip()]

species_mu = content[ini_mu:fin_mu]
species_mu = [name for name in species_mu if name.strip()]

binary_diff = content[ini_bindiff:fin_bindiff]
binary_diff = [name for name in binary_diff if name.strip()]

therm_diff = content[ini_thermdiff:]
therm_diff = [name for name in therm_diff if name.strip()]

specl = [['']*5]*(len(species_lambda))
specm = [['']*5]*(len(species_mu))
bindf = [['']*6]*(len(binary_diff))
thermdf = [['']*6]*(len(therm_diff))

for j in range(len(species_lambda)):
    specl[j] = species_lambda[j].split()
for k in range(len(species_mu)):
    specm[k] = species_mu[k].split()
for ll in range(len(binary_diff)):
    bindf[ll] = binary_diff[ll].split()
for kk in range(len(therm_diff)):
    thermdf[kk] = therm_diff[kk].split()

with open("./original_form/speciesLambda") as f:
    with open("speciesLambda", "w") as fn1:
        for line in f:
            fn1.write(line)
            if "}" in line:
                break
fn1=open("speciesLambda", "a+")
orig_stdout = sys.stdout
sys.stdout = fn1
print('')
print(u'// Scheme: {0}'.format(scheme))
print('')

for p in range(len(species_lambda)):
    print(specl[p][0])
    print('{')
    print('Lambda1	{0}	;'.format(specl[p][1]))
    print('Lambda2	{0}	;'.format(specl[p][2]))
    print('Lambda3	{0}	;'.format(specl[p][3]))
    print('Lambda4	{0}	;'.format(specl[p][4]))
    print('}')
sys.stdout = orig_stdout
fn1.close()

with open("./original_form/speciesMu") as f:
    with open("speciesMu", "w") as fn1:
        for line in f:
            fn1.write(line)
            if "}" in line:
                break
fn1=open("speciesMu", "a+")
orig_stdout = sys.stdout
sys.stdout = fn1
print('')
print(u'// Scheme: {0}'.format(scheme))
print('')

for q in range(len(species_mu)):
    print(specm[q][0])
    print('{')
    print('Mu1	{0}	;'.format(specm[q][1]))
    print('Mu2	{0}	;'.format(specm[q][2]))
    print('Mu3	{0}	;'.format(specm[q][3]))
    print('Mu4	{0}	;'.format(specm[q][4]))
    print('}')
sys.stdout = orig_stdout
fn1.close()

with open("./original_form/binaryDiff") as f:
    with open("binaryDiff", "w") as fn1:
        for line in f:
            fn1.write(line)
            if "}" in line:
                break
fn1=open("binaryDiff", "a+")
orig_stdout = sys.stdout
sys.stdout = fn1
print('')
print(u'// Scheme: {0}'.format(scheme))
print('')

for r in range(len(binary_diff)):
    print('{0}-{1}'.format(bindf[r][0], bindf[r][1]))
    print('{')
    print('Diff1	{0}	;'.format(bindf[r][2]))
    print('Diff2	{0}	;'.format(bindf[r][3]))
    print('Diff3	{0}	;'.format(bindf[r][4]))
    print('Diff4	{0}	;'.format(bindf[r][5]))
    print('}')
sys.stdout = orig_stdout
fn1.close()


with open("./original_form/thermoDiff") as f:
    with open("thermoDiff", "w") as fn1:
        for line in f:
            fn1.write(line)
            if "}" in line:
                break
fn1=open("thermoDiff", "a+")
orig_stdout = sys.stdout
sys.stdout = fn1
print('')
print(u'// Scheme: {0}'.format(scheme))
print('')

for s in range(len(therm_diff)):
    print('{0}-{1}'.format(thermdf[s][0], thermdf[s][1]))
    print('{')
    print('ThermDiff_1	{0}	;'.format(thermdf[s][2]))
    print('ThermDiff_2	{0}	;'.format(thermdf[s][3]))
    print('ThermDiff_3	{0}	;'.format(thermdf[s][4]))
    print('ThermDiff_4	{0}	;'.format(thermdf[s][5]))
    print('}')
sys.stdout = orig_stdout
fn1.close()

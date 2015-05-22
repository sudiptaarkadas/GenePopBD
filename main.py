import sys

from classes.load import Load
from classes.forensic import Forensic

input_file = sys.argv[1]

csv_file = Load(input_file)

data = csv_file.load()

r_locas = ['Locas']
r_homo = ['Homo']
r_hetero = ['Hetero']
r_per_homo = ['Per Homo']
r_per_hetero = ['Per Hetero']

r_poe = ['Power of Exclusion']
r_tpi = ['TPI']
r_allele = ['Allele']
r_pic = ['PIC           ']
r_mp = ['Matching Probability']
r_pod = ['PoD           ']

for j in range(1, len(data[0])-1,2):
  allele_l = []
  allele_r = []
  allele_g = []
  allele_all = []

  for i in range(1,len(data)):
    allele_l.append(data[i][j])
    allele_r.append(data[i][j+1])
    if data[i][j] != '?' and data[i][j+1] != '?':
      if float(data[i][j]) > float(data[i][j+1]):
        allele_g.append(data[i][j+1]+','+data[i][j])
      else:
        allele_g.append(data[i][j]+','+data[i][j+1])
    if data[i][j] != '?':
      allele_all.append(data[i][j])
    if data[i][j+1] != '?':
      allele_all.append(data[i][j+1])

  'Initiating forensic class'

  forensic = Forensic(allele_l, allele_r, allele_g, allele_all)

  per_homo_hetero = forensic.per_homo_hetero()

  print '\n'
  print 'Locus: ', data[0][j]
  print '--------------------------'
  print '\n \t Summary'
  print '\t --------------------------'
  print '\t Total = ', per_homo_hetero[0]
  print '\t Homos = ', per_homo_hetero[1], '\t Percent = ', per_homo_hetero[3] 
  print '\t Heteros = ', per_homo_hetero[2], '\t Percent = ', per_homo_hetero[4]

  'Saving locas, homo, hetero in array'

  r_locas.append(data[0][j])
  r_homo.append(per_homo_hetero[1])
  r_per_homo.append(per_homo_hetero[3])
  r_hetero.append(per_homo_hetero[2])
  r_per_hetero.append(per_homo_hetero[4])

  'Saving Power of exclusion and Typical Peternity Index in array'

  pat_stats = forensic.paternity_statistics()

  print '\n \t Paternity Statistics'
  print '\t --------------------------'
  print '\t Power of exclusion = ', pat_stats[0]
  print '\t Typical paternity index = ', pat_stats[1]

  r_poe.append(pat_stats[0])
  r_tpi.append(pat_stats[1])

  'Saving Allele and PIC in array'

  freq_calc = forensic.freq_calc()

  print '\n \t Frequency Calculations'
  print '\t --------------------------'
  print '\t Total allele = ', freq_calc[0]
  print '\t PIC = ', freq_calc[6]
  print '\t Table for frequency '
  print '\t --------------------------'
  print '\n \tAllele       \tFrequency  \tRelative Frequency\tPi2            \t\tPi4'

  for k in range(len(freq_calc[1])):
    print '\t',freq_calc[1][k], '       \t', freq_calc[2][k], '        \t', freq_calc[3][k], '     \t', freq_calc[4][k], '     \t', freq_calc[5][k]

  
  r_allele.append(freq_calc[1])
  r_pic.append(freq_calc[6])

  'Saving Maching Probability and PoD in array'

  geno_calc = forensic.geno_calc()

  print '\n \t Genotype Calculations'
  print '\t --------------------------'
  print '\t Total Genotypes = ', geno_calc[0]
  print '\t Matching Probability = ', geno_calc[4]
  print '\t Power of Discrimination = ', geno_calc[5]
  print '\t Table for frequency '
  print '\t --------------------------'
  print '\n \tGenotype     \tFrequency  \tRelative Frequency'

  for k in range(len(geno_calc[1])):
    print '\t',geno_calc[1][k], '     \t', geno_calc[2][k], '          \t', geno_calc[3][k]


  r_mp.append(geno_calc[4])
  r_pod.append(geno_calc[5])

print '\n\nShort description'
print '------------------'
print 'Homo: Total number of homogeneus'
print 'Per Homo: Percentage of homogeneus'
print 'Hetero: Total number of Heterogeneus'
print 'Per Hetero: Percentage of Heterogeneus'
print 'PoE: Power of Exclusion'
print 'TPI: Typical Paternity Index'
print '\n'

for i in range(len(r_locas)):
  print r_locas[i], '    \t', r_poe[i],  '    \t', r_pic[i], '   \t', r_mp[i], '            \t',r_pod[i], '     \t', r_tpi[i]

'Writting output file in csv format'
csv_file.down([r_locas, r_homo, r_per_homo, r_hetero, r_per_hetero, r_poe, r_tpi, r_pic, r_mp, r_pod])

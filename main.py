import sys

from classes.load import Load
from classes.forensic import Forensic

input_file = sys.argv[1]

csv_file = Load(input_file)

data = csv_file.load()

for j in range(1, len(data[0])-1,2):
  allele_l = []
  allele_r = []

  for i in range(1,len(data)):
    allele_l.append(data[i][j])
    allele_r.append(data[i][j+1])

  forensic = Forensic(allele_l, allele_r)

  per_homo_hetero = forensic.per_homo_hetero()
  
  print '\n'
  print 'Locus: ', data[0][j]
  print '--------------------------'
  print '\n \t Summary'
  print '\t --------------------------'
  print '\t Total = ', per_homo_hetero[0]
  print '\t Homos = ', per_homo_hetero[1], '\t Percent = ', per_homo_hetero[3] 
  print '\t Heteros = ', per_homo_hetero[2], '\t Percent = ', per_homo_hetero[4]

  pat_stats = forensic.paternity_statistics()

  print '\n \t Paternity Statistics'
  print '\t --------------------------'
  print '\t Power of exclusion = ', pat_stats[0]
  print '\t Typical paternity index = ', pat_stats[1]

  freq_calc = forensic.freq_calc()
  
  print '\n \t Frequency Calculations'
  print '\t --------------------------'
  print '\t Total allele = ', freq_calc[0]
  print '\t Table for frequency '
  print '\t --------------------------'
  print '\n \t Allele \t Frequency \t Percent'

  for k in range(len(freq_calc[1])):
    print '\t',freq_calc[1][k], '\t', freq_calc[2][k], '\t', freq_calc[3][k]

  'print (j+1)/2, forensic.per_homo_hetero(), forensic.foren_param(), forensic.freq_calc()'


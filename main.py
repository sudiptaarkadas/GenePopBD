import sys

from classes.load import Load
from classes.forensic import Forensic

input_file = sys.argv[1]

csv_file = Load(input_file)

data = csv_file.load()

for j in range(1, len(data[0])-1,2):
  allele_l = []
  allele_r = []
  allele_g = []

  for i in range(1,len(data)):
    allele_l.append(data[i][j])
    allele_r.append(data[i][j+1])
    allele_g.append(data[i][j]+','+data[i][j+1])

  forensic = Forensic(allele_l, allele_r, allele_g)

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
  print '\t PIC = ', freq_calc[6]
  print '\t Table for frequency '
  print '\t --------------------------'
  print '\n \t Allele \t Frequency \t Percent \t Pi2 \t Pi4'

  for k in range(len(freq_calc[1])):
    print '\t',freq_calc[1][k], '\t', freq_calc[2][k], '\t', freq_calc[3][k], '\t', freq_calc[4][k], '\t', freq_calc[5][k]


  geno_calc = forensic.geno_calc()
  
  print '\n \t Genotype Calculations'
  print '\t --------------------------'
  print '\t Total Genotypes = ', geno_calc[0]
  print '\t Matching Probability = ', geno_calc[4]
  print '\t Power of Discrimination = ', geno_calc[5]
  print '\t Table for frequency '
  print '\t --------------------------'
  print '\n \t Genotype \t Frequency \t Percent'

  for k in range(len(geno_calc[1])):
    print '\t',geno_calc[1][k], '\t', geno_calc[2][k], '\t', geno_calc[3][k]

  'print (j+1)/2, forensic.per_homo_hetero(), forensic.foren_param(), forensic.freq_calc()'

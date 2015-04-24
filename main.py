from classes.load import Load
from classes.forensic import Forensic

input_file = 'input.csv'

csv_file = Load(input_file)

data = csv_file.load()

for j in range(1, len(data[0])-1,2):
  allele_l = []
  allele_r = []

  for i in range(1,len(data)):
    allele_l.append(data[i][j])
    allele_r.append(data[i][j+1])

  forensic = Forensic(allele_l, allele_r)

  print (j+1)/2, forensic.per_homo_hetero(), forensic.foren_param(), forensic.freq_calc()


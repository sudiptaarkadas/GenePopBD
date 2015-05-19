class Forensic(object):
  '''
  Forensic parameter calculation.
  Methods:
    per_homo_hetero: calculate the homo and hetero in percent
    find_homo_hetero: calculate the number of homo and hetero in allele
    peternity_statistics: calculate PoE and TPI
    freq_calc: calculate the frequency of the allele
    geno_calc: calculate the PIC
    check: function to check the errors in input file
  '''

  def __init__(self, a, b, c):
    self.a = a
    self.b = b
    self.c = c

    self.total = len(self.a)
    self.homo = 0
    self.hetero = 0
    self.total_allele = 0
    self.per_homo = 0
    self.per_hetero = 0
    self.total_genotype = 0
    self.check()

  def per_homo_hetero(self):
    self.find_homo_hetero()

    self.per_homo = self.homo * 1.0 / self.total
    self.per_hetero = 1.0 - self.per_homo

    return self.total, self.homo, self.hetero, self.per_homo, self.per_hetero

  def find_homo_hetero(self):
    for i in range(len(self.a)):
      if self.a[i] == self.b[i]:
        self.homo += 1
    
    self.hetero = self.total - self.homo

  def paternity_statistics(self):
    PoE = (self.per_hetero**2) * (1.0 - (2.0 * self.per_hetero * (self.per_homo**2)))
    TPI = 1 / (2 * self.per_homo)

    return PoE, TPI

  def freq_calc(self):
    dom = self.a + self.b

    d = {x:dom.count(x) for x in dom}
  
    sum_occ = 0
    for x in d:
      sum_occ += d[x]

    self.total_allele = sum_occ

    alle = []
    number = []
    percent = []
    Pi2 = []
    Pi4 = []
    sum_Pi2 = 0.0
    sum_Pi4 = 0.0
    pic = 1.0
    for x in sorted(d):
      alle.append(x)
      number.append(d[x])
      percent.append(d[x]*1.0/sum_occ)
      Pi2.append((d[x]*1.0/sum_occ)**2)
      Pi4.append((d[x]*1.0/sum_occ)**4)
      sum_Pi2 = sum_Pi2 + ((d[x]*1.0/sum_occ)**2)
      sum_Pi4 = sum_Pi4 + ((d[x]*1.0/sum_occ)**4)

    pic = 1.0-sum_Pi2-(sum_Pi2**2)+sum_Pi4

    return sum_occ, alle, number, percent, Pi2, Pi4, pic

  def geno_calc(self):
    dom = self.c

    d = {x:dom.count(x) for x in dom}
  
    sum_geno = 0
    for x in d:
      sum_geno += d[x]

    self.total_genotype = sum_geno

    genotype = []
    number = []
    percent = []
    matching_probability = 0.0
    power_of_discrimination = 1.0
    for x in sorted(d):
      genotype.append(x)
      number.append(d[x])
      percent.append((d[x]*1.0/sum_geno)**2)
      matching_probability = matching_probability + ((d[x]*1.0/sum_geno)**2)

    power_of_discrimination = 1.0 - matching_probability

    return sum_geno, genotype, number, percent, matching_probability, power_of_discrimination

  'Functions for checking data file'
  def check(self):
    'Bug 1: fix start'
    miss_allele_id = []
    for i in range(self.total):
      if self.a[i] =='' or self.b[i]=='':
        miss_allele_id.append(i)
    
    for x in miss_allele_id:
      print x
      print self.a.pop(x)
      print self.b.pop(x)
      self.total = self.total - 1

    'But 1: fix end'

    if len(self.a) != len(self.b):
      print "Alleles are not equal."
      return

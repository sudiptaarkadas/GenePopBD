class Forensic(object):
  '''
  Forensic parameter calculation.
  '''

  def __init__(self, a, b):
    self.a = a
    self.b = b

    self.check()

    self.total = len(self.a)
    self.homo = 0
    self.hetero = 0

  def per_homo_hetero(self):
    self.find_homo_hetero()

    per_homo = self.homo * 100.0 / self.total
    per_hetero = 100.0 - per_homo

    return self.total, self.homo, self.hetero, per_homo, per_hetero

  def find_homo_hetero(self):
    for i in range(len(self.a)):
      if self.a[i] == self.b[i]:
        self.homo += 1
    
    self.hetero = self.total - self.homo

  def foren_param(self):
    PoE = self.hetero*self.hetero*(1.0-2.0*self.hetero*
      self.homo*self.homo)
    TPI = 0.5/self.hetero

    return PoE, TPI

  def freq_calc(self):
    dom = self.a + self.b

    d = {x:dom.count(x) for x in dom}
    
    sum_occ = 0
    for x in d:
      sum_occ += d[x]

    alle = []
    number = []
    percent = []
    for x in d:
      alle.append(x)
      number.append(d[x])
      percent.append( d[x]*100.0/sum_occ)

    return alle, number, percent

  'Functions for checking data file'
  def check(self):
    if len(self.a) != len(self.b):
      print "Alleles are not equal."
      return

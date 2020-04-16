from . import model

class SIR(Model):

  def __init__(self):

      s, i, r = sp.Symbol('s'), sp.Symbol('i'), sp.Symbol('r')

      ds = s - (0.34*s*i)*0.1;
      di = i + (0.34*s*i - 0.05*i)*0.1;
      dr = r + 0.05*i*0.1;

      self.dyns = sp.Matrix([ds, di, dr])
      self.vars = [s,i,r] #In predetermined order
      sys_dim = len(self.vars)

      num_direct = 3
      num_temps = 1

      L = np.zeros(num_direct, sys_dim)
      T = np.zeros(num_temps, sys_dim)

      L[0][0] = 1;  #[1 0 0 ]^T
      L[1][1] = 1;  #[0 1 0 ]^T
      L[2][2] = 1;  #[0 0 1 ]^T

      self.offu = np.zeros(num_direct)
      self.offp = np.zeros(num_direct)

      offp[0] = 0.8
      offm[0] = -0.79

      offp[1] = 0.2
      offm[1] = -0.19

      offp[2] = 0.0001
      offm[2] = -0.000099

      b = Bundle(T, L, offu, offl, vars, dyns)

      super().__init__(b)

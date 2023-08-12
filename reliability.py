class Reliability():
  def __init__(self, beam, variables=None):
    self.beam = beam
    self.n_points = 100
    self.set_original(True)
    self.variables = variables
    self.tetha_r = 1
    self.inverted = False

  def check_ELU(self, load, n_points=500):
    if self.inverted:
      max_moment = abs(self.beam.get_max_moment(n_points=n_points, is_inverted=self.inverted))
    else:
      max_moment = self.beam.get_max_moment(n_points=n_points)
    print(f'Momento_rd = {max_moment*1e-3:.1f} kN.m')
    if load < max_moment:
      print('Safe')
    else:
      print('Not safe')

  def set_original(self,update=False):
    if not update:
      i=0
      for section in self.beam.section:
        if isinstance(section,sa.geometry.Rebar):
          section.area = self.rebar_area[i]
          section.center[1] = self.rebar_pos[i]
          i+=1
    else:
      self.rebar_area = []
      self.rebar_pos = []
      for section in self.beam.section:
        if isinstance(section,sa.geometry.Rebar):
          self.rebar_area.append(section.area)
          self.rebar_pos.append(section.center[1])

  def update_beam(self,variable,value,reset=False):
    if(variable=='fc'):
      for section in self.beam.section:
        if isinstance(section,sa.geometry.RectSection):
          section.material.fc = value
    elif(variable=='fy'):
      for section in self.beam.section:
        if type(section) == sa.geometry.Rebar:
          section.material.fy = value
    elif(variable=='fpt'):
      for section in self.beam.section:
        if type(section) == sa.geometry.Tendon:
          section.material.ft = value
          section.material.fy = 0.85 * value # JCSS
    elif(variable=='Es'):
      for section in self.beam.section:
        if isinstance(section,sa.geometry.Rebar):
          section.material.young = value
    elif(variable=='As'):
      for section in self.beam.section:
        if type(section) == sa.geometry.Rebar:
          section.area = value * section.area
    elif(variable=='Ap'):
      for section in self.beam.section:
        if type(section) == sa.geometry.Tendon:
          section.area = value * section.area
    elif(variable=='cover_bottom'):
      for section in self.beam.section:
        if isinstance(section,sa.geometry.Rebar):
          if section.center[1] < self.beam.centroid[1]:
            section.center[1] = value + section.center[1]
    elif(variable=='cover_top'):
      for section in self.beam.section:
        if isinstance(section,sa.geometry.Rebar):
          if section.center[1] > self.beam.centroid[1]:
            section.center[1] = -value + section.center[1]
    elif(variable=='b_web'):
      self.beam.section[0].width = value
    elif(variable=='h_web'):
      if self.inverted:
        top = self.beam.section[0].boundary[1][1]
        self.beam.section[0].height = value
        self.beam.section[0].center[1] = top - value/2
      else:
        self.beam.section[0].height = value
        self.beam.section[0].center[1] = value/2
        height = self.beam.section[1].height
        self.beam.section[1].center[1] = value + height/2
        if isinstance(self.beam.section[2],sa.geometry.RectSection):
          top = self.beam.section[1].boundary[1][1]
          height = self.beam.section[2].height
          self.beam.section[2].center[1] = top + height/2
    elif(variable=='h_flange'):
      if self.inverted:
        top = self.beam.section[1].boundary[1][1]
        self.beam.section[1].height = value
        self.beam.section[1].center[1] = top - value/2
        bottom = self.beam.section[1].boundary[0][1]
        top = self.beam.section[0].boundary[1][1] 
        self.beam.section[0].center[1] -= top-bottom
      else:
        if isinstance(self.beam.section[2],sa.geometry.RectSection):
          self.beam.section[2].height = value
          bottom = self.beam.section[1].boundary[1][1]
          self.beam.section[2].center[1] = bottom + value/2
        else:
          self.beam.section[1].height = value
          bottom = self.beam.section[0].boundary[1][1]
          self.beam.section[1].center[1] = bottom + value/2
    #elif(variable=='tetha_r'):
      #self.tetha_r = value



  def get_resistance_moment(self,n_points=None, error=1e3, max_increment_e0=0.0005):
    height = self.beam.get_section_boundary()[1][1]
    bottom = self.beam.get_section_boundary()[0][1]
    et, eb = -4e-3, 15e-3
    max_curvature = (eb - et) / (height - bottom)

    if n_points == None:
      n_points=self.n_points

    try:
      if self.inverted:
        return self.beam.get_max_moment(n_points=n_points, is_inverted=self.inverted, error=1e3, max_increment_e0=0.0005)
        #return self.beam.get_moment_curvature(max_curvature, normal_force=0, n_points=500)
        #return abs(self.beam.get_max_moment(n_points, inverted=self.inverted, error=0.01))
        #return self.beam.get_max_moment_reversed_simplified()
      else:
        return self.beam.get_max_moment(n_points, error=1e3, max_increment_e0=0.0005)
        #return self.beam.get_max_moment_simplified()
    except TypeError:
        return 0


  def get_bk(self,variable):
    pos = (self.variables['names'].index(variable))
    mean = self.variables['bounds'][pos][0]
    if mean == 0:
      return 0
    std = self.variables['bounds'][pos][1]
    self.update_beam(variable,mean)
    ym = self.get_resistance_moment()
    bk = 0
    for i,d_x in enumerate(std*np.array([-3,-2,-1,1,2,3])):
      self.update_beam(variable,mean+d_x)
      y = self.get_resistance_moment()
      d_y = y - ym
      bk += abs((d_y/ym)/(d_x/mean))
    self.update_beam(variable,mean)
    cov = std/mean
    return bk*cov

  def sensitivity_bk(self):
    bk = []
    for variable in self.variables['names']:
      bk.append(self.get_bk(variable))
    return np.array(bk)/max(bk)*100

  def sensitivity_Sobol(self, n=512):
    values_sobol = sobol_sample.sample(variables, n)
    Y_sobol = []
    for value in values_sobol:
      self.set_original()
      for variable in self.variables['names']:
        pos = (self.variables['names'].index(variable))
        self.update_beam(variable,value[pos])
      Y_sobol.append(self.get_resistance_moment())
    self.set_original()
    #self.beam._compute_centroid()
    Y_sobol = np.array(Y_sobol)
    Si = sobol.analyze(variables, Y_sobol)
    for i in range(len(self.variables['names'])):
      fig, graph = plt.subplots(1,2)
      graph[0].plot(values_sobol.T[i], Y_sobol/1e3,'o')
      graph[1].hist(values_sobol.T[i])
      graph[0].set(xlabel =self.variables['names'][i])
    return Si['S1'], Si['ST']

  def monte_carlo_simplified(self,n_LHS=100000):
    param_values = latin.sample(self.variables, n_LHS)
    moment_simplified =[]
    
    for values in param_values:
        self.set_original()
        for j,value in enumerate(values):
          self.update_beam(self.variables['names'][j],value)
        self.beam._compute_centroid()
        moment_simplified.append(self.get_resistance_moment())
    
    moment_simplified = np.array(moment_simplified)
    moment_simplified = moment_simplified[moment_simplified>0]
    return moment_simplified

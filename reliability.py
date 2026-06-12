import numpy as np
import secan as sa
from SALib.analyze import sobol
from SALib.sample import latin
from SALib.sample import sobol as sobol_sample

class Reliability():
  def __init__(self, beam, variables=None):
    self.beam = beam
    self.n_points = 100
    self.error = 1e3
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
      for section, original in zip(self.beam.section, self.original_sections):
        section.center = original['center'].copy()
        if isinstance(section,sa.geometry.RectSection):
          section.width = original['width']
          section.height = original['height']
        if isinstance(section,sa.geometry.Rebar):
          section.area = original['area']
          section.material.young = original['young']
          section.material.fy = original['fy']
        if isinstance(section,sa.geometry.Tendon):
          section.material.ft = original['ft']
        if hasattr(section.material, 'fc'):
          section.material.fc = original['fc']
      self.beam._compute_centroid()
    else:
      self.original_sections = []
      self.rebar_area = []
      self.rebar_pos = []
      for section in self.beam.section:
        original = {
          'center': section.center.copy()
        }
        if isinstance(section,sa.geometry.RectSection):
          original['width'] = section.width
          original['height'] = section.height
        if isinstance(section,sa.geometry.Rebar):
          original['area'] = section.area
          original['young'] = section.material.young
          original['fy'] = section.material.fy
          self.rebar_area.append(section.area)
          self.rebar_pos.append(section.center[1])
        if isinstance(section,sa.geometry.Tendon):
          original['ft'] = section.material.ft
        if hasattr(section.material, 'fc'):
          original['fc'] = section.material.fc
        self.original_sections.append(original)
      self.original_centroid = self.beam.centroid.copy()

  def update_beam(self,variable,value,reset=False):
    handlers = {
      'fc': self._update_fc,
      'fy': self._update_fy,
      'fpt': self._update_fpt,
      'Es': self._update_Es,
      'As': self._update_As,
      'Ap': self._update_Ap,
      'cover_bottom': self._update_cover_bottom,
      'cover_top': self._update_cover_top,
      'b_web': self._update_b_web,
      'h_web': self._update_h_web,
      'h_flange': self._update_h_flange
    }
    if variable in handlers:
      handlers[variable](value)

  def _update_fc(self,value):
    for section in self.beam.section:
      if isinstance(section,sa.geometry.RectSection):
        section.material.fc = value

  def _update_fy(self,value):
    for section in self.beam.section:
      if type(section) == sa.geometry.Rebar:
        section.material.fy = value

  def _update_fpt(self,value):
    for section in self.beam.section:
      if type(section) == sa.geometry.Tendon:
        section.material.ft = value
        section.material.fy = 0.85 * value # JCSS

  def _update_Es(self,value):
    for section in self.beam.section:
      if isinstance(section,sa.geometry.Rebar):
        section.material.young = value

  def _update_As(self,value):
    for i,section in enumerate(self.beam.section):
      if type(section) == sa.geometry.Rebar:
        section.area = value * self.original_sections[i]['area']

  def _update_Ap(self,value):
    for i,section in enumerate(self.beam.section):
      if type(section) == sa.geometry.Tendon:
        section.area = value * self.original_sections[i]['area']

  def _update_cover_bottom(self,value):
    for i,section in enumerate(self.beam.section):
      if isinstance(section,sa.geometry.Rebar):
        if self.original_sections[i]['center'][1] < self.original_centroid[1]:
          section.center[1] = value + self.original_sections[i]['center'][1]

  def _update_cover_top(self,value):
    for i,section in enumerate(self.beam.section):
      if isinstance(section,sa.geometry.Rebar):
        if self.original_sections[i]['center'][1] > self.original_centroid[1]:
          section.center[1] = -value + self.original_sections[i]['center'][1]

  def _update_b_web(self,value):
    self.beam.section[0].width = value

  def _update_h_web(self,value):
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

  def _update_h_flange(self,value):
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



  def get_resistance_moment(self,n_points=None):
    height = self.beam.get_section_boundary()[1][1]
    bottom = self.beam.get_section_boundary()[0][1]
    et, eb = -4e-3, 15e-3
    max_curvature = (eb - et) / (height - bottom)

    if n_points == None:
      n_points=self.n_points

    try:
      if self.inverted:
        return self.beam.get_max_moment(n_points=n_points, is_inverted=self.inverted, error=self.error)
        #return self.beam.get_moment_curvature(max_curvature, normal_force=0, n_points=500)
        #return abs(self.beam.get_max_moment(n_points, inverted=self.inverted, error=0.01))
        #return self.beam.get_max_moment_reversed_simplified()
      else:
        return self.beam.get_max_moment(n_points, error=self.error)
        #return self.beam.get_max_moment_simplified()
    except TypeError:
        return 0


  def get_bk(self,variable):
    pos = (self.variables['names'].index(variable))
    mean = self.variables['bounds'][pos][0]
    if mean == 0:
      return 0
    std = self.variables['bounds'][pos][1]
    self.set_original()
    self.update_beam(variable,mean)
    self.beam._compute_centroid()
    ym = self.get_resistance_moment()
    bk = 0
    for i,d_x in enumerate(std*np.array([-3,-2,-1,1,2,3])):
      self.set_original()
      self.update_beam(variable,mean+d_x)
      self.beam._compute_centroid()
      y = self.get_resistance_moment()
      d_y = y - ym
      bk += abs((d_y/ym)/(d_x/mean))
    self.set_original()
    cov = std/mean
    return bk*cov

  def sensitivity_bk(self):
    bk = []
    for variable in self.variables['names']:
      bk.append(self.get_bk(variable))
    self.set_original()
    return np.array(bk)/max(bk)*100

  def sensitivity_Sobol(self, n=512, seed=None, calc_second_order=False):
    values_sobol = sobol_sample.sample(
      self.variables,
      n,
      calc_second_order=calc_second_order,
      seed=seed
    )
    Y_sobol = []
    for value in values_sobol:
      self.set_original()
      for variable in self.variables['names']:
        pos = (self.variables['names'].index(variable))
        self.update_beam(variable,value[pos])
      self.beam._compute_centroid()
      Y_sobol.append(self.get_resistance_moment())
    self.set_original()
    self.beam._compute_centroid()
    Y_sobol = np.array(Y_sobol)
    Si = sobol.analyze(
      self.variables,
      Y_sobol,
      calc_second_order=calc_second_order,
      seed=seed
    )
    return Si['S1'], Si['ST']

  def monte_carlo_simplified(self,n_LHS=100000, seed=None):
    param_values = latin.sample(self.variables, n_LHS, seed=seed)
    moment_simplified =[]
    
    for values in param_values:
        self.set_original()
        for j,value in enumerate(values):
          self.update_beam(self.variables['names'][j],value)
        self.beam._compute_centroid()
        moment_simplified.append(self.get_resistance_moment())
    self.set_original()
    self.beam._compute_centroid()
    
    moment_simplified = np.array(moment_simplified)
    moment_simplified = moment_simplified[moment_simplified != 0]
    return moment_simplified

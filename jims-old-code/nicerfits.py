from math import *
import biggles, Numeric, cPickle as pickle
import gbwkf
import gbwkftau

allthat = pickle.load(file("/home/mccann/antithesis/novemberdata.p"))
u1runs = allthat["u1runs"]
u2runs = allthat["u2runs"]
u3runs = allthat["u3runs"]
u1data = allthat["u1data"]
u2data = allthat["u2data"]
u3data = allthat["u3data"]

class FitRecord: pass
fitrecord = pickle.load(file("/home/mccann/antithesis/fit_results/novemberfits_noapr03_3_1.0.p"))

class RunSummary : pass
initialrunlist = pickle.load(open("/home/mccann/synthesis/lumieff/initialrunlist.p"))
runsummary = pickle.load(open("/home/mccann/synthesis/lumieff/runsummary.p"))
runsummary[123828].energy = 4.72992
runsummary[123832].energy = 4.72990

def mygbwkf(mass, fullgam, rmsbeam, yint, phi, w):
  "yint = 0.018, 0.018, 0.018; phi=0"
  if w > mass + 200.:
    return 0.076/(w-mass)
  return gbwkf.gbwkf(mass, fullgam, rmsbeam, yint, phi, w-mass)

def mygbwkftau(mass, fullgam, rmsbeam, yint, phi, w):
  "yint = 0.20, 0.37, 0.27; phi = 0"
  if w > mass + 200.:
    return 0.076/(w-mass)
  return gbwkftau.gbwkf(mass, fullgam, rmsbeam, yint, phi, w-mass)

def u1func(area, rmsbeam, back, fullgam, yint, phi, btautau, tauyint, tauphi, twophofrac, w):
  tmp = 0.
  tmp += area * 0.9793 * mygbwkf(9460.30, fullgam, rmsbeam, yint, phi, w)
  tmp += area * 0.578 * btautau * mygbwkftau(9460.30, fullgam, rmsbeam, tauyint, tauphi, w)
  tmp += back * (1.-twophofrac) * 9000.**2 / w**2
  tmp += back * twophofrac * log(w**2/9000.**2)
  return tmp

def u2func(area, rmsbeam, back, fullgam, yint, phi, btautau, tauyint, tauphi, twophofrac, u1area, w):
  tmp = 0.
  tmp += area * 0.9618 * mygbwkf(10023.26, fullgam, rmsbeam, yint, phi, w)
  tmp += area * 0.578 * btautau * mygbwkftau(10023.26, fullgam, rmsbeam, tauyint, tauphi, w)
  tmp += back * (1.-twophofrac) * 9000.**2 / w**2
  tmp += back * twophofrac * log(w**2/9000.**2)
  tmp += u1area * mygbwkf(9460.30, 0., 0., 0., 0., w)
  return tmp

def u3func(area, rmsbeam, back, fullgam, yint, phi, btautau, tauyint, tauphi, twophofrac, u1area, u2area, w):
  tmp = 0.
  tmp += area * 0.9641 * mygbwkf(10355.2, fullgam, rmsbeam, yint, phi, w)
  tmp += area * 0.578 * btautau * mygbwkftau(10355.2, fullgam, rmsbeam, tauyint, tauphi, w)
  tmp += back * (1.-twophofrac) * 9000.**2 / w**2
  tmp += back * twophofrac * log(w**2/9000.**2)
  tmp += u1area * mygbwkf(9460.30, 0., 0., 0., 0., w)
  tmp += u2area * mygbwkf(10023.26, 0., 0., 0., 0., w)
  return tmp

def u1func_bkgndonly(area, rmsbeam, back, fullgam, yint, phi, btautau, tauyint, tauphi, twophofrac, w):
  tmp = 0.
  tmp += area * 0.578 * btautau * mygbwkftau(9460.30, fullgam, rmsbeam, tauyint, tauphi, w)
  tmp += back * (1.-twophofrac) * 9000.**2 / w**2
  tmp += back * twophofrac * log(w**2/9000.**2)
  return tmp

def u2func_bkgndonly(area, rmsbeam, back, fullgam, yint, phi, btautau, tauyint, tauphi, twophofrac, u1area, w):
  tmp = 0.
  tmp += area * 0.578 * btautau * mygbwkftau(10023.26, fullgam, rmsbeam, tauyint, tauphi, w)
  tmp += back * (1.-twophofrac) * 9000.**2 / w**2
  tmp += back * twophofrac * log(w**2/9000.**2)
  tmp += u1area * mygbwkf(9460.30, 0., 0., 0., 0., w)
  return tmp

def u3func_bkgndonly(area, rmsbeam, back, fullgam, yint, phi, btautau, tauyint, tauphi, twophofrac, u1area, u2area, w):
  tmp = 0.
  tmp += area * 0.578 * btautau * mygbwkftau(10355.2, fullgam, rmsbeam, tauyint, tauphi, w)
  tmp += back * (1.-twophofrac) * 9000.**2 / w**2
  tmp += back * twophofrac * log(w**2/9000.**2)
  tmp += u1area * mygbwkf(9460.30, 0., 0., 0., 0., w)
  tmp += u2area * mygbwkf(10023.26, 0., 0., 0., 0., w)
  return tmp

def mean(xlist):
  """Takes a list of values and returns the mean."""
  s = 0.
  for x in xlist:
    s += x
  return s/float(len(xlist))

def rms(xlist):
  """Takes a list of values and returns the root mean square (not standard deviation)."""
  s2 = 0.
  for x in xlist:
    s2 += x**2
  return sqrt(s2/float(len(xlist)))

def stdev(xlist):
  """Takes a list of values and returns the standard deviation."""
  s = 0.
  s2 = 0.
  for x in xlist:
    s += x
    s2 += x**2
  return sqrt(s2/float(len(xlist)) - (s/float(len(xlist)))**2)

def wmean(xlist):
  """Takes a list of (value, error) pairs and returns (weighted mean, its error).  Values with non-positive errors are ignored (dropped from the list)."""
  s = 0.
  w = 0.
  for (x,e) in xlist:
    if e > 0:
      wi = 1.0/float(e)**2
      s += x*wi
      w += wi
  return s/w, sqrt(1.0/w)

def adddata(p, runs, data, shift):
  x = []
  y = []
  yerr = []
  allpeak = []
  themap = {}
  for r, (e, h, herr) in zip(runs, data):
    if r != None and runsummary[r].kind == "p":
      allpeak.append((e, h, herr))
    else:
      thetag = int(round(e*2000.*10.))
      if thetag not in themap:
        themap[thetag] = []
      themap[thetag].append((e,h,herr))
  for thetag, thelist in themap.items():
    some_energy = mean(map(lambda (e,h,herr): e, thelist)) * 2. + shift/1000.
    some_h, some_herr = wmean(map(lambda (e,h,herr): (h,herr), thelist))
    x.append(some_energy)
    y.append(some_h)
    yerr.append(some_herr)
  if allpeak != []:
    allpeak_energy = mean(map(lambda (e,h,herr): e, allpeak)) * 2. + shift/1000.
    allpeak_h, allpeak_herr = wmean(map(lambda (e,h,herr): (h,herr), allpeak))
    x.append(allpeak_energy)
    y.append(allpeak_h)
    yerr.append(allpeak_herr)
  p.add(biggles.Points(x, y, symboltype="filled circle", symbolsize=0.8))
  p.add(biggles.SymmetricErrorBarsY(x, y, yerr))
  return None

def addfunc(p, f, low, high, points=100., linetype="solid", linewidth=1.):
  x = Numeric.arange(low, high+(high-low)/points, (high-low)/points)
  y = Numeric.arange(low, high+(high-low)/points, (high-low)/points)
  for i, xi in enumerate(x):
    y[i] = f(xi)
  tmp = biggles.Curve(x/1000., y, linetype=linetype, linewidth=linewidth)
  p.add(tmp)
  return tmp

def adddata_pull(p, data, shift, f):
  x = []
  y = []
  for (e, h, herr) in data:
    x.append((e*2000.+shift)/1000.)
    y.append((h - f(e*2000.+shift))/herr)
  p.add(biggles.Points(x, y, symboltype="filled circle", symbolsize=0.8))
  return y

p = biggles.Table(1,3)
pull = biggles.Table(1,3)

myarea, myrmsbeam, myback, myjan16, myjan30, myfeb06, myfeb13, myfeb20, myfeb27, mymar06, mymar13, myapr08, myapr09, myapr10, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myrjan, myrfeb, myrapr2 = fitrecord[1].values
thefunc = lambda w: u1func(myarea, myrmsbeam, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, w)
thefunc_bkgnd = lambda w: u1func_bkgndonly(myarea, myrmsbeam, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, w)

q = biggles.FramedPlot()
adddata(q, [None], u1data["high"], 0.)
addfunc(q, thefunc, u1data["high"][0][0]*2000.-5., u1data["high"][0][0]*2000.+5.)
addfunc(q, thefunc_bkgnd, u1data["high"][0][0]*2000.-5., u1data["high"][0][0]*2000.+5., linetype="dashed")
q.x.range = (u1data["high"][0][0]*2000.-5.)/1000., (u1data["high"][0][0]*2000.+5.)/1000.
q.x.ticks = [u1data["high"][0][0]*2.]
q.x.subticks = [u1data["high"][0][0]*2.]
q.x.ticklabels = ["%4.2f" % (u1data["high"][0][0]*2.)]
q.x.draw_subticks = 0
q.x2.draw_ticklabels = 0
q.y.range = 7.6, 8.1
q.y.ticks = [7.6, 7.8, 8.]
q.y2.draw_ticklabels = 0
q.y.draw_subticks = 0

p[0,0] = biggles.FramedPlot()
adddata(p[0,0], [None], u1data["cont"], 0.)
adddata(p[0,0], [None], u1data["high"], 0.)
adddata(p[0,0], u1runs["jan16"], u1data["jan16"], myjan16)
adddata(p[0,0], u1runs["jan30"], u1data["jan30"], myjan30)
adddata(p[0,0], u1runs["feb06"], u1data["feb06"], myfeb06)
adddata(p[0,0], u1runs["feb13"], u1data["feb13"], myfeb13)
adddata(p[0,0], u1runs["feb20"], u1data["feb20"], myfeb20)
adddata(p[0,0], u1runs["feb27"], u1data["feb27"], myfeb27)
adddata(p[0,0], u1runs["mar06"], u1data["mar06"], mymar06)
adddata(p[0,0], u1runs["mar13"], u1data["mar13"], mymar13)
adddata(p[0,0], u1runs["apr08"], u1data["apr08"], myapr08)
adddata(p[0,0], u1runs["apr09"], u1data["apr09"], myapr09)
adddata(p[0,0], u1runs["apr10"], u1data["apr10"], myapr10)
addfunc(p[0,0], thefunc, 9420., 9510.)
addfunc(p[0,0], thefunc_bkgnd, 9445., 9510., linetype="dashed")
p[0,0].add(biggles.PlotInset((0.7, 0.5), (0.9, 0.9), q))
p[0,0].x.range = 9.42, 9.51
p[0,0].x.ticks = [9.42, 9.45, 9.48, 9.51]
p[0,0].x.ticklabels = ["9.42", "9.45", "9.48", "9.51"]
p[0,0].x.subticks = Numeric.arange(9.42, 9.51+0.01, 0.01)
p[0,0].x2.draw_ticklabels = 0

pull[0,0] = biggles.FramedPlot()
adddata_pull(pull[0,0], u1data["cont"], 0., lambda w: u1func(myarea, myrmsbeam, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, w))
adddata_pull(pull[0,0], u1data["high"], 0., lambda w: u1func(myarea, myrmsbeam, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, w))
adddata_pull(pull[0,0], u1data["jan16"], myjan16, lambda w: u1func(myarea, myrjan, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, w))
adddata_pull(pull[0,0], u1data["jan30"], myjan30, lambda w: u1func(myarea, myrjan, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, w))
adddata_pull(pull[0,0], u1data["feb06"], myfeb06, lambda w: u1func(myarea, myrjan, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, w))
adddata_pull(pull[0,0], u1data["feb13"], myfeb13, lambda w: u1func(myarea, myrjan, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, w))
adddata_pull(pull[0,0], u1data["feb20"], myfeb20, lambda w: u1func(myarea, myrjan, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, w))
adddata_pull(pull[0,0], u1data["feb27"], myfeb27, lambda w: u1func(myarea, myrfeb, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, w))
adddata_pull(pull[0,0], u1data["mar06"], mymar06, lambda w: u1func(myarea, myrfeb, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, w))
adddata_pull(pull[0,0], u1data["mar13"], mymar13, lambda w: u1func(myarea, myrfeb, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, w))
adddata_pull(pull[0,0], u1data["apr08"], myapr08, lambda w: u1func(myarea, myrapr2, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, w))
adddata_pull(pull[0,0], u1data["apr09"], myapr09, lambda w: u1func(myarea, myrapr2, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, w))
adddata_pull(pull[0,0], u1data["apr10"], myapr10, lambda w: u1func(myarea, myrapr2, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, w))
pull[0,0].add(biggles.LineY(0.))
pull[0,0].x.range = 9.42, 9.51
pull[0,0].x.ticks = [9.42, 9.45, 9.48, 9.51]
pull[0,0].x.draw_ticklabels = 0
pull[0,0].x.subticks = Numeric.arange(9.42, 9.51+0.01, 0.01)
pull[0,0].y.range = -4., 4.
pull[0,0].y.ticks = [-3., 0., 3.]
pull[0,0].y.subticks = [-4., -3., -2., -1., 0., 1., 2., 3., 4.]
pull[0,0].y2.draw_ticklabels = 0

myarea, myrmsbeam, myback, mymay29, myjun11, myjun12, myjul10, myjul24, myaug07, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area = fitrecord[2].values
thefunc = lambda w: u2func(myarea, myrmsbeam, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, w)
thefunc_bkgnd = lambda w: u2func_bkgndonly(myarea, myrmsbeam, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, w)

q = biggles.FramedPlot()
adddata(q, [None], u2data["high"], 0.)
addfunc(q, thefunc, 10080., 10090.)
addfunc(q, thefunc_bkgnd, 10080., 10090., linetype="dashed")
q.x.range = 10.080, 10.090
q.x.ticks = [u2data["high"][0][0]*2.]
q.x.subticks = [u2data["high"][0][0]*2.]
q.x.ticklabels = ["%5.2f" % (u2data["high"][0][0]*2.)]
q.x.draw_subticks = 0
q.x2.draw_ticklabels = 0
q.y.range = 6.9, 7.5
q.y.ticks = [7., 7.2, 7.4]
q.y2.draw_ticklabels = 0
q.y.draw_subticks = 0

p[0,1] = biggles.FramedPlot()
adddata(p[0,1], [None], u2data["cont"], 0.)
adddata(p[0,1], [None], u2data["high"], 0.)
adddata(p[0,1], u2runs["may29"], u2data["may29"], mymay29)
adddata(p[0,1], u2runs["jun11"], u2data["jun11"], myjun11)
adddata(p[0,1], u2runs["jun12"], u2data["jun12"], myjun12)
adddata(p[0,1], u2runs["jul10"], u2data["jul10"], myjul10)
adddata(p[0,1], u2runs["jul24"], u2data["jul24"], myjul24)
adddata(p[0,1], u2runs["aug07"], u2data["aug07"], myaug07)
addfunc(p[0,1], thefunc, 9980., 10100.)
addfunc(p[0,1], thefunc_bkgnd, 10012, 10100., linetype="dashed")
p[0,1].add(biggles.PlotInset((0.7, 0.5), (0.9, 0.9), q))
p[0,1].x.range = 9.980, 10.100
p[0,1].x.ticks = [9.98, 10.02, 10.06, 10.10]
p[0,1].x.ticklabels = ["9.98", "10.02", "10.06", "10.10"]
p[0,1].x.subticks = Numeric.arange(9.98, 10.10+0.01, 0.01)
p[0,1].x2.draw_ticklabels = 0

pull[0,1] = biggles.FramedPlot()
adddata_pull(pull[0,1], u2data["cont"], 0., lambda w: u2func(myarea, myrmsbeam, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, w))
adddata_pull(pull[0,1], u2data["high"], 0., lambda w: u2func(myarea, myrmsbeam, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, w))
adddata_pull(pull[0,1], u2data["may29"], mymay29, lambda w: u2func(myarea, myrmsbeam, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, w))
adddata_pull(pull[0,1], u2data["jun11"], myjun11, lambda w: u2func(myarea, myrmsbeam, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, w))
adddata_pull(pull[0,1], u2data["jun12"], myjun12, lambda w: u2func(myarea, myrmsbeam, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, w))
adddata_pull(pull[0,1], u2data["jul10"], myjul10, lambda w: u2func(myarea, myrmsbeam, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, w))
adddata_pull(pull[0,1], u2data["jul24"], myjul24, lambda w: u2func(myarea, myrmsbeam, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, w))
adddata_pull(pull[0,1], u2data["aug07"], myaug07, lambda w: u2func(myarea, myrmsbeam, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, w))
pull[0,1].add(biggles.LineY(0.))
pull[0,1].x.range = 9.980, 10.100
pull[0,1].x.ticks = [9.98, 10.02, 10.06, 10.10]
pull[0,1].x.draw_ticklabels = 0
pull[0,1].x.subticks = Numeric.arange(9.98, 10.10+0.01, 0.01)
pull[0,1].y.range = -4., 4.
pull[0,1].y.ticks = [-3., 0., 3.]
pull[0,1].y.subticks = [-4., -3., -2., -1., 0., 1., 2., 3., 4.]
pull[0,1].y2.draw_ticklabels = 0

myarea, myrmsbeam, myback, mynov28, mydec05, mydec12, mydec19, mydec26, myjan02, myjan09, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, myu2area, myrnov28, myrdec05, myrdec12, myrdec19, myrdec26, myrjan02, myrjan09 = fitrecord[3].values
thefunc = lambda w: u3func(myarea, myrmsbeam, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, myu2area, w)
thefunc_bkgnd = lambda w: u3func_bkgndonly(myarea, myrmsbeam, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, myu2area, w)

q = biggles.FramedPlot()
adddata(q, [None], u3data["high"], 0.)
addfunc(q, thefunc, u3data["high"][0][0]*2000.-5., u3data["high"][0][0]*2000.+5.)
addfunc(q, thefunc_bkgnd, u3data["high"][0][0]*2000.-5., u3data["high"][0][0]*2000.+5., linetype="dashed")
q.x.range = (u3data["high"][0][0]*2000.-5.)/1000., (u3data["high"][0][0]*2000.+5.)/1000.
q.x.ticks = [u3data["high"][0][0]*2.]
q.x.subticks = [u3data["high"][0][0]*2.]
q.x.ticklabels = ["%5.2f" % (u3data["high"][0][0]*2.)]
q.x.draw_subticks = 0
q.x2.draw_ticklabels = 0
q.y.range = 6.6, 6.9
q.y.ticks = [6.6, 6.7, 6.8, 6.9]
q.y2.draw_ticklabels = 0
q.y.draw_subticks = 0

p[0,2] = biggles.FramedPlot()
adddata(p[0,2], [None], u3data["cont"], 0.)
adddata(p[0,2], [None], u3data["high"], 0.)
adddata(p[0,2], u3runs["nov28"], u3data["nov28"], mynov28)
adddata(p[0,2], u3runs["dec05"], u3data["dec05"], mydec05)
adddata(p[0,2], u3runs["dec12"], u3data["dec12"], mydec12)
adddata(p[0,2], u3runs["dec19"], u3data["dec19"], mydec19)
adddata(p[0,2], u3runs["dec26"], u3data["dec26"], mydec26)
adddata(p[0,2], u3runs["jan02"], u3data["jan02"], myjan02)
adddata(p[0,2], u3runs["jan09"], u3data["jan09"], myjan09)
addfunc(p[0,2], thefunc, 10320., 10410.)
addfunc(p[0,2], thefunc_bkgnd, 10342.5, 10410., linetype="dashed")
p[0,2].add(biggles.PlotInset((0.7, 0.5), (0.9, 0.9), q))
p[0,2].x.range = 10.320, 10.410
p[0,2].x.ticks = [10.32, 10.35, 10.38, 10.41]
p[0,2].x.ticklabels = ["10.32", "10.35", "10.38", "10.41"]
p[0,2].x.subticks = Numeric.arange(10.32, 10.41+0.01, 0.01)
p[0,2].x2.draw_ticklabels = 0

pull[0,2] = biggles.FramedPlot()
adddata_pull(pull[0,2], u3data["cont"], 0., lambda w: u3func(myarea, myrmsbeam, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, myu2area, w))
adddata_pull(pull[0,2], u3data["high"], 0., lambda w: u3func(myarea, myrmsbeam, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, myu2area, w))
adddata_pull(pull[0,2], u3data["nov28"], mynov28, lambda w: u3func(myarea, myrnov28, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, myu2area, w))
adddata_pull(pull[0,2], u3data["dec05"], mydec05, lambda w: u3func(myarea, myrdec05, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, myu2area, w))
adddata_pull(pull[0,2], u3data["dec12"], mydec12, lambda w: u3func(myarea, myrdec12, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, myu2area, w))
adddata_pull(pull[0,2], u3data["dec19"], mydec19, lambda w: u3func(myarea, myrdec19, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, myu2area, w))
adddata_pull(pull[0,2], u3data["dec26"], mydec26, lambda w: u3func(myarea, myrdec26, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, myu2area, w))
adddata_pull(pull[0,2], u3data["jan02"], myjan02, lambda w: u3func(myarea, myrjan02, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, myu2area, w))
adddata_pull(pull[0,2], u3data["jan09"], myjan09, lambda w: u3func(myarea, myrjan09, myback, myfullgam, myyint, myphi, mybtautau, mytauyint, mytauphi, mytwophofrac, myu1area, myu2area, w))
pull[0,2].add(biggles.LineY(0.))
pull[0,2].x.range = 10.320, 10.410
pull[0,2].x.ticks = [10.32, 10.35, 10.38, 10.41]
pull[0,2].x.draw_ticklabels = 0
pull[0,2].x.subticks = Numeric.arange(10.32, 10.41+0.01, 0.01)
pull[0,2].y.range = -4., 4.
pull[0,2].y.ticks = [-3., 0., 3.]
pull[0,2].y.subticks = [-4., -3., -2., -1., 0., 1., 2., 3., 4.]
pull[0,2].y2.draw_ticklabels = 0

# p[0,0].y1.label = r"Selected events / nb$^{-1}$"
# pull[0,0].y1.label = " "

p[0,0].x1.label = r" "
p[0,1].x1.label = r"$E_{CM} (GeV)$"
p[0,2].x1.label = r" "

## p.aspect_ratio = 1./3.25
## p.write_eps("nicerfits.eps")

## p.aspect_ratio = 2./3.25
## p.write_eps("nicerfits_narrow.eps")

p.aspect_ratio = 0.75/3.25
pull.aspect_ratio = 0.25/3.25
p.write_eps("nicerfits_bottom.eps")
pull.write_eps("nicerfits_top.eps")
pull[0,0].aspect_ratio = 0.25/3.25*3.
pull[0,1].aspect_ratio = 0.25/3.25*3.
pull[0,2].aspect_ratio = 0.25/3.25*3.
pull[0,0].write_eps("nicerfits_top1.eps")
pull[0,1].write_eps("nicerfits_top2.eps")
pull[0,2].write_eps("nicerfits_top3.eps")

#!/usr/bin/env python

import time
from astropy.time import Time,TimeDelta
from astroquery.jplhorizons import Horizons
from optparse import OptionParser

parser = OptionParser(usage="%prog: [options] [outfile (leave blank for stdout]")
parser.add_option("-l", "--location", dest="loc", type=str, default="-72",
                  help="Observer location (can be catalog ID or name) [default=%default]")
parser.add_option("-t", "--target", dest="targetid", type=str, default="Voyager 1",
                  help="Target ID (can be catalog ID or name) [default=%default]")
parser.add_option("-y", "--target-type", dest="targettype", type=str, default="majorbody",
                  help="Target catalog type: can be majorbody, smallbody, designation, name, asteroid_name, or comet_name [default=%default]")
parser.add_option("-s", "--step-size", dest="stepsize", type=str, default="1 hour",
                  help="Step size of output -- e.g. '1 hour', '10.5 minutes' [default=%default]")
parser.add_option("-d", "--duration", dest="duration", type=int, default=24,
                  help="Number of steps to output [default=%default]")
parser.add_option("-b", "--start-time", dest="starttime", type=str, default=None,
                  help="Start time (leave blank for current time)")
parser.add_option("-f", "--time-format", dest="timefmt", type=str, default="mjd",
                  help="Format for start-time option (valid choices are 'jd', 'mjd', 'isot', 'unix', 'fits') [default=%default]")

(options,args) = parser.parse_args()

if options.starttime is not None:
    if options.timefmt in ("isot", "fits"):
        start = Time(options.starttime, format=options.timefmt)
    else:
        start = Time(float(options.starttime), format=options.timefmt)
else:
    start = Time.now()

start.format = "jd"
res = options.stepsize.split()
resolution_secs = float(res[0])*{"second": 1,
                               "minute": 60,
                               "hour": 3600,
                               "day": 36400}[res[1]]

dt = TimeDelta(resolution_secs, format="sec")
times = [start + dt*i for i in range(options.duration)]
jdtimes = [t.value for t in times]

print("Going out to Horizons...")
query = Horizons(id=options.targetid,
                 location=options.loc,
                 epochs=jdtimes,
                 id_type=options.targettype)
eph = query.ephemerides()
print("Finished.")

pointingdata = []
for i in eph:
    invrange = 0.0 #FIXME calc this -- what units?
    mjdtime = Time(i['datetime_jd'], format="jd")
    mjdtime.format = "mjd"
    pointingdata.append("{:.0f}  {:.05f}  {:.05f}  {:e}".format(mjdtime.value*1e9, i['AZ'], i['EL'], invrange))

if len(args) > 0:
    filename=args[0]
    f = open(filename, 'w')
    for line in pointingdata:
        f.write(line + '\n')
    f.close()
else:
    for line in pointingdata:
        print line

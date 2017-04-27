import ConfigParser

config = ConfigParser.ConfigParser()

config.read("./input/parameter_all.ini")

# read values from a section--section is defined through the title in []
nside = config.getint('parameters','nside')
ifo1 = config.get('parameters','ifo1')
ifo2 = config.get('parameters','ifo2')
GPStime_start = config.getint('parameters','GPStime_start')
GPStime_end = config.getint('parameters','GPStime_end')
segDuration = config.getint('parameters','segDuration')
fLow = config.getint('parameters','fLow')
fHigh = config.getint('parameters','fHigh')
deltaF = config.getfloat('parameters','deltaF')

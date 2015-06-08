"""
This will help manage multiple screens/plots without
overgarbagin the main code
"""
import numpy as np
import xrt.backends.raycing.screens as rsc
#import xrt.backends.raycing.run as rr
import xrt.backends.raycing as raycing
import xrt.plotter as xrtp

# Insert screens in build_beamline() function
# Usage: beamLine.myScreens = screens.createScreens()
# p_range = min, max
# N = number of screens spanned on that distance

def createScreens(beamLine, ranges, howmany):
    positions = np.linspace(ranges[0], ranges[1], howmany)
    screens = {}
    for it in positions:
        screen = rsc.Screen(beamLine,'Screen_at_'+str(it), (0,it,0))
        screens[it] = screen
        print 'Screen at ' + str(it)
    # Save ranges as beamline property, no need to keep globally
    beamLine.s_positions = positions
    beamLine.screens  = screens

# Next screens need to be exposed in run_process and
# dictionary must be prepared
def exposeScreens(beamLine, beamToBeSeen):
    partDict = {}
    for it in beamLine.s_positions:
        partDict['Screen_at_'+str(int(it))] =\
                beamLine.screens[it].expose(beamToBeSeen)
    return partDict

# Prepare plots for raycing
def createPlots(beamLine):
    xLimits = [-beamLine.Dout/1.9, beamLine.Dout/1.9]
    zLimits = xLimits
    cLimits = [0, beamLine.nRefl]
    plots = []
    # Plots after the exit
    for it in beamLine.s_positions:
#        print 'current range: ' + str(it)
#        print 'Screen_at_'+str(it)
        name = 'Screen_at_'+str(int(it))
        plot = xrtp.XYCPlot(name, (1,),
            xaxis=xrtp.XYCAxis(r'$x$', 'mm', data=raycing.get_x,\
                    bins=256, ppb=2, limits=xLimits),
            yaxis=xrtp.XYCAxis(r'$z$', 'mm', data=raycing.get_z,\
                    bins=256, ppb=2, limits=zLimits),
            caxis=xrtp.XYCAxis('Reflections', 'number',\
                    data=raycing.get_reflection_number,\
                    bins=256, ppb=2, limits=cLimits),
            beamState=name, title=name,\
                    aspect='auto')
        plot.baseName = name
        plot.saveName = 'png/' + plot.baseName + '.png'
#        plot.persistentName = 'pickle/' + plot.baseName + '.pickle'
        plots.append(plot)
#    print 'leaving createPlots()'

    return plots

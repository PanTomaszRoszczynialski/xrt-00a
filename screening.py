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
    isUsed  = {}
    for it in positions:
        screen = rsc.Screen(beamLine,'Screen_at_'+str(int(it)), (0,it,0))
        screens[it] = screen
        isUsed[it]  = False
        print 'Screen at ' + str(it)
    # Save ranges as beamline property, no need to keep globally
    beamLine.s_positions = positions
    beamLine.isUsed = isUsed
    beamLine.screens  = screens

# Special function specyfying whicch screens will be used
def setUsed(beamLine, _range):
    for it in beamLine.s_positions:
        if it >= _range[0] and it < _range[1]:
            beamLine.isUsed[it] = True

# Next: screens need to be exposed in run_process and
# dictionary must be prepared
# FIXME: run_process is apparently being run inside run_ray_traycing,
# so createScreens doesn't know yet which screens are being used
cunter = 0
def exposeScreens(beamLine, beamToBeSeen, _range):

    # Debug utilities
    global cunter
    cunter += 1
    print 'It\'s being run: ' + str(cunter) + ' time!'

    partDict = {}
    for it in beamLine.s_positions:
        if it >= _range[0] and it < _range[1]:

            partDict[beamLine.screens[it].name] =\
                    beamLine.screens[it].expose(beamToBeSeen)
            beamLine.isUsed[it] = True

    return partDict

# Prepare plots for raycing
def createPlots(beamLine):
    xLimits = [-beamLine.Dout/1.9, beamLine.Dout/1.9]
    zLimits = xLimits
    cLimits = [0, beamLine.nRefl]
    plots = []
    # Plots after the exit
    for it in beamLine.s_positions:
        if beamLine.isUsed[it]:

            name = 'Screen_at_'+str(int(it))
            print 'Creating plot' + name

            # Plot Construction
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
            # END Plot Construction

            plot.baseName = name
            plot.saveName = 'png/' + plot.baseName + '.png'
#            plot.persistentName = 'pickle/' + plot.baseName + '.pickle'
            plots.append(plot)

    return plots

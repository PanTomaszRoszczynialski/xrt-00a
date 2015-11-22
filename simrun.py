
# Main file for polycapillary optics xrt simulations
# FIXME - for now this is only to come up
# with some kind of interface allowing control
# over many parameters in a simple way

# We need to stay within xrt0.99 paradigm:
# build beamline, override run_process
# and run_ray_tracing

import xrt.backends.raycing as raycing
import xrt.backends.raycing.run as rr
import multiprocessing as mp

# Lower expectations for home computers
# TODO put this in settings
if mp.cpu_count() <= 2:
    repeats = 8
    processes = 1
    print "Running on a slow machine"
else:
    repeats = 100
    processes = 8
    print "Running on a fast machine"

def build_beamline():
    # [0] - Instansiate xrt.BeamLine
    # [1] - Create a source of light
    # [2] - Create a lens object
    # [3] - Create an object 
    # [4] - Create xrt.Screens 

    # [0] - Straighforward, not much to set
    beamLine = raycing.BeamLine(height=0)

    # [1] - Source is abstract with only y-position
    # defined: x and z are defined for each capillary
    # separately in the run_process when shine() is called

    # [2] - This should be a oneliner
    # lens = poly.lens_a

    # [3] - This should be a oneliner as well
    # thing = transmissive.element

    # [4] - This might need to be distributed throughout
    # the whole function?

def run_process(beamLine, shineOnly1stSource=False):
    # [0] - Propagate through the Lens
    # [1] - Propagate through objects
    # [2] - Expose xrt.Screens
    # return out_dict

# This is necessary
rr.run_process = run_process

def main():
    beamLine = build_beamline()

    # Create xrt.Plots
    # TODO - this is tricky and it might be prefered to
    # get rid of those plots and acquire data independently
    plots = []

    xrtr.run_ray_tracing(plots,
                         beamLine = beamLine,
                         repeats = repeats,
                         processes = processes)

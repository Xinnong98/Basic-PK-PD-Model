#
# Magnolia Python script created on 2020-07-23T01:30:34.311
#

# Note: Magnolia doesn't have a capability to do ribbon plots
# yet, so we'll visualize the distribution of trajectoried by 
# plotting the 5, 25, 50, 75 and 95th percentile lines.

# import any needed python packages
import math
import random

# Set the seed for repeatability
random.seed(123456)

# Note: python 3.8 and above has a standard stats library
# that can be used to compute percentile for a collection of data.
# Currently, Magnolia uses Jython, which support Python 2.7 libraries.
# So we'll define a little function here to return the requested 
# percentile for an array of data points
def pctls(data, pct):
    
    # Check validity of range parameter
    if(pct < 1 or pct > 99):
        return math.float("nan")
    
    # Convert the percentile into an index into the data array,
    # which should already be sorted
    idx = int(round(((pct/100.0)*len(data))))
    
    if(idx >= 0 and idx < len(data)):
        return data[idx]
    
    # Error
    return math.float("nan")

# Get a reference to the model object
# Note that the model has to be open in the Magnolia editor
m = models.get("Li_DOX_PBPK_human.csl")

# Number of MC iterations
nits = 500

# Containers for trajectories to be saved
tvals = []  # time point values
cvals1 = []  # RDV plasma concentrations values  plasma
cvals2 = []  # A plasma concentration values   intestine
cvals3 = []  # N plasma concentration values   kidney

# add the outputs of interest to the PREPARE list
# (list of outputs to be logged at each time points)
m.prepare("t")
m.prepare("C_pl")
m.prepare("C_in")
m.prepare("C_li")

# Main MC loop
# Just gather the trajectories here...
# post-processing to calculate percentiles
# and create plots will be done below
for it in range(1, nits):
    
    # Vary the model parameters using normal dists
    # Just use a 20% CV for this example
    cv = 0.2
    
    # TODO: use more appropriate distributions above and add distributions for other parameters
    
    # Physiological values
    #m.vblood     = math.log(random.lognormvariate(5.6, cv*5.6))
    # m.PS_ad= math.log(random.lognormvariate(3.79, cv*3.79))
   

    # Compound-specific properties
    #m.cl_li   = math.log(random.lognormvariate(585, cv*585))  cl_li is huge and the result would be out of range if included so I didn't incorporate it
    m.ps_pa   = math.log(random.lognormvariate(9.546, cv*9.546))
    m.ps_sp   = math.log(random.lognormvariate(3.865, cv*3.865))
    m.ps_ad   = math.log(random.lognormvariate(3.295, cv*3.295))
    m.kplu   = math.log(random.lognormvariate(102, cv*102))
    m.kphrt   = math.log(random.lognormvariate(33.2, cv*33.2))
    m.kpbrn   = math.log(random.lognormvariate(0.245, cv*0.245))
    m.kpmu   = math.log(random.lognormvariate(15, cv*15))
    m.kpkd   = math.log(random.lognormvariate(90.9, cv*90.9))
    m.kpin   = math.log(random.lognormvariate(36.9, cv*36.9))
    m.kpli   = math.log(random.lognormvariate(168, cv*168))
    m.kpre   = math.log(random.lognormvariate(101, cv*101))
    m.f_pa_u   = math.log(random.lognormvariate(0.00399, cv*0.00399))
    m.f_sp_u   = math.log(random.lognormvariate(0.00122, cv*0.00122))
    m.f_ad_u   = math.log(random.lognormvariate(0.0444, cv*0.0444))
    
    
    # Run the simulation
    m.run()
    
    # Collect the outputs of interest
    tvals.append(m.history("t"))
    cvals1.append(m.history("c_pl"))
    cvals2.append(m.history("c_in"))
    cvals3.append(m.history("c_li"))


# Visualize the distribution of concentration trajectories using
# confidence intervals of predictions at each time point

# Iterate over time points
tt = tvals[0];

# Container for the percentile curves which will be plotted later
C_pll05 = []
C_pll25 = []
C_pll50 = []
C_pll75 = []
C_pll95 = []

C_inl05 = []
C_inl25 = []
C_inl50 = []
C_inl75 = []
C_inl95 = []

C_lil05 = []
C_lil25 = []
C_lil50 = []
C_lil75 = []
C_lil95 = []

for i in range(0, len(tt)):
    
    # Create the list of conc values at this time point
    # using Python list comprehension
    cc = [d[i] for d in cvals1]
    cc2 = [d[i] for d in cvals2]
    cc3 = [d[i] for d in cvals3]
    
    # Sort the list of conc values to make it easy to find percentile
    cc.sort()
    cc2.sort()
    cc3.sort()
    
    # Now calculate percentiles
    C_pll05.append(pctls(cc, 5))
    C_pll25.append(pctls(cc, 25))
    C_pll50.append(pctls(cc, 50))
    C_pll75.append(pctls(cc, 75))
    C_pll95.append(pctls(cc, 95))
    
    C_inl05.append(pctls(cc2, 5))
    C_inl25.append(pctls(cc2, 25))
    C_inl50.append(pctls(cc2, 50))
    C_inl75.append(pctls(cc2, 75))
    C_inl95.append(pctls(cc2, 95))
    
    C_lil05.append(pctls(cc3, 5))
    C_lil25.append(pctls(cc3, 25))
    C_lil50.append(pctls(cc3, 50))
    C_lil75.append(pctls(cc3, 75))
    C_lil95.append(pctls(cc3, 95))
    
# Finally, plot the percentiles.  Use different colors
# for the 5/95, 25/75 and 50 percentiles
h = plot.line(tt, C_pll05, "b-", "5th pctl")
plot.append(h, tt, C_pll25, "r-", "25th pctl")
plot.append(h, tt, C_pll50, "k-", "50th pctl")
plot.append(h, tt, C_pll75, "r-", "75th pctl")
plot.append(h, tt, C_pll95, "b-", "95th pctl")
plot.title(h, "Plasma Concentration Variability")
plot.xlabel(h, "Time (h)")
plot.ylabel(h, "Plasma Concentration  (ug/mL)")

# for the 5/95, 25/75 and 50 percentiles
h = plot.line(tt, C_inl05, "b-", "5th pctl")
plot.append(h, tt, C_inl25, "r-", "25th pctl")
plot.append(h, tt, C_inl50, "k-", "50th pctl")
plot.append(h, tt, C_inl75, "r-", "75th pctl")
plot.append(h, tt, C_inl95, "b-", "95th pctl")
plot.title(h, "Intestine Concentration Variability")
plot.xlabel(h, "Time (h)")
plot.ylabel(h, "Intestine Concentration  (ug/mL)")

# for the 5/95, 25/75 and 50 percentiles
h = plot.line(tt,  C_lil05, "b-", "5th pctl")
plot.append(h, tt, C_lil25, "r-", "25th pctl")
plot.append(h, tt, C_lil50, "k-", "50th pctl")
plot.append(h, tt, C_lil75, "r-", "75th pctl")
plot.append(h, tt, C_lil95, "b-", "95th pctl")
plot.title(h, "Liver Concentration Variability")
plot.xlabel(h, "Time (h)")
plot.ylabel(h, "Kidney Concentration  (ug/mL)")

# Create an identical plotm, but log scale
#h = plot.line(tt, pctl05, "b-", "5th pctl")
#plot.append(h, tt, pctl25, "r-", "25th pctl")
#plot.append(h, tt, pctl50, "k-", "50th pctl")
#plot.append(h, tt, pctl75, "r-", "75th pctl")
#plot.append(h, tt, pctl95, "b-", "95th pctl")
#plot.title(h, "RDV Plasma Concentration Variability")
#plot.xlabel(h, "Time (h)")
#plot.ylabel(h, "Concentration  (mg/L)")
#plot.logy(h, True)
# Adjust the y axis limits
#plot.ymin(h, 1.5e-3)
#plot.ymax(h, 4.0)






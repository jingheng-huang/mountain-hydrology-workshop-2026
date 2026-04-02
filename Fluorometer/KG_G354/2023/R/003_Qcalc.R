# - reading FL30 raw data
# - transforming data into readable format with signal in mV
require(readxl)
require(lubridate)

# Arguments
path = "/Users/pohle/Dropbox/Central_Asia/Fluorometer/2023/G354/data/"

Standard = 50        # g/l - The concentration of the color tracer (Standard solution)
TracerVolume = 0.1  # l - How many liters! were injected into the river?

# Input file (Excel)
filename_excel = 'f1664_G354.xlsx'


sheet_name = '20230713_0950'  
TracerVolume = 0.02 

sheet_name = '20230713_1020'  
TracerVolume = 0.05

sheet_name = '20230713_1100'  
TracerVolume = 0.1

sheet_name = '20230713_1110'  
TracerVolume = 0.2









# ------------ NO NEED TO CHANGE ANYTHING HEREAFTER --------------- #
gc() # This cleans the memory - against JAVA Problems ...
setwd(path)
source('../R/FUN_install_packages.R')
options(java.parameters = "-Xmx4000m")  # allow for 2GB RAM - for JAVA
M = Standard * TracerVolume * 1000 # to mg
# Source all functions for cleaning and preparation
source('../../../R/FUN_curve-processing.R')

# ------- READ DATA ------- #
# d_ = read.xlsx(sheetName = sheet_name,  # this one needs to be changed to read in the different sheets,
#                file = filename_excel, 
#                header = F)
d_ = read_xlsx_(filename_excel, sheet_name)

# ------- PROCESS DATA ---- #
# ---- remove NAs 
d_ = rm_nas(d_)

# ---- Make date-time
tmp_date = make_date_time(d_)

# ---- Combine data
d = data.frame(date=tmp_date, Tracer1=d_$X4)

# ---- Noise removal
d = remove_noise(d)

# ---- Calibration
calibrate(d,1,200)

# ---- Relationship mV to ppb
d = relate_mV_ppb(d)


# ---- Discharge Calculation --------- #
# Calculate discharge with the global method
# ppb = 1 µg/L
# Q = M/A
# - M [mg] = 1000 µg
# - A [µg*s / L]
# ->
# L*1000/s = m3/s
# ---------------------------------------- #


# ---- Curve 1 
calc_discharge(d)


# ---- Curve 2 
# If second (overlapping) breakthrough curve exists:
print("Is there a second breakthrough curve (with overlap)?")

# Was the same amount of tracer used? --> Change if yes!
# ----- Change the following two Arguments if need be ---------------------------- #
Standard = Standard          # g/l - The concentration of the color tracer (Standard solution)
TracerVolume = TracerVolume  # l - How many liters! were injected into the river?
# -------------------------------------------------------------------------------- #
M = Standard * TracerVolume * 1000 

# calc_discharge_overlap(d, dss)



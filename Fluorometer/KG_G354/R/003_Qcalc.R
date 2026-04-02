require(readxl)
require(lubridate)
require(MESS)

# - reading FL30 raw data
# - transforming data into readable format with signal in mV

# Arguments
path = "/Users/pohle/Dropbox/Fribourg/2026/workshops/CA-hydroworkshop/materials/Fluorometer/KG_G354"

Standard = 50        # g/l - The concentration of the color tracer (Standard solution)
TracerVolume = 0.1  # l - How many liters! were injected into the river?

# ---- 2021 -----
# Input file (Excel)
filename_excel = '2021 - Measurements_G354.xlsx'

# ---- 
sheet_name = '20210928_1052'    # Use this for all measurements
TracerVolume = .02
sheet_name = '20210928_1220_cleaned'    # Use this for all measurements
TracerVolume = .025
sheet_name = '20210928_1310'    # Use this for all measurements
TracerVolume = .04
sheet_name = '20210928_1409'    # Use this for all measurements
TracerVolume = .06
sheet_name = '20210928_1510_cleaned'    # Use this for all measurements
TracerVolume = .12
sheet_name = '20210928_1600'    # Use this for all measurements
TracerVolume = .12

# ---- 2022-----
# Input file (Excel)
filename_excel = '2022 - Measurements_G354.xlsx'

# ---- 
sheet_name = '20220817_1100'    # Use this for all measurements
TracerVolume = .02
sheet_name = '20220817_1115'    # Use this for all measurements
TracerVolume = .03
sheet_name = '20220817_1440'    # Use this for all measurements
TracerVolume = .2
sheet_name = '20220817_1500'    # Use this for all measurements
TracerVolume = .4
sheet_name = '20220817_1700'    # Use this for all measurements
TracerVolume = .3
sheet_name = '20220817_1900'    # Use this for all measurements
TracerVolume = .3
sheet_name = '20220818_0900'    # Use this for all measurements
TracerVolume = .1


# ---- 2023 -----
# Input file (Excel)
filename_excel = '2023 - Measurements_G354.xlsx'

# ---- 
sheet_name = '20230713_0950'    # Use this for all measurements
TracerVolume = .02
sheet_name = '20230713_1020_repaired'    # Use this for all measurements
TracerVolume = .05
sheet_name = '20230713_1100'    # Use this for all measurements
TracerVolume = .1
sheet_name = '20230713_1110'    # Use this for all measurements
TracerVolume = .2

# ---- 2024 -----
# Input file (Excel)
filename_excel = '2024 - Measurements_G354.xlsx'

# ---- 
sheet_name = '20240820_0900'    # Use this for all measurements
TracerVolume = .050
sheet_name = '20240820_0933'    # Use this for all measurements
TracerVolume = .0265
sheet_name = '20240820_1000'    # Use this for all measurements
TracerVolume = .025
sheet_name = '20240820_1030'    # Use this for all measurements
TracerVolume = .041
sheet_name = '20240820_1103'    # Use this for all measurements
TracerVolume = .031
sheet_name = '20240820_1155_cleaned'
TracerVolume = .041
sheet_name = '20240820_1215_cleaned'
TracerVolume = .052

sheet_name = '20240820_1320'
TracerVolume = 0.0565
sheet_name = '20240820_1418_cleaned'
TracerVolume = 0.1
sheet_name = '20240820_1510_cleaned'
TracerVolume = 0.2
sheet_name = '20240820_1600_repair'
TracerVolume = 0.2 # ??????????????????????????
sheet_name = '20240820_1638'
TracerVolume = 0.2 # ??????????????????????????
sheet_name = '20240821_0615'
TracerVolume = 0.150 

sheet_name = '20240821_0630'
TracerVolume = 0.07 
sheet_name = '20240821_0700'
TracerVolume = 0.07 
sheet_name = '20240821_0730'
TracerVolume = 0.07 



# ------------ NO NEED TO CHANGE ANYTHING HEREAFTER --------------- #
# gc() # This cleans the memory - against JAVA Problems ...
setwd(path)
# source('../R/FUN_install_packages.R')
M = Standard * TracerVolume * 1000 # to mg
# Source all functions for cleaning and preparation
source('../0_R-FUN/FUN_curve-processing.R')

# ------- READ DATA ------- #
# d_ = read.xlsx(sheetName = sheet_name,  # this one needs to be changed to read in the different sheets,
#                file = filename_excel, 
#                header = F)
library(lubridate)
d_ = read_xlsx_(filename = filename_excel,sheet = sheet_name)
# ---- remove NAs 
d_ = rm_nas(d_)
plot(d_$X1, d_$X4, type='l')
require(zoo)
ttt = na.approx(d_$X4)
d_$X4 = ttt
lines(d_$X1, ttt, col='red')
d_ = na.omit(d_)
# ---- Make date-time
tmp_date = make_date_time(d_)

# ---- Combine data
d = data.frame(date=tmp_date, Tracer1=d_$X4)

# ---- Noise removal
d = remove_noise(d)

# ---- Calibration
calibrate(d, 1, 600)

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



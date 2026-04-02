# - reading HOBO data
# - calculating water level (removing atmospheric pressure from water pressure)
# - plot corrected water level
# -------------------------------------
# logger time GMT+02:00
# start at Aug 20, 04:00 GMT+02:00 == 08:00 local time GMT+06:00
# (Kyrgyzstan does not observe daylight saving time)
# --->
# 2023-12-13: Removed the need for time shift. 
# Instead EXPORT AS GMT+6 or whatever timezone you are in, in the HOBOware
time_difference_hours = 0

# Arguments
path = "/Users/pohle/Dropbox/Central_Asia/Fluorometer/2023/G354/data/"
# For Windows it will be something like this:
# path = "C:/summer_school_Suyek/data/"

# Input files
filename_atmo_pressure = 'atmo-21411022.csv'
filename_water_pressure_clean = 'water-clean-21411021.csv'
filename_water_pressure_original = 'water-21411021.csv'   # Excel is messing up the date! Take origninal dates from original file

# Output files
filename_output = 'HOBO_water-corrected.csv'

# Constant:
kPa_to_cm <- 1/0.0980665  # to transform kPa to cm water level




# ------------ NO NEED TO CHANGE ANYTHING HEREAFTER --------------- #
setwd(path)
source('../R/FUN_install_packages.R')   # this checks if packages are installed

# This reads all the data
dat_atmo = read.csv(filename_atmo_pressure, skip=2, header=F, stringsAsFactors = F)
dat_wl_clean = read.csv(filename_water_pressure_clean, skip=2, header=F, stringsAsFactors = F)
dat_wl_orig = read.csv(filename_water_pressure_original, skip=2, header=F, stringsAsFactors = F)


date_at = as.POSIXct(dat_atmo$V2, format='%m/%d/%Y %I:%M:%S %p', tz='GMT')
# date_wl = as.POSIXct(dat_wl_orig$V2, format='%m/%d/%Y %I:%M:%S %p', tz='GMT')
# Some mess up.. with the original hobo file. A corrupt end of file apparently destroys the dates.
# It seems like this and the atmo logger started at the same time. So I use the atmo logger time to overwrite the water level ogger

# However, this requires ADDING A ----->             ,      < ---------- after %Y !!!!!!!
# However, this requires ADDING A ----->             ,      < ---------- after %Y !!!!!!!
date_wl = as.POSIXct(dat_wl_orig$V2, format='%m/%d/%Y %I:%M:%S %p', tz='GMT')




wl = dat_wl_clean$V3   # water + atmospheric pressure
at = dat_atmo$V3       # atmospheric pressure
ta = dat_atmo$V4       # atmospheric temperature
tw = dat_wl_clean$V4   # water temperature

# merge and align
df_wl = data.frame(date= date_wl, level=wl, water_temperature=tw)
df_at = data.frame(date= date_at, atmo=at, atmo_temperature=ta)
wl_at = merge(df_wl, df_at, by='date')
plot(wl_at$date, wl_at$wl_clean, type='l')
range(wl_at$date)
# na_ind = which(wl_at$level>75)
# wl_at[na_ind,] = NA
ind_2023 = which(format(wl_at$date, '%Y') == '2023'  & format(wl_at$date, '%m') == '07')
wl_2023 = wl_at[ind_2023,]
wl_2023 = na.omit(wl_2023)
plot(wl_2023$date, wl_2023$level, type='l')
par(new=T)
plot(wl_2023$date, wl_2023$water_temperature, type='l', col='turquoise', ann=F, axes=F, ylim=c(-.1, .1))
abline(h=0, col='red', lty=3)
axis(4, col='turquoise')
#
plot(wl_2023$date, wl_2023$wl_clean, type='l')
lines(wl_2023$date, wl_2023$atmo - 56, type='l', col='blue')

wl_at$wl_clean = wl_at$level-wl_at$atmo
plot(wl_at$date,wl_at$wl_clean, type='l', xlab = 'Time', ylab='Water level [dm]')
plot(wl_at$date[1:30],wl_at$wl_clean[1:30], type='l', xlab = 'Time', ylab='Water level [dm]')

write.csv(wl_at,file = filename_output,quote = F,row.names = F )



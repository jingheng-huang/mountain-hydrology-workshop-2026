# rating curve between discharge measurements and HOBO water level


# Arguments
path = "/Users/pohle/Dropbox/Central_Asia/HOBO/summer_school_Suyek/data/"

# ------ INPUT FILES -------
# HOBO water level logger file (corrected for atmospheric pressure)
filename_water_pressure_clean = 'HOBO_water-corrected_suyek.csv'
# FL30 derived discharges 
filename_fl30 = "Overview measurements.csv"

# ------ OUTPUT FILES -------
# Discharge - calibrated water level
filename_discharge_out = 'Q_rating-curve_suyek.csv'


# ------------ NO NEED TO CHANGE ANYTHING HEREAFTER --------------- #
setwd(path)
source('../R/FUN_install_packages.R')
options(java.parameters = "-Xmx4000m")  # allow for 2GB RAM - for JAVA
# Source all functions for cleaning and preparation
source('../R/FUN_curve-processing.R')


# This reads all the data
# HOBO
dat_wl = read.csv(filename_water_pressure_clean, skip=0, header=T, stringsAsFactors = F)
dat_wl$date = as.POSIXct(dat_wl$date, format='%Y-%m-%d %H:%M:%S', tz='GMT')
df_wl = data.frame(date=dat_wl$date, WL=dat_wl$wl_clean)


as.POSIXct(as.numeric(datetime), origin=as.POSIXct("1970-01-01", tz="Europe/Paris"),
           tz="Europe/Paris")

# FL30 discharge
dat_q = read.csv(filename_fl30)
date_q = as.POSIXct(paste(dat_q$date, dat_q$time), format='%d.%m.%y %H:%M', tz='GMT')
df_q = data.frame(date = date_q, Q=dat_q[,3])

# merge and align
# - xts to apply hourly
q_xts = xts(df_q$Q, order.by = df_q$date)
q_xts = na.omit(q_xts)

wl_xts = xts(df_wl$WL, order.by = df_wl$date)
# Interpolate the WL to minutes 
t_min_q = min(time(q_xts))
t_max_q = max(time(q_xts))

# WL min max time
match_ind = which(df_wl$date >= t_min_q & df_wl$date <= t_max_q)
df_sub = df_wl[match_ind,]

# Check water level at day of FL30 measurements

date_s = as.POSIXct('2023-07-06 00:00:00')
date_e = as.POSIXct('2023-07-07 00:00:00')
# t_ind_start_ = which.min(abs(df_wl$date - t_min_q))
# t_ind_end_ = which.min(abs(df_wl$date - t_max_q))


# date_s_ = df_wl$date[t_ind_start_]
# date_e_ = df_wl$date[t_ind_end_]
crop_time_s = 0
crop_time_e = 0
# xts 
plot(wl_xts[(t_ind_start_-crop_time_s -crop_time_s ):(t_ind_end_ -crop_time_e)])



#
x_time = df_wl$date[(t_ind_start_ - crop_time_s):(t_ind_end_ )]
x_time_limit = range(x_time)

plot(x_time, df_wl$WL[(t_ind_start_ - crop_time_s):(t_ind_end_ )],
     type='l', xlim=x_time_limit,
     axes=T)

plot(x_time, df_wl$WL[(t_ind_start_ - crop_time_s):(t_ind_end_ )],
     type='l', xlim=x_time_limit,
     axes=F)
axis(1)
par(new=T)
plot(df_q$date, df_q$Q, xlim=x_time_limit, axes=F, ann=F)
# axis.Date(1,df_q$date)
axis(1)
axis(4)



t_mins = seq(t_min_q, t_max_q, by='mins')  
df_all = data.frame(date=t_mins, val=NA)
df_wl_min = merge(df_all,df_sub, on='date', all.x = T)
df_wl_min = df_wl_min[,names(df_wl_min)!='val']
xts_wl_min = xts(df_wl_min$WL, order.by = df_wl_min$date)

xts_wl_min_approx = na.approx(xts_wl_min)
qwl_xts = merge(xts_wl_min_approx,q_xts,all = F)
# plot(as.data.frame(qwl_xts), xlab = 'Water level [cm]', ylab='Discharge [m3/s]')



# -------- RATING CURVE ----------
df = as.data.frame(qwl_xts)
lev = c(0,df$xts_wl_min_approx)
mes = c(0,df$q_xts)
# plot(lev,mes)


require(pracma)
p <- polyfit(lev,mes,n = 2)
p2 <- lm(mes~lev)
xnew= seq(0,38)
yf <- polyval(p, xnew)
# yf <- predict(p,newdata = data.frame(xnew))
plot(lev,mes, xlab = 'Water level [cm]', ylab='Discharge [m3/s]')
lines(xnew, yf, col="red")

# making the series into Q
df_wl$Q <- polyval(p, as.numeric(df_wl$WL))
df_wl$Q[df_wl$Q<0] <- NA
plot(df_wl$date, df_wl$Q, type='l', ylab='Discharge [m3/s]', xlab='Time')

# save output:
write.csv(df_wl, file = filename_discharge_out, row.names = F)
# ----------------------------------------------------------------
# Spring/Summer 2023
month = format(df_wl$date, '%m')
year = format(df_wl$date, '%Y')

ind_2023 = which(year == '2023' & as.numeric(month) >=  4)
df_wl_2023 = df_wl[ind_2023,]
# remove NA dates
ind_na = which(!is.na(df_wl_2023$Q))
df_wl_2023_clean = df_wl_2023[ind_na,]

plot(df_wl_2023_clean$date, df_wl_2023_clean$Q, type='l', ylab='Discharge [m3/s]', xlab='Time')


# Example calculations

# Sum of discharge in June -------------------------------------
start = as.POSIXct('2023-06-01', tz='GMT')
end = as.POSIXct('2023-06-30', tz='GMT')
ind_jun_2023 = which(df_wl_2023_clean$date >= start & df_wl_2023_clean$date <= end)
df_jun_2023 = df_wl_2023_clean[ind_jun_2023,]
xts_jun_2023 = xts(df_jun_2023$Q, order.by = df_jun_2023$date)
plot(xts_jun_2023, ylab='Discharge [m3/s]', xlab='Time')

# from m3/s to m3/june: *60 (seconds) *60 (minutes) *24 (hours) *30 (days)
Vwater = mean(xts_jun_2023)*60*60*24*30  # m3
print(paste0('Discharge volume in June: ', round(Vwater*1000,2) ,' liters'))
AreaGlacier = 1000000 # around 1km2
print(paste0('Specific discharge in June: ', round(Vwater/AreaGlacier*1000,2) ,' mm'))


# Sum of discharge in August 2022 -------------------------------------
start = as.POSIXct('2022-08-01', tz='GMT')
end = as.POSIXct('2022-08-31', tz='GMT')
ind_aug_2022 = which(df_wl$date >= start & df_wl$date <= end)
df_aug_2022 = df_wl[ind_aug_2022,]
xts_aug_2022 = xts(df_aug_2022$Q, order.by = df_aug_2022$date)
plot(xts_aug_2022, ylab='Discharge [m3/s]', xlab='Time')

# from m3/s to m3/june: *60 (seconds) *60 (minutes) *24 (hours) *30 (days)
Vwater = mean(xts_aug_2022, na.rm=T)*60*60*24*30  # m3
print(paste0('Discharge volume in Aug: ', round(Vwater*1000,2) ,' liters'))
AreaGlacier = 1000000 # around 1km2
print(paste0('Specific discharge in Aug: ', round(Vwater/AreaGlacier*1000,2) ,' mm'))


# rating curve between discharge measurements and HOBO water level


# Arguments
path = "/Users/pohle/Dropbox/Central_Asia/Fluorometer/2023/G354/data/"

# ------ INPUT FILES -------
# HOBO water level logger file (corrected for atmospheric pressure)
filename_water_pressure_clean = 'HOBO_water-corrected.csv'
# FL30 derived discharges 
filename_fl30 = "Overview measurements 2023-12-15.csv"

# ------ OUTPUT FILES -------
# Discharge - calibrated water level
filename_discharge_out = 'Q_rating-curve_g354.csv'


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
val_max_wl = max(df_wl$WL, na.rm=T)
plot(df_wl$date, df_wl$WL, type='l', xlab='Time', ylab='Water level [dm]')
tail(df_wl)

# FL30 discharge
dat_q = read.csv(filename_fl30)
date_q = as.POSIXct(paste(dat_q$date, dat_q$time), format='%d.%m.%y %H:%M', tz='GMT')
df_q = data.frame(date = date_q, Q=dat_q[,3])
df_q = na.omit(df_q)

# merge and align
# - xts to apply hourly
q_xts = xts(df_q$Q, order.by = df_q$date)
q_xts = na.omit(q_xts)

# Interpolate the WL to minutes 
t_min_q = min(time(q_xts))
t_max_q = max(time(q_xts))

# WL min max time
t_min_q = t_min_q - minute(t_min_q)*60
match_ind = which(df_wl$date >= t_min_q & df_wl$date <= t_max_q)
df_sub = df_wl[match_ind,]

t_mins = seq(t_min_q, t_max_q, by='mins')  
df_all = data.frame(date=t_mins, val=NA)
df_wl_min = merge(df_all,df_sub, on='date', all.x = T)
df_wl_min = df_wl_min[,names(df_wl_min)!='val']
xts_wl_min = xts(df_wl_min$WL, order.by = df_wl_min$date)

xts_wl_min_approx = na.approx(xts_wl_min)
qwl_xts = merge(xts_wl_min_approx,q_xts,all = F)
# plot(as.data.frame(qwl_xts), xlab = 'Water level [dm]', ylab='Discharge [m3/s]')
plot(xts_wl_min_approx)


# -------- RATING CURVE ----------
df = as.data.frame(qwl_xts)
lev = c(0,df$xts_wl_min_approx)
mes = c(0,df$q_xts)
# plot(lev,mes)


require(pracma)
p <- polyfit(lev,mes,n = 5)
p2 <- lm(mes~lev)
xnew= seq(0,val_max_wl*1.1,length.out=50)
yf <- polyval(p, xnew)

# The theoretical values for the quadratic/polynomial model
# par(mfrow=c(2,1))
plot(xnew, yf, xlab = 'Water level [dm]', ylab='Discharge [m3/s]')
# yf <- predict(p,newdata = data.frame(xnew))
plot(lev,mes, xlab = 'Water level [dm]', ylab='Discharge [m3/s]')
lines(xnew, yf, col="red")


# making the series into Q
df_wl$Q <- polyval(p, as.numeric(df_wl$WL))
df_wl$Q[df_wl$Q<0] <- NA
plot(df_wl$date, df_wl$Q, type='l', ylab='Discharge [m3/s]', xlab='Time')

par(mar=c(4,4,1,4))
plot(df_wl$date[match_ind], df_wl$Q[match_ind], type='l', ylab='Discharge [m3/s]', xlab='Time',
     ylim=range(c(df_wl$Q[match_ind], as.vector(q_xts) ), na.rm = T))
points(time(q_xts), q_xts, col='red')
par(new=T)
plot(df_wl$date[match_ind], df_wl$WL[match_ind],type='l', lty=2, col='blue', axes=F, ann=F)
axis(4, col='blue')
mtext(text='Water level [dm]',side = 4,line = 2.6, col='blue')
legend('topright', legend=c('Q(calc)', "Q(measured)", "Water level"), 
       col=c(1,2,'blue'), lty=c(1,-1,2), pch=c(-1,1,-1))

# save output:
write.csv(df_wl, file = filename_discharge_out, row.names = F)


xts_ql_out = xts(df_wl$Q, order.by = df_wl$date)
plot(xts_ql_out)
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
AreaGlacier = 1.58*10^6 # around 1km2


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
print(paste0('Specific discharge in June: ', round(Vwater/AreaGlacier*1000,2) ,' mm'))

# Sum of discharge in July 2023 -------------------------------------
start = as.POSIXct('2023-07-01', tz='GMT')
end = as.POSIXct('2023-07-31', tz='GMT')
ind_jul_2022 = which(df_wl$date >= start & df_wl$date <= end)
df_jul_2022 = df_wl[ind_jul_2022,]
xts_jul_2022 = xts(df_jul_2022$Q, order.by = df_jul_2022$date)
plot(xts_jul_2022, ylab='Discharge [m3/s]', xlab='Time')

# from m3/s to m3/june: *60 (seconds) *60 (minutes) *24 (hours) *30 (days)
Vwater = mean(xts_jul_2022, na.rm=T)*60*60*24*30  # m3
print(paste0('Discharge volume in jul: ', round(Vwater*1000,2) ,' liters'))
print(paste0('Specific discharge in jul: ', round(Vwater/AreaGlacier*1000,2) ,' mm'))


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
print(paste0('Specific discharge in Aug: ', round(Vwater/AreaGlacier*1000,2) ,' mm'))


# Sum of discharge in 2023 -------------------------------------
start = as.POSIXct('2023-06-01', tz='GMT')
end = as.POSIXct('2023-08-31', tz='GMT')
ind_aug_2022 = which(df_wl$date >= start & df_wl$date <= end)
df_aug_2022 = df_wl[ind_aug_2022,]
xts_aug_2022 = xts(df_aug_2022$Q, order.by = df_aug_2022$date)
plot(xts_aug_2022, ylab='Discharge [m3/s]', xlab='Time')
test = apply.daily(xts_aug_2022, mean, na.rm=T)
plot(test, main='',ylab='Discharge [m3/s]', xlab='Time')
# from m3/s to m3/june: *60 (seconds) *60 (minutes) *24 (hours) *30 (days)
Vwater = mean(xts_aug_2022, na.rm=T)*60*60*24*30  # m3
print(paste0('Discharge volume in Aug: ', round(Vwater*1000,2) ,' liters'))
print(paste0('Specific discharge in Aug: ', round(Vwater/AreaGlacier*1000,2) ,' mm'))


# Sum of discharge in 2022 -------------------------------------
start = as.POSIXct('2022-07-01', tz='GMT')
end = as.POSIXct('2022-10-15', tz='GMT')
ind_aug_2022 = which(df_wl$date >= start & df_wl$date <= end)
df_aug_2022 = df_wl[ind_aug_2022,]
xts_aug_2022 = xts(df_aug_2022$Q, order.by = df_aug_2022$date)
plot(xts_aug_2022, ylab='Discharge [m3/s]', xlab='Time')
# test = apply.daily(xts_aug_2022, mean, na.rm=T)
# plot(test, ylab='Discharge [m3/s]', xlab='Time', main='')
# from m3/s to m3/june: *60 (seconds) *60 (minutes) *24 (hours) *30 (days)
Vwater = mean(xts_aug_2022, na.rm=T)*60*60*24*30  # m3
print(paste0('Discharge volume in Aug: ', round(Vwater*1000,2) ,' liters'))
print(paste0('Specific discharge in Aug: ', round(Vwater/AreaGlacier*1000,2) ,' mm'))


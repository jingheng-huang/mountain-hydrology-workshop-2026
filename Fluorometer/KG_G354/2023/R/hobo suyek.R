# -----  PART I -------
# - reading data
# - calculating water level
# - plot water level

# load packages:
require('MESS') # for AUC

kPa_to_cm <- 1/0.0980665 #

setwd('~/Dropbox/Fribourg/2021/fieldwork/hobo/golubin/')

dat_atmo = read.csv('10778185_atmo2_clean.csv', skip=2, header=F, stringsAsFactors = F)
dat_wl = read.csv('10778183_water_clean.csv', skip=2, header=F, stringsAsFactors = F)
dat_cd = read.csv('10292196_conduc.csv', skip=2, header=F, stringsAsFactors = F)

# Sys.getlocale (category = "LC_ALL")

date_at = as.POSIXct(dat_atmo$V2, format='%m/%d/%Y %I:%M:%S %p', tz='GMT')
date_wl = as.POSIXct(dat_wl$V2, format='%m/%d/%Y %I:%M:%S %p', tz='GMT')
date_cd = as.POSIXct(dat_cd$V2, format='%m/%d/%Y %I:%M:%S %p', tz='GMT')

# which(!date_at %in% date_wl)
# date_at[16697]
# tail(date_at)
# tail(date_wl)
wl = dat_wl$V3 
at = dat_atmo$V3
tw = dat_wl$V4
wl_clean = wl-at
cd = dat_cd$V3


which(is.na(date_at))
date_at[5220:5229]
require(xts)
wl_xts  = xts(wl_clean, order.by = date_at) *kPa_to_cm
wt_xts  = xts(tw, order.by = date_at)
cd_xts  = xts(cd, order.by = date_cd)
plot(wl_xts)

# -----  PART II -------
# - select data of interest
# - select by date in the format YYYY-MM-DD
# - example below for melting periods in 2020 nad 2021 (May to October)
ind2020 = which(time(wl_xts) > as.POSIXct('2020-05-01') & time(wl_xts) < as.POSIXct('2020-10-10'))
ind2021 = which(time(wl_xts) > as.POSIXct('2021-05-01') & time(wl_xts) < as.POSIXct('2021-10-10'))
# ind2020cd = which(time(cd_xts) > as.POSIXct('2020-05-01') & time(cd_xts) < as.POSIXct('2020-10-10'))
# ind2021cd = which(time(cd_xts) > as.POSIXct('2021-05-01') & time(cd_xts) < as.POSIXct('2021-10-10'))

# plot e.g. the 2020 data
# Water level
plot(wl_xts[ind2020], ylim=c(0,66))
par(new=T)
# and temperature
lines(wt_xts[ind2020], col='blue')
# and conductivity
lines(cd_xts[ind2020cd], col = 'red')


# -------------# -------------# -------------# -------------# -------------# -------------# -------------
# recession model from Jobard and Dzikowski 2016
# (2) Q(t) = Q0 / (1+a*t)^3
# Q(t) - the discharge in m3·s−1 at any time t
# Q0   - the discharge in m3·s−1 at time t0, corresponding to the time 2 h after peak discharge;
# a    - recession parameter in h−1
# t    - time in hours from t0
# (3) V = AUC[Q(t)][t_i : t_{i+1}]  +  Q_{0_i} * AUC[(1/[1+a_i(t-t_{0_i})]^3)][t_{i+1} : inf] -  Q_{0_{i-1}} * AUC[(1/[1 + a{i-1} * (t-t{0_i})]^3)][t_i : inf] 
# Q_{0_i}  - parameters of the recession law fit to the depletion corresponding to the flood on day i–1; 
# a_i      - parameters of the recession law fit to the depletion corresponding to the flood on day i–1; 
# t_{0_i}  - time corresponding to Q_0_i
# t_i      - start time of the extrapolated recession with parameters Q_{0_{i-1}} and a_{i-1}


# -------------# -------------# -------------# -------------# -------------# -------------# -------------

plot(wl_xts[seq(from=16641,to=16691)])
measurements_time = c(as.POSIXct('2021-07-20 08:00:00', tz='GMT'),
                      as.POSIXct('2021-07-18 17:00:00', tz='GMT'),
                      as.POSIXct('2021-07-19 10:00:00', tz='GMT'))#,
                      
measurements_time_num = as.numeric(measurements_time)
mes_val = c(0.45,   # the measurements of the fluorometer test
            1,
            0.74)
 
head(wl_xts)
START = 16641  # this needs to be adjusted once the real data is available
END = 16695
x = as.POSIXct(time(wl_xts[,1][seq(from=START,to=END)]))
y = as.vector(wl_xts[,1][seq(from=START,to=END)])

pdf('q_vs_level_gulobin.pdf',width=4, height = 3)
par(mar=c(4,4,2,4))
plot(x,y, type='l', axes=F,ann=F)
axis(2)
box()
mtext('Water level [cm]',side = 2, 2.6)
axis(1, at =as.numeric(x),labels =  format(x,'%H:%M'),las=2)
par(new=T)

abline(v=measurements_time_num, col='red', lty=2, lwd=2)

par(new=T)
plot(measurements_time_num,mes_val,xlim=range(x) ,ann=F, axes=F, pch=16, cex=1.4, col='red')
axis(4)
mtext(expression('Discharge ['*m^3*s^{-1}*']'),side = 4, 2.6, col='red')

dev.off()

# -------- RATING CURVE ----------
rc_time = which(as.numeric(time(wl_xts))%in%measurements_time)
lev = c(0,as.numeric(wl_xts[rc_time,1]))
mes = c(0,mes_val[c(2,3,1)])
# plot(lev,mes)


require(pracma)
p <- polyfit(lev,mes,n = 2)
p2 <- lm(mes~lev)
xnew= seq(0,38)
yf <- polyval(p, xnew)
# yf <- predict(p,newdata = data.frame(xnew))
plot(lev,mes)
lines(xnew, yf, col="red")

# making the series into Q
yf_all <- polyval(p, as.numeric(wl_xts))
yf_all <- yf_all
yf_all[yf_all<0] <- NA
plot(yf_all,type='l')

wl_xts_adj = xts(yf_all, time(wl_xts))
wl_xts_adj <- wl_xts_adj


plot(wl_xts_adj[ind2021], ylim=c(0,2.5))
x1 = as.POSIXct(time(wl_xts_adj[,1][seq(from=START,to=END)]))
y1 = as.vector(wl_xts_adj[,1][seq(from=START,to=END)])



ind2021[1]
START2=15500
plot(wl_xts_adj[START2:END], ylim=c(0,6))
x1 = as.POSIXct(time(wl_xts_adj[,1][seq(from=START2,to=END)]))
y1 = as.vector(wl_xts_adj[,1][seq(from=START2,to=END)])


pdf('q_golubin_calibrated.pdf',width=4, height = 3)
par(mar=c(4,4,2,4))
plot(x1,y1, type='l', axes=F,ann=F)
axis(2)

dat_seq = seq(as.POSIXct('2021-05-01'), as.POSIXct('2021-09-01'),'days')
axis(1, at = dat_seq,labels = format(dat_seq, '%m-%d'),las=2)
box()
mtext(expression('Discharge ['*m^3*s^{-1}*']'),side = 2, 2.6)
mtext('2021',side = 3, -1.6)
dev.off()


#

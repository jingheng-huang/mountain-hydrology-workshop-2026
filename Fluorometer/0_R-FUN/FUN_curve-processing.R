# median filter to smooth the signal 
median_filter_width = 5  # should be lower for very peaky breakthrough curves 
mmed <- function(x,n=5){runmed(x,n)} #Median

# NA remove -----
rm_nas = function(x){
  tmp_rm = which(is.na(x$X1))
  if(length(tmp_rm)>0){
    x <- x[-tmp_rm,]  
  }
  return(x)
}

# JAVA heap space cleaner
jgc <- function()
{
  .jcall("java/lang/System", method = "gc")
}  



# Read XLSX
read_xlsx_ = function(filename,sheet){
  # filename = filename_excel
  # sheet = sheet_name
  d_ = readxl::read_xlsx(path = filename,sheet = sheet,col_names = F)
  names(d_) = paste0("X", seq(NCOL(d_)))
  
  # Extract date and time parts
  correct_date <- as.Date(d_$X1)           # Extract the date part from date column
  correct_time <- format(d_$X3, "%H:%M:%S") # Extract the time part as a string
  
  # Combine date and time into a single datetime
  d_$X1 <- as.POSIXct(paste(correct_date, correct_time), format = "%Y-%m-%d %H:%M:%S")
  
  return(d_)
}



# make date-time -----
make_date_time = function(x){
  # 1 - make sure the manual time shift is not creating a 1 day gap:
  tmp_tm_orig = format(x$X2, '%H:%M:%S')  # 3rd column local time
  tmp_tm_orig_h = as.numeric(format(x$X2, '%H'))  # 3rd column local time
  if(any(diff(tmp_tm_orig_h) < -10)){
    print('Time-series is extending over midnight; checking if simple date (day) adjustment is sufficient to fix it.')
    tmp_yr = format(x$X1, '%Y-%m-%d')  # 1st column year
    tmp_day= as.numeric(format(x$X1, '%d'))
    max_day = range(tmp_day)[2]
    min_day = range(tmp_day)[1]
    tmp_days_to_change = which(tmp_day == max_day)
    tmp_tm = format(x$X3, '%H:%M:%S')  # 3rd column local time
    to_rep = sprintf("-%02d", max_day)
    to_repwith = sprintf("-%02d", min_day)
    
    tmp_yr[tmp_days_to_change] = gsub (pattern =to_rep, replacement = to_repwith,x =  tmp_yr[tmp_days_to_change])
    tmp_date = as.POSIXct(paste(tmp_yr, tmp_tm))
    if(max(diff(tmp_date)) > (10*60)){
      print('Something went wrong. Check please. Time difference still more than 10 hours.')
    }else{
      print('Adjusting the day of the following day (after midnight) seems to have worked.')
    }
  }else{
    # 2
    print('Time series looks fine.')
    tmp_yr = format(x$X1, '%Y-%m-%d')  # 1st column year
    tmp_tm = format(x$X3, '%H:%M:%S')  # 3rd column local time
    tmp_date = as.POSIXct(paste(tmp_yr, tmp_tm))
  }
  return(tmp_date)
}

# -------- Noise removal ------- 
remove_noise = function(x){
  # x = d
  plot(x$date, x$Tracer1, type='l', ylab='Signal [mV]')
  abline(h=0, lty=2, col='red')
  print("Select first point (1/2) in the plot (click) to identify background noise ... Press [Esc] to interrupt!")
  pn_12_1 = identify(x$date, x$Tracer1, n = 1, col='darkgreen', pch=16)
  points(x$date[pn_12_1],x$Tracer1[pn_12_1], col='darkgreen', pch=1, cex=2, lwd=3)
  print("Select second point (2/2) in the plot (click) to identify background noise ... Press [Esc] to interrupt!")
  pn_12_2 = identify(x$date, x$Tracer1, n = 1, col='darkgreen', pch=16)
  points(x$date[pn_12_2],x$Tracer1[pn_12_2], col='darkgreen', pch=1, cex=2, lwd=3)
  
  
  mn = mean(x$Tracer1[pn_12_1:pn_12_2])
  print(paste0("Background noise signal is: ", round(mn,2), ' mV'))
  x$Tracer1 <- x$Tracer1 - mn
  print(paste0("Background noise removed."))
  
  # median filter 1
  x$Tracer1_mmed9 = mmed(x$Tracer1, median_filter_width)
  # plot(d$date[1:200], d$Tracer1_mmed9[1:200] , col=1, type='l')
  plot(x$date, x$Tracer1_mmed9 , col=1, type='l')
  abline(h=0,col='red')
  return(x)
}

# -------- Calibration -------
calibrate = function(x,start,stop,cs1,cs2,cs3, bucketvolume=5){
  # x = d
  print(paste0("---- Identifying calibration levels ----"))
  # seq_cal = seq(1,length(x$Tracer1_mmed9))  # all points
  seq_cal = seq(start,stop)                      # some points at beginning - if difficult to see the calib
  plot(x$date[seq_cal], x$Tracer1_mmed9[seq_cal] , col=1, type='l', xlab='Time', ylab='Signal [mV]')
  abline(h=0,col='red')
  
  c1 = c2 = c3 = FALSE
  p1_12_2 = p2_12_2 = p3_12_2 = NULL
  m1 = m2 = m3 = NULL
  
  # ------ Calibration solution 1 -------
  print("Select first point (1/2) in the plot (click) to identify first calibration level (same plateau)... Press [Esc] to interrupt!")
  p1_12_1 = identify(x$date[seq_cal], x$Tracer1_mmed9[seq_cal], n = 1, col='red')
  points(x$date[p1_12_1],x$Tracer1_mmed9[p1_12_1], col='red', pch=1, cex=2, lwd=3)
  print("Select second point (2/2) in the plot (click) to identify first calibration level (same plateau)... Press [Esc] to interrupt!")
  p1_12_2 = identify(x$date[seq_cal], x$Tracer1_mmed9[seq_cal], n = 1, col='red')
  points(x$date[p1_12_2],x$Tracer1_mmed9[p1_12_2], col='red', pch=1, cex=2, lwd=3)
  
  if(length(p1_12_2) > 0){
    m1 = mean(x$Tracer1_mmed9[seq_cal][p1_12_1:p1_12_2])
    print(paste0("Calibration signal (first level) is: ", round(m1,2), ' mV'))
    c1 = T
  }else{
    print('No first calibration level determined.')
  }
  
  
  # ------ Calibration solution 2 -------
  print("Select first point (1/2) in the plot (click) to identify second calibration level (same plateau)... Press [Esc] to interrupt!")
  p2_12_1 = identify(x$date[seq_cal], x$Tracer1_mmed9[seq_cal], n = 1, col='orange')
  points(x$date[p2_12_1],x$Tracer1_mmed9[p2_12_1], col='orange', pch=1, cex=2, lwd=3)
  
  print("Select second point (2/2) in the plot (click) to identify second calibration level (same plateau)... Press [Esc] to interrupt!")
  p2_12_2 = identify(x$date[seq_cal], x$Tracer1_mmed9[seq_cal], n = 1, col='orange')
  points(x$date[p2_12_2],x$Tracer1_mmed9[p2_12_2], col='orange', pch=1, cex=2, lwd=3)
  
  if(length(p2_12_2) > 0){
    m2 = mean(x$Tracer1_mmed9[seq_cal][p2_12_1:p2_12_2])
    print(paste0("Calibration signal (second level) is: ", round(m2,2), ' mV'))
    c2 = T
  }else{
    print('No second calibration level determined.')
  }
  
  # ------ Calibration solution 3 -------
  print("Select first point (1/2) in the plot (click) to identify third calibration level (same plateau)... Press [Esc] to interrupt!")
  p3_12_1 = identify(x$date[seq_cal], x$Tracer1_mmed9[seq_cal], n = 1, col='purple')
  points(x$date[p3_12_1],x$Tracer1_mmed9[p3_12_1], col='purple', pch=1, cex=2, lwd=3)
  
  print("Select second point (2/2) in the plot (click) to identify third calibration level (same plateau)... Press [Esc] to interrupt!")
  p3_12_2 = identify(x$date[seq_cal], x$Tracer1_mmed9[seq_cal], n = 1, col='purple')
  points(x$date[p3_12_2],x$Tracer1_mmed9[p3_12_2], col='purple', pch=1, cex=2, lwd=3)
  
  
  if(length(p3_12_2) > 0){
    m3 = mean(x$Tracer1_mmed9[seq_cal][p3_12_1:p3_12_2])
    print(paste0("Calibration signal (third level) is: ", round(m3,2), ' mV'))
    c3 = T
  }else{
    print('No third calibration level determined.')
  }
  
  
  # --------------------- Checking calibration -------------
  n_c = sum(c(c1,c2,c3))
  print(paste0('Number of calibration levels: ',n_c))
  
  # real concentration cs = c(40, 79, 119, 158)
  # 2025-07-06:
  # These are now input parameters!
  if(length(cs1)<1){
    print('Missing mass of calibration solution 1. Please do so in the setup! Stopping.')
    break
  }
  if(length(cs2)<1){
    print('Missing mass of calibration solution 2. Using only one point for calibration.')
  }
  if(length(cs3)<1){
    print('Missing mass of calibration solution 3. Using only one point for calibration.')
  }
  if(length(cs4)<1){
    print('Missing mass of calibration solution 4. Using only one point for calibration.')
  }
  # For 1 calibration point
  if(n_c == 1 & c1){
    bucketvolume_ = bucketvolume + m_calib_solution1/1000
    cs1_ = cs1 *10^6 / bucketvolume_   #c(49.5)
    cs = c(cs1_)
    mVs = c(m1)
    print(paste0('Using one calibration point for concetration of ', cs, ' ppb and a signal of ',round(mVs,2), ' mV'))
  }else if(n_c == 2 & c1 & c2){
    # For 2 calibration point
    bucketvolume_ = bucketvolume + m_calib_solution1/1000
    cs1_ = cs1 *10^6 / bucketvolume_   
    bucketvolume_ = bucketvolume_ + m_calib_solution2/1000
    cs2_ = (cs1 *10^6 + cs2 *10^6) / (bucketvolume_)
    cs = c(cs1_, cs2_)   
    mVs = c(m1,m2)
    print(paste0('Using two calibration points for concetrations of ', paste(cs, collapse = ' & '),
                 ' ppb and signals of ',paste(round(mVs,2), collapse = ' & '), ' mV'))
  }else if(n_c == 3 & c1 & c2 & c3){
    # For 3 calibration point
    # cs = c(49.5, 98.04, 145.63)
    bucketvolume_ = bucketvolume + m_calib_solution1/1000
    cs1_ = cs1 *10^6 / bucketvolume_   
    bucketvolume_ = bucketvolume_ + m_calib_solution2/1000
    cs2_ = (cs1 *10^6 + cs2 *10^6) / (bucketvolume_)
    bucketvolume_ = bucketvolume_ + m_calib_solution3/1000
    cs3_ = (cs1 *10^6 + cs2 *10^6 + cs3 *10^6) / (bucketvolume_)
    
    cs = c(cs1_, cs2_, cs3_)  
    
    mVs = c(m1,m2,m3) 
    print(paste0('Using three calibration points for concetrations of ', paste(cs, collapse = ' & '),
                 ' ppb and signals of ',paste(round(mVs,2), collapse = ' & '), ' mV'))
  }else if(n_c == 2 & c1 & c3){
    # For 2 calibration points (rejecting the point at cs2 for regression)
    # cs = c(49.5,  145.63)
    bucketvolume_ = bucketvolume + m_calib_solution1/1000
    cs1_ = cs1 *10^6 / bucketvolume_   
    bucketvolume_ = bucketvolume_ + m_calib_solution2/1000
    cs2_ = (cs1 *10^6 + cs2 *10^6) / (bucketvolume_)
    bucketvolume_ = bucketvolume_ + m_calib_solution3/1000
    cs3_ = (cs1 *10^6 + cs2 *10^6 + cs3 *10^6) / (bucketvolume_)
    
    cs = c(cs1_, cs3_)  
    
    mVs = c(m1,m3) 
    print(paste0('Using two calibration points for concetrations of ', paste(cs, collapse = ' & '),
                 ' ppb and signals of ',paste(round(mVs,2), collapse = ' & '), ' mV'))
  }else if(n_c == 2 & c2 & c3){
    # For 2 calibration points (rejecting the point at cs2 for regression)
    # cs = c(49.5,  145.63)
    bucketvolume_ = bucketvolume + m_calib_solution1/1000
    cs1_ = cs1 *10^6 / bucketvolume_   
    bucketvolume_ = bucketvolume_ + m_calib_solution2/1000
    cs2_ = (cs1 *10^6 + cs2 *10^6) / (bucketvolume_)
    bucketvolume_ = bucketvolume_ + m_calib_solution3/1000
    cs3_ = (cs1 *10^6 + cs2 *10^6 + cs3 *10^6) / (bucketvolume_)
    
    cs = c(cs2_, cs3_)  
    
    mVs = c(m2,m3) 
    print(paste0('Using two calibration points for concetrations of ', paste(cs, collapse = ' & '),
                 ' ppb and signals of ',paste(round(mVs,2), collapse = ' & '), ' mV'))
  }else{
    if(exists("LM1")){
      print(paste0('Could not find calibration points. Forgot? Will use existing relationship from before...'))
    }else{
      print(paste0('Could not find any calibration. Forgot? This is neccessary! Using the latest calibration I can find ...'))  
      load("_LM1.RData")
      load("_cs.RData")
      load("_mVs.RData")
    }
  }
  if(length(mVs) > 0){
    assign("mVs",mVs,envir = .GlobalEnv)
    assign("cs",cs,envir = .GlobalEnv)
  }
}


# Create linear model to relate mV to ppb --------
relate_mV_ppb = function(x){
  opar <- par()
  par(mfrow=c(2,1), mar=c(4,4,.1,.1))
  plot(mVs, cs, type='p', xlab = 'Signal [mV]' ,  #, xlim = c(0,2500), ylim=c(0,250)
       ylab='Concentration [µg/l = ppb]')
  lines(mVs, cs)
  LM1 = lm(cs~mVs-1)
  summary(LM1)
  predLM = predict(LM1, newdata = data.frame(mVs = mVs))
  lines(mVs, predLM, col='blue')
  abline(LM1, col='red', lty=2)
  txt = paste0("R2: ", round(summary(LM1)$r.squared,4))
  mtext(txt,side = 3,line = -2,adj = 0.05)
  
  # ----- conversion  ----- #
  x$adj = predict(LM1, newdata = data.frame(mVs = x$Tracer1))
  plot(x$date, x$adj, type='l', ylab='Concentration [µg/l = ppb]', xlab='Time')
  abline(h=cs[1], col='blue', lty=3, lwd=.5)
  abline(h=cs[2], col='blue', lty=3, lwd=.5)
  abline(h=cs[3], col='blue', lty=3, lwd=.5)
  
  # overwrite "adj" with median filtered value
  x$adj = mmed( x$adj, median_filter_width)
  lines(x$date, x$adj, col='red')
  par(opar)
  assign("LM1",LM1,envir = .GlobalEnv)
  # save model in case heap-space error
  save(LM1,file = '_LM1.RData')
  save(mVs,file = '_mVs.RData')
  save(cs,file = '_cs.RData')
  return(x)
}


# Calculate discharge ------
calc_discharge = function(x, start=NULL, stop=NULL){
  # ----- breakthrough  ----- Run #1
  # start = 700
  # stop = NULL
  # start = 1500
  # stop=2500
  # x=d
  if(length(start)>0){
    if(length(stop)==0){
      stop=length(x$date)
      if(start>stop){
        print("Start point too big (higher than entire series!). Choose earlier point!")
        break
      }
    }
    plot(x$date[start:stop], x$adj[start:stop], type = 'l', ylab='Concentration [ppb]', xlab='Time')
    
    print("Select start point (1/2) in the plot (click) of breakthrough curve ... Press [Esc] to interrupt!")
    pb_12_1 = identify(x$date[start:stop], x$adj[start:stop], n = 1, col='brown')
    points(x$date[start:stop][pb_12_1],x$adj[start:stop][pb_12_1], col='brown', pch=1, cex=2, lwd=3)
    print("Select end point (2/2) in the plot (click) of breakthrough curve ... Press [Esc] to interrupt!")
    pb_12_2 = identify(x$date[start:stop], x$adj[start:stop], n = 1, col='brown')
    points(x$date[start:stop][pb_12_2],x$adj[start:stop][pb_12_2], col='brown', pch=1, cex=2, lwd=3)
    # re-plot
    dss = x[start:stop,][pb_12_1:pb_12_2,]
  }else{
    plot(x$date, x$adj, type = 'l', ylab='Concentration [ppb]', xlab='Time')
    
    print("Select start point (1/2) in the plot (click) of breakthrough curve ... Press [Esc] to interrupt!")
    pb_12_1 = identify(x$date, x$adj, n = 1, col='brown')
    points(x$date[pb_12_1],x$adj[pb_12_1], col='brown', pch=1, cex=2, lwd=3)
    print("Select end point (2/2) in the plot (click) of breakthrough curve ... Press [Esc] to interrupt!")
    pb_12_2 = identify(x$date, x$adj, n = 1, col='brown')
    points(x$date[pb_12_2],x$adj[pb_12_2], col='brown', pch=1, cex=2, lwd=3)
    
    # re-plot
    dss = x[pb_12_1:pb_12_2,]
  }
  
  
  plot(dss$date, dss$adj, type = 'l', ylab='Concentration [ppb]', xlab='Time')
  abline(h=0, col='red', lty=2)
  
  # ---- log fit -----
  # choose last 10 or 20% of values to fit asymptotic log curve.
  
  n_dss = length(dss$date)
  n_fit = round(n_dss*0.2,0)
  v_fit = dss$adj[(n_dss-n_fit):n_dss]
  fit_df = data.frame(v = log(v_fit),x =seq(length(v_fit)))
  
  plot(dss$adj, type='l', xlim = c(0,20*n_fit+n_dss), ylim=c(0,10), ylab='Concentration [ppb]', xlab='Time')
  lm_fit = lm(v~x, data = fit_df)
  nwdat = data.frame(x=seq(20*n_fit))
  lm_pred_log = predict(lm_fit,newdata = nwdat)
  lm_pred = exp(lm_pred_log)
  lines(seq(length(lm_pred))+n_dss-n_fit,lm_pred, col='red')
  
  
  # without extrapolation
  dt_sec = cumsum(c(0,diff(dss$date)))
  A = auc(dt_sec, dss$adj)
  Q1 = M/A
  
  
  # [1] 0.01978055
  # with extrapolation
  adj_extr = c(dss$adj, lm_pred[n_fit:length(lm_pred)])
  date_adj_extr = rep(diff(dss$date)[1], length(adj_extr))
  # cumsum only with numeric / not datetime
  date_adj_extr = as.numeric(date_adj_extr)
  dt_sec_extr =  cumsum(c(date_adj_extr))
  A = auc(dt_sec_extr, adj_extr)
  
  Q2 = M/A
  Q2
  
  plot(dss$date, dss$adj, type = 'l', ylab='Concentration [ppb]', xlab='Time')
  pp_x = c(dss$date, rev(dss$date))
  pp_y = c(dss$adj, rep(0,length(dss$adj)))
  print(paste0("Parameters: Tracer mass=",M,' mg  ; Standard=',Standard,' g/l  ;  Tracer volume=',TracerMass,' g'))
  print(paste0("Calculated discharge is: ",round(Q1,5)," m^3/s"))
  print(paste0("Calculated discharge (with extrapolation) is: ",round(Q2,5)," m^3/s"))
  polygon(pp_x,pp_y, col='#3388AA66')
  assign("dss",dss,envir = .GlobalEnv)
  assign("n_fit",n_fit,envir = .GlobalEnv)
  assign("lm_fit",lm_fit,envir = .GlobalEnv)
  assign("Q_no_extrapol",Q1,envir = .GlobalEnv)
  assign("Q_with_extrapol",Q2,envir = .GlobalEnv)
  # 0.01991656 m3/s
}


auc2 = function (x, y, from = min(x, na.rm = TRUE), to = max(x, na.rm = TRUE), 
                type = c("linear", "spline"), absolutearea = FALSE, subdivisions = 100, 
                ...) 
{
  type <- match.arg(type)
  stopifnot(length(x) == length(y))
  stopifnot(!is.na(from))
  if (length(unique(x)) < 2) 
    return(NA)
  if (type == "linear") {
    if (absolutearea == FALSE) {
      values <- approx(x, y, xout = sort(unique(c(from, 
                                                  to, x[x > from & x < to]))), ...)
      res <- 0.5 * sum(diff(values$x) * (values$y[-1] + 
                                           values$y[-length(values$y)]))
    }
    else {
      o <- order(x)
      ox <- x[o]
      oy <- y[o]
      idx <- which(diff(oy >= 0) != 0)
      newx <- c(x, x[idx] - oy[idx] * (x[idx + 1] - x[idx])/(y[idx + 
                                                                 1] - y[idx]))
      newy <- c(y, rep(0, length(idx)))
      values <- approx(newx, newy, xout = sort(unique(c(from, 
                                                        to, newx[newx > from & newx < to]))), ...)
      res <- 0.5 * sum(diff(values$x) * (abs(values$y[-1]) + 
                                           abs(values$y[-length(values$y)])))
    }
  }
  else {
    if (absolutearea) 
      myfunction <- function(z) {
        abs(splinefun(x, y, method = "natural")(z))
      }
    else myfunction <- splinefun(x, y, method = "natural")
    res <- integrate(myfunction, lower = from, upper = to, 
                     subdivisions = subdivisions)$value
  }
  res
}

# Calculate discharge with overlapping curves -----
calc_discharge_overlap = function(x, xss){
  # - Preparation for curve 2 - #
  # do the fit also for the rest of the time series to adjust potential second
  # injection for remaining tracer from first injection in the stream !
  tnd = length(xss$date)  # get last point of Inj1
  tndlw = length(xss$date) - n_fit  # starting point of extrapolation
  tind = which(x$date == xss$date[tnd])
  tindlw = which(x$date == xss$date[tndlw])
  tnright = length(x$date) - tind
  tnrightlw = length(x$date) - tindlw
  nwdat = data.frame(x=seq(tnrightlw))
  lm_pred_log = predict(lm_fit,newdata = nwdat)
  lm_pred = exp(lm_pred_log)
  # xx = d$date[(tindlw+1):length(d$date)]
  # yy = lm_pred
  # lines(xx,yy, col='red' )
  # - fill original data until prediction starts - and then fill with predictions
  x$Tracer1bgc1 = x$adj
  x$Tracer1bgc1[(tindlw+1):length(x$Tracer1bgc1)] = lm_pred
  x$Tracer1cv2_adj = x$adj - x$Tracer1bgc1
  plot(x$Tracer1cv2_adj)
  # !Preparation for curve 2 - #
  
  # !------ Curve 2 ------- !
  
  # ----- breakthrough  ----- Run #1
  plot(x$date, x$Tracer1cv2_adj, type = 'l', ylab='Concentration [ppb]', xlab='Time')
  
  print("Select start point (1/2) in the plot (click) of breakthrough curve ... Press [Esc] to interrupt!")
  pb_12_1 = identify(x$date, x$Tracer1cv2_adj, n = 1, col='brown')
  points(x$date[pb_12_1],x$Tracer1cv2_adj[pb_12_1], col='brown', pch=1, cex=2, lwd=3)
  print("Select end point (2/2) in the plot (click) of breakthrough curve ... Press [Esc] to interrupt!")
  pb_12_2 = identify(x$date, x$Tracer1cv2_adj, n = 1, col='brown')
  points(x$date[pb_12_2],x$Tracer1cv2_adj[pb_12_2], col='brown', pch=1, cex=2, lwd=3)
  
  
  # plot(x$date, x$Tracer1cv2_adj, type = 'l')
  # pbt_12 = identify(x$date, x$Tracer1cv2_adj, n = 2)
  
  # re-plot
  dss = x[pb_12_1:pb_12_2,]
  
  plot(dss$date, dss$Tracer1cv2_adj, type='l', ylab='Concentration [ppb]', xlab='Time')
  abline(h=0, col='red', lty=2)
  
  # ---- log fit -----
  # choose last 10 or 20% of values to fit asymptotic log curve.
  
  n_dss = length(dss$date)
  n_fit = round(n_dss*0.2,0)
  v_fit = dss$Tracer1cv2_adj[(n_dss-n_fit):n_dss]
  fit_df = data.frame(v = log(v_fit),x =seq(length(v_fit)))
  
  plot(dss$Tracer1cv2_adj, type='l', xlim = c(0,20*n_fit+n_dss), ylim=c(0,10))
  lm_fit = lm(v~x, data = fit_df)
  nwdat = data.frame(x=seq(20*n_fit))
  lm_pred_log = predict(lm_fit,newdata = nwdat)
  lm_pred = exp(lm_pred_log)
  lines(seq(length(lm_pred))+n_dss-n_fit,lm_pred, col='red')
  
  # CALC Q
  # without extrapolation
  dt_sec = cumsum(c(0,diff(dss$date)))
  A = auc2(dt_sec, dss$Tracer1cv2_adj)
  Q1 = M/A
  
  # with extrapolation
  adj_extr = c(dss$Tracer1cv2_adj, lm_pred[n_fit:length(lm_pred)])
  date_adj_extr = rep(diff(dss$date)[1], length(adj_extr))
  # cumsum only with numeric / not datetime
  date_adj_extr = as.numeric(date_adj_extr)
  dt_sec_extr =  cumsum(c(date_adj_extr))
  A = auc(dt_sec_extr, adj_extr)
  
  Q2 = M/A

  plot(dss$date, dss$Tracer1cv2_adj, type = 'l', ylab='Concentration [ppb]', xlab='Time')
  pp_x = c(dss$date, rev(dss$date))
  pp_y = c(dss$Tracer1cv2_adj, rep(0,length(dss$Tracer1cv2_adj)))
  print(paste0("Parameters: Tracer mass=",M,' mg  ; Standard=',Standard,' g/l  ;  Tracer volume=',TracerVolume,' l'))
  print(paste0("Calculated discharge is: ",round(Q1,5)," m^3/s"))
  print(paste0("Calculated discharge (with extrapolation) is: ",round(Q2,5)," m^3/s"))
  polygon(pp_x,pp_y, col='#AA883366')
  assign("dss",dss,envir = .GlobalEnv)
  assign("Q2_no_extrapol",Q1,envir = .GlobalEnv)
  assign("Q2_with_extrapol",Q2,envir = .GlobalEnv)
  
}

# function to read the analogue digital converted data from the FL30

bit_to_mv = function(data_file, sep=NULL){
  # data_file = 'tests/_20210701-03-f1664.txt'
  # data_file = '/Users/pohle/Dropbox/Central_Asia/Fluorometer/2024/Abramov/data/f1664_2024_Abr.txt'
  # data_file = paste0(path, filename_raw)
  # data_file = '/Users/pohle/Dropbox/FL30/CALIBRAT_1993/backup_2025-07-06/test/part1.txt'
  # sep=','
  # data_file = paste0(path, filename_raw)
  # # ---------------- RAW ------------------------- #
  # bd = read.table('tests/_20210701-03-f1664.txt')
  # data_file = "/Users/pohle/Dropbox/Central_Asia/Fluorometer/2023/Zulmart/data/f900 2_2.csv"
  # data_file = "/Users/pohle/Dropbox/Central_Asia/HOBO/summer_school_Suyek/data/FL30-binary_f1664.txt"
  if(length(sep) == 0){
    bd = read.table(data_file)  
  }else if(sep==','){
    bd = read.csv(data_file, header = F)
  }
  # }else{
  #   bd = read.table(data_file)  
  # }
  

  # str(bd)
  # head(bd)
  # nrow(bd)
  # which(is.na(as.numeric(bd$V9)))
  # bd[8831:8834,]
  # bd <- bd[-1,]  # remove first line (blank)
 
  # date and time
  dt_string = paste(bd$V2, bd$V3)
  datetime = as.POSIXct(dt_string, format = '%d:%m:%y %H:%M:%S', tz='UTC')
  
  
  # signal:
  # s = 2^16 * v1 + v2
  # with v1 (column7), and v2 (column8)
  
  
  # Baseline
  vt21 = bd$V5
  vt22 = bd$V6
  st2 = 2^16*vt21 + vt22
  smvt2 = round(st2 * 2500/2^24, 2)
  
  # Tracer 1
  v1 = bd$V7
  v2 = bd$V8
  if(!is.numeric(v2)){
    print("Some problem with the binary file - non numeric data in column 8 - please check \n Transforming data into numeric for now - Be Aware!!!")
    v2 = as.numeric(bd$V8)
  }
  s = 2^16*v1 + v2
  smv = round(s * 2500/2^24, 2)
  
  
  # Voltage
  # battery = bit*(float)((v0<<16)+v1)*6.53/1076.8
  sv = (2^16*bd$V9 + bd$V10 ) *6.53/1076.8
  # plot(sv)
  svmv = round(sv *2500/2^24, 2)
  # plot(svmv)

  # Temperature
  st = 2^16*bd$V11 + bd$V12
  stmv = round((st * 2500/2^24)/10, 2)
  # plot(stmv/10)
  
  # output df
  out_df = data.frame(date=datetime, Tracer1 = smv, Battery=svmv, Baseline=smvt2, Temperature=stmv)
  return(out_df)  
}

bit_to_mv_30sec = function(data_file){
  # data_file = '/Users/pohle/Dropbox/Central_Asia/Fluorometer/2022/raw/f1664_30sec-Abramov.txt'
  # data_file = '/Users/pohle/Dropbox/Central_Asia/Fluorometer/2022/raw/f1664_Abramov.txt'
  # ---------------- RAW ------------------------- #
  # bd = read.table('tests/_20210701-03-f1664.txt')
  bd = read.table(data_file)
  bd <- bd[-1,]  # remove first line (blank)
  # head(bd)
  # date and time
  dt_string = paste(bd$V2, bd$V3)
  datetime = as.POSIXct(dt_string, format = '%d:%m:%y %H:%M:%S', tz='UTC')
  
  
  # signal:
  # s = 2^16 * v1 + v2
  # with v1 (column7), and v2 (column8)
  
  
  # Baseline
  vt21 = bd$V5
  vt22 = bd$V6
  st2 = 2^16*vt21 + vt22
  smvt2 = round(st2 * 2500/2^24, 2)
  
  # Tracer 1
  v1 = bd$V7
  v2 = bd$V8
  s = 2^16*v1 + v2
  smv = round(s * 2500/2^24, 2)
  # plot(smv)
  
  
  # Voltage
  # battery = bit*(float)((v0<<16)+v1)*6.53/1076.8
  sv = (2^16*bd$V9 + bd$V10 ) *6.53/1076.8
  # plot(sv)
  svmv = round(sv *2500/2^24, 2)
  # plot(svmv)
  
  # Temperature
  st = 2^16*bd$V11 + bd$V12
  stmv = round((st * 2500/2^24)/10, 2)
  # plot(stmv/10)
  
  # --- test ----- #
  # There are multiple columns - check what they are in the 30 sec file
  head(bd)
  v1 = bd$V9
  v2 = bd$V10
  s = 2^16*v1 + v2
  xxx = round(s * 2500/2^24, 2)
  plot(xxx)
  
  v1 = bd$V11
  v2 = bd$V12
  s = 2^16*v1 + v2
  xxx = round(s * 2500/2^24, 2)
  plot(xxx)
  
  v1 = bd$V13
  v2 = bd$V14
  s = 2^16*v1 + v2
  xxx = round(s * 2500/2^24, 2)
  plot(xxx)
  
  v1 = bd$V15
  v2 = bd$V16
  s = 2^16*v1 + v2
  xxx = round(s * 2500/2^24, 2)
  plot(xxx)
  
  # same as V9 + V10
  v1 = bd$V17
  v2 = bd$V18
  s = 2^16*v1 + v2
  xxx = round(s * 2500/2^24, 2)
  plot(xxx)
  
  v1 = bd$V19
  v2 = bd$V20
  s = 2^16*v1 + v2
  xxx = round(s * 2500/2^24, 2)
  plot(xxx)
  
  v1 = bd$V21
  v2 = bd$V22
  s = 2^16*v1 + v2
  xxx = round(s * 2500/2^24, 2)
  plot(xxx)
  
  v1 = bd$V23
  v2 = bd$V24
  s = 2^16*v1 + v2
  xxx = round(s * 2500/2^24, 2)
  plot(xxx)
  
  v1 = bd$V25
  v2 = bd$V26
  s = 2^16*v1 + v2
  xxx = round(s * 2500/2^24, 2)
  plot(xxx)
  
  v1 = bd$V27
  v2 = bd$V28
  s = 2^16*v1 + v2
  xxx = round(s * 2500/2^24, 2)
  plot(xxx)
  # output df
  out_df = data.frame(date=datetime, Tracer1 = smv, Battery=svmv, Baseline=smvt2, Temperature=stmv)
  return(out_df)  
}


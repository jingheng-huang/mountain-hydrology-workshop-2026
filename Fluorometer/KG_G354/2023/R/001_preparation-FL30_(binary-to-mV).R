# - reading FL30 raw data
# - transforming data into readable format with signal in mV

# Arguments
path = "/Users/pohle/Dropbox/Central_Asia/Fluorometer/2023/G354/data/"
# For Windows it will be something like this:
# path = "C:/summer_school_Suyek/data/"

# Input files
filename_raw = 'f1664_g354.txt'

# Output files
filename_mV = 'f1664.mv'



# ------------ NO NEED TO CHANGE ANYTHING HEREAFTER --------------- #
setwd(path)
source('../R/FUN_install_packages.R')
options(java.parameters = "-Xmx2000m")  # allow for 2GB RAM - for JAVA
source("../R/FUN_binary_to_mv.R")

aaa = bit_to_mv(data_file = paste0(path, filename_raw),sep=' ')
tail(aaa)
head(aaa)
opar <- par()
par(mfrow=c(2,1), mar=c(4,4,.2,.2))
plot(aaa$Tracer1, type='l', xlab='Time [s]', ylab='Signal [mV]')

plot(aaa$date,aaa$Tracer1, type='l', xlab='Date', ylab='Signal [mV]')
par(opar)

write.table(aaa,file = paste0(path,filename_mV),
            quote = F, row.names = F)

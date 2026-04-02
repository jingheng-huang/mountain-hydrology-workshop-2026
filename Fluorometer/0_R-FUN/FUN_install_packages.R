# This script installs all required packages
list.of.packages <- c("MESS", "pracma", "xlsx", "xts", "lubridate")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){install.packages(new.packages)}

# load packages:
require('MESS') # for AUC
require("xts")
require("xlsx")
require("pracma")
require("lubridate")


library(Seurat)
library(tidyverse)
library(future)
plan()

options(future.globals.maxSize =1024*1024^2)
sham1 <- 
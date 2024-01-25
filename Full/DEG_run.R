library(DEGman)
library("Seurat")
library("dplyr")
library("hdf5r")
library(philentropy)
library(foreach)
library(doParallel)
source('distribution.R')

hc1= Read10X_h5('/home/cuiyaxuan/spatialLIBD/151673/151673_filtered_feature_bc_matrix.h5') #### to your path and project name
label=read.csv("/home/cuiyaxuan/Zero/DenoiseST-master2/Full/label.csv",header = T,row.names = 1) # cluster label
k=1 ##### Default parameters
dis<-distri(hc1,label,k)


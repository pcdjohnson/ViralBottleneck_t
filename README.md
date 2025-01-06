# ViralBottleneck
## Description
This package is used for estimating viral transmission bottleneck size using different methods.
The main steps to use the program are:
1. Create Transmission Object
2. Calculate Transmission bottleneck size 
   
There are the he [manual](manual_and_tutorial/ViralBottleneck_manual_0.1.0.pdf) and [turtorial](manual_and_tutorial/Tutorial.pdf).

The test data for "ViralBottleneck" is in [test_dataset folder](test_dataset).

## Download
To install the package in R use the following commands: 
```
library(devtools)
install_github("BowenArchaman/ViralBottleneck_t")
library(ViralBottleneck)
```

## Requirements
- R 4.2.2
- methods
- utils
- ggplot2
- pbapply
- rmutil

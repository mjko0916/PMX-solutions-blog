# A step-by-step guide to goodness of fit figures of NONMEM models in R using ggplot2

## **GOF(Goodness Of Fit)**  
- Population model predictions vs observations  
- Individual model predictions vs observations  
- Conditional weighted residuals vs population predictions    
- Conditional weighted residuals vs time  
- Individual + population model prediction and observations over time per individual

## libraries
```{r}
library(xpose4)
library(tidyverse)
library(ggforce)
```

## data
- ID
- TIME
- DV : observations
- IPRED or IPRE : individual prediction
- CPREDI : FO conditional prediction with full interaction
- CRESI : FO conditional residuals with full interaction
- CWRESI : FO conditional weighted residuals with full interaction
- PRED : population prediction  

## 1. Population/individual model predictions versus observations
```{r}
df <- read.table(./SDTABModel001.txt, skip=1, sep="")
head(df)
```




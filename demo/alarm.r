library(catnet)

## set the stored random generator kind and its seed
data(alarmseed)
saved.rng <- RNGkind()
RNGkind(kind="Mersenne-Twister")

## load the alarm data sample and ground-thruth network
data(alarm)
data(alarmnet)
alarmnet

##alarmseed <- .Random.seed
##save(alarmseed, file="alarmseed.rda")

eval <- cnEvaluate(object=alarmnet, data=alarm, perturbations=NULL, maxParentSet=3)
eval
##pdf("alarmeval.pdf")
cnPlot(eval)
##dev.off()

bstnet<- cnFind(eval, cnComplexity(alarmnet))
bstnet
cnCompare(alarmnet, bstnet)

set.seed(alarmseed)
sanets1 <- cnSearchSA(data=alarm, perturbations=NULL, 
                      maxParentSet=2, parentSizes=NULL, maxComplexity=600,
                      parentsPool=NULL, fixedParentsPool=NULL, edgeProb=NULL,
                      selectMode="BIC", 
                      tempStart=1, tempCoolFact=1, tempCheckOrders=100, maxIter=100,
                      orderShuffles=0, stopDiff=0,
                      numThreads = 2, 
                      priorSearch=NULL, echo = FALSE)
sanets1
sanet1 <- cnFindBIC(sanets1)
cnCompare(alarmnet, sanet1)
sanet1@meta <- "sanet1, SA learning, Stage I"
sanet1
cnPlot(sanet1, "sanet1")

set.seed(alarmseed)
sanets <- cnSearchSA(data=alarm, perturbations=NULL, 
                     maxParentSet=2, parentSizes=NULL, maxComplexity=600,
                     parentsPool=NULL, fixedParentsPool=NULL, edgeProb=NULL,
                     selectMode="BIC", 
                     tempStart=1e-2, tempCoolFact=0.9, tempCheckOrders=20, maxIter=200,
                     orderShuffles=4, stopDiff=1e-3,
                     numThreads = 4, 
                     priorSearch=sanets1, echo=TRUE)
sanets
sanet2 <- cnFindBIC(sanets)
cnCompare(alarmnet, sanet2)
sanet2@meta <- "sanet2, SA learning, Stage II"
sanet2
cnPlot(sanet2, "sanet2")

## restore the random generator
RNGkind(kind=saved.rng[1])

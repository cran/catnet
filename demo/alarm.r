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
cnPlot(eval)

bstnet<- cnFind(eval, cnComplexity(alarmnet))
bstnet
cnCompare(alarmnet, bstnet)

set.seed(alarmseed)
sanets <- cnSearchSA(data=alarm, perturbations=NULL, 
			maxParentSet=2, maxComplexity=600,
			parentsPool=NULL, fixedParentsPool=NULL,
			selectMode="BIC", 
			tempStart=10, tempCoolFact=1, tempCheckOrders=100, maxIter=100,
			orderShuffles=0, stopDiff=0.0000000001,
			priorSearch=NULL)
sanets
sanet1 <- cnFindBIC(sanets, nrow(alarm))
cnCompare(alarmnet, sanet1)
sanet1@meta <- "sanet1, SA learning, Stage I"
sanet1
cnPlot(sanet1, "sanet1")

set.seed(alarmseed)
sanets <- cnSearchSA(data=alarm, perturbations=NULL, 
			maxParentSet=2, maxComplexity=600,
			parentsPool=NULL, fixedParentsPool=NULL,
			selectMode="BIC", 
			tempStart=10, tempCoolFact=0.95, tempCheckOrders=16, maxIter=256,
			orderShuffles=4, stopDiff=0.0000000001,
			priorSearch=sanets)
sanets
sanet2 <- cnFindBIC(sanets, nrow(alarm))
cnCompare(alarmnet, sanet2)
sanet2@meta <- "sanet2, SA learning, Stage II"
sanet2
cnDot(sanet2, "sanet2")

set.seed(alarmseed)
sacnets <- cnSearchSAcluster(data=alarm, perturbations=NULL, 
			maxParentSet=2, maxComplexity=600,
			parentsPool=NULL, fixedParentsPool=NULL,
			tempStart=100, tempCoolFact=0.9, tempCheckOrders=16, maxIter=256,
			orderShuffles=-4, stopDiff=0.00000001,
			priorSearch=NULL,
			clusterNodes=4)
sacnets
sanet3<- cnFindBIC(sacnets, nrow(alarm))
cnCompare(alarmnet, sanet3)
sanet3@meta <- "sanet3, SA cluster learning"
cnPlot(sanet3, "sanet3")
sanet3

## restore the random generator's kind
RNGkind(kind=saved.rng[1])

as.character(Sys.getenv("R_DOTVIEWER"))

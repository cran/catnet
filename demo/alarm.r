library(catnet)

## load the alarm data sample and ground-thruth network
data(alarm)
data(alarmnet)
alarmnet

eval <- cnEvaluate(object=alarmnet, data=alarm, perturbations=NULL, maxParentSet=3, echo=FALSE)
##eval <- cnSearchOrder(data=alarm, perturbations=NULL, maxParentSet=3, parentSizes=NULL,maxComplexity=cnComplexity(alarmnet), nodeOrder=cnOrder(alarmnet))
eval
##pdf("alarmeval.pdf")
cnPlot(eval)
##dev.off()

bstnet<- cnFind(eval, cnComplexity(alarmnet))
bstnet
cnCompare(alarmnet, bstnet)


cnSetSeed(123)
sanets1 <- cnSearchSA(data = alarm, perturbations = NULL, maxParentSet = 2, 
     maxComplexity = 600, parentsPool = NULL, fixedParents = NULL, 
     selectMode = "BIC", tempStart = 10, tempCoolFact = 1, tempCheckOrders = 100, 
     maxIter = 100, orderShuffles = 0, stopDiff = 1e-10, numThreads=4, priorSearch = NULL)
sanets1

sanet1 <- cnFindBIC(sanets1)
cnCompare(alarmnet, sanet1)

sanet1@meta <- "sanet1, SA stage I"
sanet1
cnPlot(sanet1, "sanet1")

cnSetSeed(456)
sanets <- cnSearchSA(data = alarm, perturbations = NULL, maxParentSet = 2, 
     maxComplexity = 600, parentsPool = NULL, fixedParents = NULL, 
     selectMode = "BIC", tempStart = 1e-2, tempCoolFact = 0.9, 
     tempCheckOrders = 20, maxIter = 100, orderShuffles = 4, stopDiff = 1e-2, 
     numThreads=4, priorSearch = sanets1, echo=TRUE)
sanets

sanet2 <- cnFindBIC(sanets)
cnCompare(alarmnet, sanet2)

sanet2@meta <- "sanet2, SA stage II"
sanet2
cnPlot(sanet2, "sanet2")


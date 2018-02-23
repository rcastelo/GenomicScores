## for GScores class
setGeneric("name", function(x) standardGeneric("name"))
setGeneric("type", function(x) standardGeneric("type"))
setGeneric("scores", function(object, ranges, ...) standardGeneric("scores"))
setGeneric("gscores", function(x, ranges, ...) standardGeneric("gscores"))
setGeneric("qfun", function(object, ...) standardGeneric("qfun"))
setGeneric("dqfun", function(object, ...) standardGeneric("dqfun"))
setGeneric("populations", function(x) standardGeneric("populations"))
setGeneric("defaultPopulation", function(x) standardGeneric("defaultPopulation"))
setGeneric("defaultPopulation<-", function(x, value) standardGeneric("defaultPopulation<-"))
setGeneric("nsites", function(x) standardGeneric("nsites"))

## for MafDb class
setGeneric("mafByOverlaps", function(x, ...) standardGeneric("mafByOverlaps"))
setGeneric("mafById", function(x, ...) standardGeneric("mafById"))

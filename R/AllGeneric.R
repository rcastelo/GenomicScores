## for GScores class
setGeneric("name", function(x) standardGeneric("name"))
setGeneric("type", function(x) standardGeneric("type"))
setGeneric("scores", function(object, ranges, ...) standardGeneric("scores"))
setGeneric("qfun", function(object) standardGeneric("qfun"))
setGeneric("dqfun", function(object) standardGeneric("dqfun"))

## for MafDb class
setGeneric("mafByOverlaps", function(x, ...) standardGeneric("mafByOverlaps"))
setGeneric("mafById", function(x, ...) standardGeneric("mafById"))
setGeneric("populations", function(x) standardGeneric("populations"))

# STATS GBM extension command

#Licensed Materials - Property of IBM
#IBM SPSS Products: Statistics General
#(c) Copyright IBM Corp. 2014
#US Government Users Restricted Rights - Use, duplication or disclosure 
#restricted by GSA ADP Schedule Contract with IBM Corp.

# author=  'jkp, IBM'
# version=  '1.0.3'
# history
# 03-30-2013 original version
# 04-10-2013 ignore errors in monkey patching for older R/gbm versions
# 05-24-2013 work around problem in older R that caused Relative Importance table
#            to have no variable names in it

helptext = 'STATS GBM DISTRIBUTION=GAUSSIAN, LAPLACE, TDIST, BERNOULLI,
HUBERIZED, MULTINOMIAL, ADABOOST, POISSON, COXPH, or QUANTILE
ALPHA=quantile TDF=df
DEPENDENT=depvar INDEPENDENT=independent variable list
INTERACTIONS=interaction level OFFSET=offset variable
MONOTONE = list of monotone specifications
/OPTIONS NTREES=integer CVFOLDS=integer
SHRINKAGE=shrinkage factor for steps MINNODESIZE=minimum node size
BAGFRAC = fraction
TRAINFRAC = fraction
CVSTRAT = TRUE or FALSE
MISSING = EXCLUDE or INCLUDE
/SAVE MODEL=file specification WORKSPACE= CLEAR or RETAIN
/OUTPUT MARGINALPLOTS = variable list MARGINALPLOTCOUNT=number
BOOSTPLOT=YES or NO BOOSTPLOTMETHOD=oob or test or cv
RELIMP = YES or NO
/HELP.

All items are optional except DISTRIBUTION, DEPENDENT,
and INDEPENDENT

DISTRIBUTION specifies the distribution to use.
GAUSSIAN - squared error
LAPLACE - absolute error
TDIST - t distribution, requires TDF specifying degrees of freedom
BERNOULLI - logistic regression with dependent 0 or 1.  
    Variable level should be scale.
HUBERIZED - huberized hinge loss with dependent 0 or 1
MULTINOMIAL - discrete dependent with more than two categories
    Variable level should be categorical (nominal or ordinal)
ADABOOST - AdaBoost exponential loss with dependent 0 or 1
POISSON - dependent is count
COXPH - right censored dependent
QUANTILE - quantile of dependent.  Requires ALPHA specifying the quantile
as a fraction.

DEPENDENT specifies the dependent variable
INDEPENDENT specifies one or more independent variables
OFFSET - specifies an offset variable for the dependent variable in the equation

DEPENDENT specifies the dependent variable.  The nature of the variable
determines which distributions might be appropriate.  Cases with missing
values for the dependent variable are discarded.

INDEPENDENT specifies a list of independent variables.  Missing values
are allowed.

INTERACTIONS specifies what interaction terms are included in the model.
1 = no interactions, 2 = all two-way interactions, etc.  Default is 1.

OFFSET names an offset variable.

MONOTONE specifies monotonicity of the effect for each independent
variable.  If used, one value must be specified for each variable.
Use p for a positive effect, n for negative, and z for no assumption.

NTREES specifies the number of trees or iterations.  Larger numbers
generally give better accuracy but take more time.  Default is 100.

CVFOLDS specifies the number of cross-validation folds.  Default is 0.

SHRINKAGE specifies the learning rate between 0 and 1.  Values close to
zero tend to do better, but they may require more iterations.

MINNODESIZE is the minimum number of cases in a node.  Default is 10.

TRAINFRAC specifies the fraction of the cases used for the training set.
Default is 1, i.e., no holdout set.

BAGFRAC specifies the fraction of the training set to be used to compute
an out-of-bag estimate of the improvement as the number of trees
increases.  Default is .5.

CVSTRAT specifies for cross validation whether it should be stratified
by class.  Default is yes.
MISSING specifies inclusion or exclusion of user missing values.  Default
is exclude.

MODELFILE specifies a file for saving the estimation results.
WORKSPACE specifies whether or not to keep the workspace contents
in memory.  The saved model can be used to make predictions for
new data with the STATS GBMPRED command.

MARGINALPLOTS specifies up to three independent variable whose
marginal effect (after integrating out other variables) is plotted.
MARGINALPLOTCOUNT specifies the number of independent variables to 
plot up to three.  The first n variables in the independent list
are plotted.  You can specify either of these keywords but not
both.

BOOSTPLOT specifies a plot of the error against the number
of iterations.  The calculation is based on the BOOSTPLOTMETHOD.
That can be OOB (out of bag), TEST (holdout sample), or
CV (cross validation).  OOB is always valid but is said
to be less accurate than the other two.  The default is OOB.

RELIMP specifies whether to display and plot the relative
importance of the independent variables.

/HELP displays this text and does nothing else.
'

library(gbm)

# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}

# monkey patching lattice to print
xyplot <- function(...) {
print(lattice::xyplot(...))
return(NULL)
}

levelplot <- function(...) {
print(lattice::levelplot(...))
return(NULL)
}

stripplot <- function(...) {
print(lattice::stripplot(...))
return(NULL)
}

# There is a bug in gbm in that it does not always wrap calls to lattice functions
# in print, and lattice charts do not print by default.  So we monkeypatch the
# relevant lattice functions by wrapping in the functions above to force printing

# With R2.12 and the corresponding gbm package, some of these names do not
# occur, so we ignore any error here.
library(lattice)
tryCatch({
    unlockBinding(sym="xyplot", env=parent.env(environment(gbm)))
    assign("xyplot", xyplot, envir= parent.env(environment(gbm)))
    unlockBinding(sym="levelplot", env=parent.env(environment(gbm)))
    assign("levelplot", levelplot, envir= parent.env(environment(gbm)))
    unlockBinding(sym="stripplot", env=parent.env(environment(gbm)))
    assign("stripplot", stripplot, envir= parent.env(environment(gbm)))
}, error=function(e) {return})

doGbm <- function(distribution, dep, indep, interactions=1, offset=NULL, alpha=NULL,
    tdf=NULL, monotone=NULL, ntrees=100, cvfolds=0, shrinkage=.001,
    minnodesize=10, bagfrac=.5, trainfrac=1., cvstrat=TRUE, marginalplots=NULL,
    marginalplotcount=NULL, modelfile=NULL, missingvalues="exclude", 
    boostplot=FALSE, boostplotmethod="oob",
    relimp=TRUE, workspace="clear") {
    
    setuplocalization("STATS_GBM")		

    tryCatch(library(gbm), error=function(e){
        stop(gtxtf("The R %s package is required but could not be loaded.","gbm"),call.=FALSE)
        }
    )

    indlen = length(indep)
    dist = distribution
    inputdistribution = distribution
    if (distribution == "tdist") {
        if(is.null(tdf))
            stop(gtxt("Degrees of freedom must be specified for the t distribution"), 
                call.=FALSE)
        dist = list(name=distribution, df=tdf)
        distribution = sprintf("%s, d.f.=%d", distribution, tdf)
    } else if (distribution == "quantile") {
        if(is.null(alpha))
            stop(gtxt("Quantile value must be specified for quantile distribution"), 
                call.=FALSE)
        dist = list(name=distribution, alpha=alpha)
        distribution = sprintf("%s, quantile=%.3f", distribution, alpha)
    }
    if (!is.null(monotone)) {
        if (length(monotone) != indlen)
            stop(gtxt("The number of monotone specifications is different from the number of
                independent variables"), call.=FALSE)
        # map monotone syntax values to 1, -1, 0 form
        monotone = unlist(monotone)
        monotoneInts = unlist(c(p=1, n=-1, z=0)[monotone], use.names=FALSE)
    } else {
        monotoneInts = NULL
        monotone = gtxt("Unspecified")
    }
    if (!is.null(marginalplots) && !is.null(marginalplotcount))
        stop(gtxt("Cannot specify both marginalplots list and marginalplotcount"), call.=FALSE)
    # for marginalplotcount, select up to the number of independent variables up to 3
    if (!is.null(marginalplotcount) && marginalplotcount > 0)
        marginalplots = indep[1:min(length(indep), marginalplotcount)]
    if (length(marginalplots) > 3)
        stop(gtxt("No more than three variables can be specified to plot"),
            call.=FALSE)
    if (length(intersect(indep, marginalplots)) != length(marginalplots))
        stop(gtxt("All variables listed for plots must also appear as independent variables"),
            call.=FALSE)
    
    allvars = c(dep, indep)
    model = paste(indep, collapse="+")
    if (!is.null(offset)) {
        model = paste("offset(", offset, ")+", model, collapse="")
        allvars = c(allvars, offset)
    }
    model = paste(dep, model, sep="~")
    keepUserMissing = ifelse(missingvalues == "exclude", FALSE, TRUE)
    dta <- spssdata.GetDataFromSPSS(allvars, missingValueToNA = TRUE, 
        keepUserMissing=keepUserMissing, factorMode = "levels")
    if (distribution == "multinomial" && !is.factor(dta[[1]]))
        stop(gtxt("Dependent variable must be categorical for the multinomial distribution"),
            call.=FALSE)
    if (distribution == "bernoulli" && is.factor(dta[[1]]))
        stop(gtxt("Dependent variable must be scale measurement level for the Bernoulli distribution"),
            call.=FALSE)
    predel = nrow(dta)
    # remove cases where dependent variable is missing
    dta = dta[complete.cases(dta[,1]),]
    casesdel = predel - nrow(dta)
    if (nrow(dta) == 0)
        stop(gtxt("All cases are missing for the dependent variable", call.=FALSE))

    # gbm will complain if class.stratify.cv is supplied for noncategorical model :-)
    # avoiding do.call here in order to avoiding making another copy of the data

    if (inputdistribution %in% c("bernoulli", "multinomial")) {
        res = tryCatch(gbm(distribution=dist, formula=as.formula(model), data=dta, 
            var.monotone=monotoneInts, n.trees = ntrees, interaction.depth=interactions,
            n.minobsinnode=minnodesize, shrinkage=shrinkage,
            bag.fraction=bagfrac, train.fraction = trainfrac,
            cv.folds = cvfolds, keep.data=FALSE, verbose=FALSE, class.stratify.cv = cvstrat),
            error=function(e) {stop(e$message, call.=FALSE)})
    }
    else {
        res = tryCatch(gbm(distribution=dist, formula=as.formula(model), data=dta, 
            var.monotone=monotoneInts, n.trees = ntrees, interaction.depth=interactions,
            n.minobsinnode=minnodesize, shrinkage=shrinkage,
            bag.fraction=bagfrac, train.fraction = trainfrac,
            cv.folds = cvfolds, keep.data=FALSE, verbose=FALSE),
            error=function(e) {stop(e$message, call.=FALSE)})
    }

    summarylabels=c(gtxt("Distribution"),
        gtxt("Dependent Variable"),
        gtxt("Independent Variables"),
        gtxt("Interaction Level"),
        gtxt("Offset"),
        gtxt("Monotone Specification"),
        gtxt("Number of Trees"),
        gtxt("Shrinkage"),
        gtxt("Cross-Validation Folds"),
        gtxt("Minimum Node Size"),
        gtxt("Bag Fraction"),
        gtxt("Training Fraction"),
        gtxt("Cases used for Fitting"),
        gtxt("Cross Validation Stratified by Class"),
        gtxt("User Missing Values"),
        gtxt("Cases Discarded Due to System Missing Data"),
        gtxt("Save Model As File"),
        gtxt("Date Fit")
    )
    if (!(inputdistribution %in% c("bernoulli", "multinomial"))) {
        classstrat = gtxt("NA")
    } else {
        classstrat = ifelse(cvstrat, gtxt("Yes"), gtxt("No"))
    }
    if (is.null(modelfile)) {
        modelfile = gtxt("Not saved")
    }
    summaryvalues = c(distribution,
        dep,
        paste(indep, collapse=" "),
        interactions,
        ifelse(is.null(offset), "--None--", offset),
        paste(monotone, collapse=" "),
        res$n.trees,
        res$shrinkage,
        res$cv.folds,
        res$n.minobsinnode,
        res$bag.fraction,
        res$train.fraction,
        res$nTrain,
        classstrat,
        missingvalues,
        casesdel,
        modelfile,
        as.character(Sys.time())
    )
    names(summaryvalues) = summarylabels
    settingsdf = data.frame(cbind(summaryvalues))
    colnames(settingsdf) = gtxt("Values")
    StartProcedure(gtxt("Generalized Boosted Regression"), "STATSGBM")
    spsspivottable.Display(settingsdf, 
        title = gtxt("Settings"),
        templateName = "GBMSUMMARY",
        outline=gtxt("Summary"),
        caption = gtxt("Results calculated by the R gbm procedure")
    )

    if (relimp) {
        relimpdf = tryCatch(summary.gbm(res, main=gtxt("Variable Relative Importance")), 
            error = function(e) {
                print(gtxt("Unable to plot relative importance"))
                return(tryCatch(summary.gbm(res, plotit=FALSE)))
            }
        )
        relimpnames = relimpdf[,"var"]
        relimpdf = relimpdf['rel.inf']
        row.names(relimpdf) = relimpnames

        names(relimpdf) = gtxt("Relative Influence")
        spsspivottable.Display(relimpdf,
            title = gtxt("Variable Relative Importance"),
            templateName = "GBMRELIMP",
            outline = gtxt("Relative Importance"),
            caption = gtxt("Importance normalized to sum to 100")
        )
    }
    if (length(marginalplots) > 0) {
            title = paste(gtxt("Marginal Effects of Variables"), 
                paste(marginalplots, collapse=" "),collapse="")
        tryCatch(plot(res, i.var=unlist(marginalplots), main=title),
            error = function(e) {
                print(gtxt("Marginal Plots"))
                print(e)})
    }
    bestiter = NULL
    if (boostplot) {
        if (boostplotmethod == "oob")
            boostplotmethod = toupper(boostplotmethod)
        # No access to title for this plot
        bestiter = tryCatch(gbm.perf(res, oobag.curve=TRUE, method=boostplotmethod),
            error = function(e) {print(e)})
        if (is.numeric(bestiter)) {  # could be an error message
            spsspivottable.Display(data.frame("Value"=bestiter, row.names=gtxt("Best Iteration")),
            title=gtxt("Best Number of Iterations"),
            templateName="GBMBESTITER", caption=sprintf(gtxt("Method:%s"), boostplotmethod))
        } else {
            bestiter = NULL
        }
    }
    # save model results for future use in scoring
    modelproperties = list()
    vdict = spssdictionary.GetDictionaryFromSPSS()
    # save dictionary entry (column)Q for dependent variable
    modelproperties["depvar"] = vdict[match(dep, vdict["varName",])]
    modelproperties["offset"] = ifelse(is.null(offset), "<None>", offset)
    modelproperties["missingvalues"] = keepUserMissing
    modelproperties["bestiter"] = bestiter
    if (!is.null(modelfile)) {
        save(res, settingsdf, modelproperties, file=modelfile, precheck=FALSE)
    }
    # clean up workspace, keeping only what is necessary for predictions
    # if workspace == "retain".
    # The variables are put in the global environment so that they will be retained
    if (workspace == "clear") {
        tryCatch(rm(list=ls()), warning = function(e) {return(NULL)})
    } else {
        rm(list = setdiff(ls(), list("res", "settingsdf", "modelproperties")))
        assign("res", res, envir=.GlobalEnv)
        assign("settingsdf", settingsdf, envir=.GlobalEnv)
        assign("modelproperties", modelproperties, envir=.GlobalEnv)
    }
    spsspkg.EndProcedure()
}

# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 

gtxt <- function(...) {
    return(gettext(...,domain="STATS_GBM"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_GBM"))
}

Run <- function(args) {
#Execute the STATS GBM extension command
    
    cmdname = args[[1]]
    args = args[[2]]
    oobj = spsspkg.Syntax(list(
    spsspkg.Template("DISTRIBUTION", subc="", ktype="str", var="distribution",
        vallist=list("gaussian","laplace","tdist","bernoulli","huberized",
            "multinomial", "adaboost", "poisson", "coxph",
            "quantile")),
    spsspkg.Template("DEPENDENT", subc="",  ktype="existingvarlist", var="dep"),
    spsspkg.Template("INDEPENDENT", subc="", ktype="existingvarlist", var="indep", islist=TRUE),
    spsspkg.Template("INTERACTIONS", subc="",  ktype="int", var="interactions", 
        vallist=list(1)),
    spsspkg.Template("OFFSET", subc="", ktype="existingvarlist", var="offset"),
    spsspkg.Template("ALPHA", subc="", ktype="float", var="alpha", vallist=list(0,1.)),
    spsspkg.Template("TDF", subc="", ktype="int", var="tdf", vallist=list(1)),
    spsspkg.Template("MONOTONE", subc="", ktype="str", var="monotone", 
        vallist=list("p", "n", "z"), islist=TRUE),
    spsspkg.Template("NTREES", subc="OPTIONS", ktype="int", var="ntrees", vallist=list(1)),
    spsspkg.Template("CVFOLDS", subc="OPTIONS", ktype="int", var="cvfolds", vallist=list(0)),
    spsspkg.Template("SHRINKAGE", subc="OPTIONS", ktype="float", var="shrinkage", 
        vallist=list(0,1.)),
    spsspkg.Template("MINNODESIZE", subc="OPTIONS", ktype="int", var="minnodesize", 
        vallist=list(1)),
    spsspkg.Template("BAGFRAC", subc="OPTIONS", ktype="float", var="bagfrac",
        vallist=list(0,1.)),
    spsspkg.Template("TRAINFRAC", subc="OPTIONS", ktype="float", var="trainfrac",
        vallist=list(0,1.)),
    spsspkg.Template("CVSTRAT", subc="OPTIONS", ktype="bool", var="cvstrat"),
    spsspkg.Template("MISSING", subc="OPTIONS", ktype="str", var="missingvalues", 
    vallist=list("include","exclude")),
    spsspkg.Template("MODELFILE", subc="SAVE", ktype="literal", var="modelfile"),
    spsspkg.Template("WORKSPACE", subc="SAVE", ktype="str", var="workspace",
        vallist = list("retain", "clear")),
    spsspkg.Template("MARGINALPLOTS", subc="OUTPUT", ktype="existingvarlist", 
        var="marginalplots", islist=TRUE),
    spsspkg.Template("MARGINALPLOTCOUNT", subc="OUTPUT", ktype="int", 
        var="marginalplotcount", vallist=list(0,3)),
    spsspkg.Template("BOOSTPLOT", subc="OUTPUT", ktype="bool", var="boostplot"),
    spsspkg.Template("BOOSTPLOTMETHOD", subc="OUTPUT", ktype="str", var="boostplotmethod",
        vallist=list("oob", "test", "cv")),
    spsspkg.Template("RELIMP", subc="OUTPUT", ktype="bool", var="relimp"),
    spsspkg.Template("HELP", subc="", ktype="bool")
    ))


# A HELP subcommand overrides all else
    if ("HELP" %in% attr(args,"names")) {
        #writeLines(helptext)
        helper(cmdname)
    }
    else {
        res <- spsspkg.processcmd(oobj, args, "doGbm")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}
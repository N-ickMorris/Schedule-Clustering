# -----------------------------------------------------------------------------------
# ---- Set Up -----------------------------------------------------------------------
# -----------------------------------------------------------------------------------

# set the path of where the input files are
mywd = "C:/ ... /Solution-Clusters"

# open up a graphics window
windows()

# -----------------------------------------------------------------------------------
# ---- Packages ---------------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# data handling
require(data.table)
require(tm)
require(gtools)
require(stringdist)
require(zoo)
require(missRanger)
require(TTR)
require(Matrix)

# plotting
require(ggplot2)
require(gridExtra)
require(GGally)
require(ggpairs)
require(scales)
require(scatterplot3d)
require(gridGraphics)
require(corrplot)
require(VIM)
require(plot3D)
require(grDevices)

# modeling
require(fitdistrplus)
require(bestNormalize)
require(fpc)
require(caret)
# require(ranger)
require(cluster)
require(car)
# require(nortest)
# require(neuralnet)
require(h2o)
require(MLmetrics)

# parallel computing
require(foreach)
require(doSNOW)

# choose how many workers to use when parallel computing
workers = max(1, floor((2/3) * detectCores()))
workers = max(1, detectCores() - 2)

}

# -----------------------------------------------------------------------------------
# ---- Functions --------------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# ---- prints the data types of each column in a data frame -------------------------

types = function(dat)
{
  require(data.table)
  
  # make dat into a data.table
  dat = data.table(dat)
  
  # get the column names
  column = names(dat)
  
  # get the class of the columns
  dataType = sapply(1:ncol(dat), function(i) class(unlist(dat[, i, with = FALSE])))
  
  # compute the number of levels for each column
  levels = sapply(1:ncol(dat), function(i) ifelse(dataType[i] == "factor", length(levels(droplevels(unlist(dat[, i, with = FALSE])))), 0))
  
  # compute the number of unique values for each column
  uniqueValues = sapply(1:ncol(dat), function(i) length(unique(unname(unlist(dat[, i, with = FALSE])))))
  
  # compute the portion of missing data
  missing = sapply(1:ncol(dat), function(i) nrow(na.omit(dat[, i, with = FALSE], invert = TRUE)) / nrow(dat))
  
  # build the output table 
  output = data.table(column, id = 1:length(column), dataType, levels, uniqueValues, missing)
  
  # order output by dataType
  output = output[order(dataType)]
  
  return(output)
}

# ---- converts all columns to a character data type --------------------------------

tochar = function(dat)
{
  require(data.table)
  
  # make dat into a data.frame
  dat = data.table(dat)
  
  # get the column names
  column = names(dat)
  
  # get the values in the columns and convert them to character data types
  values = lapply(1:ncol(dat), function(i) as.character(unname(unlist(dat[, i, with = FALSE]))))
  
  # combine the values back into a data.frame
  dat = data.table(do.call("cbind", values), stringsAsFactors = FALSE)
  
  # give dat its column names
  setnames(dat, column)
  
  return(dat)
}

# ---- a qualitative color scheme ---------------------------------------------------

qcolor = function(n, a = 1)
{
  require(grDevices)
  require(scales)
  return(alpha(colorRampPalette(c("#e41a1c", "#0099ff", "#4daf4a", "#984ea3", "#ff7f00", "#ff96ca", "#a65628"))(n), 
               a))
}

# ---- the ggplot2 color scheme -----------------------------------------------------

ggcolor = function(n, a = 1)
{
  require(grDevices)
  require(scales)
  return(alpha(hcl(h = seq(15, 375, length = n + 1), 
                   l = 65, c = 100)[1:n], 
               a))
}

# ---- prints out a dat file object in ampl syntax ----------------------------------

ampl = function(dat, object = "param", name = "c")
{
  tochar = function(dat)
  {
    require(data.table)
    
    # make dat into a data.frame
    dat = data.table(dat)
    
    # get the column names
    column = names(dat)
    
    # get the values in the columns and convert them to character data types
    values = lapply(1:ncol(dat), function(i) as.character(unname(unlist(dat[, i, with = FALSE]))))
    
    # combine the values back into a data.frame
    dat = data.table(do.call("cbind", values), stringsAsFactors = FALSE)
    
    # give dat its column names
    setnames(dat, column)
    
    return(dat)
  }
  
  # make sure the data is a data frame object
  dat = tochar(dat)
  
  # every parameter/set object in an ampl dat file must end with a semicolon
  # so set up 1 semicolon to give to dat
  semicolon = c(";", rep(" ", ncol(dat) - 1))
  
  # add this semicolon as the last row of the data frame
  result = data.frame(rbind(dat, semicolon))
  
  # every parameter/set object in an ample dat file must begin with the name of the object and what it equals
  # for example: param c := 
  # so set up a header to give to dat
  header = c(paste(object, name, ":="), rep(" ", ncol(dat) - 1))
  
  # update the column names of dat to be the header we created
  colnames(result) = header
  
  # print out the result without any row names
  # print out the result left adjusted
  # print(result, right = FALSE, row.names = FALSE)
  
  return(result)	
}

# ---- compares the quantiles of emprical data against the quantiles of any statistical distribution 

ggqq = function(x, distribution = "norm", ..., conf = 0.95, probs = c(0.25, 0.75), alpha = 0.33, basefont = 20, main = "", xlab = "\nTheoretical Quantiles", ylab = "Empirical Quantiles\n")
{
  require(ggplot2)
  
  # compute the sample quantiles and theoretical quantiles
  q.function = eval(parse(text = paste0("q", distribution)))
  d.function = eval(parse(text = paste0("d", distribution)))
  x = na.omit(x)
  ord = order(x)
  n = length(x)
  P = ppoints(length(x))
  df = data.frame(ord.x = x[ord], z = q.function(P, ...))
  
  # compute the quantile line
  Q.x = quantile(df$ord.x, c(probs[1], probs[2]))
  Q.z = q.function(c(probs[1], probs[2]), ...)
  b = diff(Q.x) / diff(Q.z)
  coef = c(Q.x[1] - (b * Q.z[1]), b)
  
  # compute the confidence interval band
  zz = qnorm(1 - (1 - conf) / 2)
  SE = (coef[2] / d.function(df$z, ...)) * sqrt(P * (1 - P) / n)
  fit.value = coef[1] + (coef[2] * df$z)
  df$upper = fit.value + (zz * SE)
  df$lower = fit.value - (zz * SE)
  
  # plot the qqplot
  p = ggplot(df, aes(x = z, y = ord.x)) + 
    geom_point(color = "blue", alpha = alpha) +
    geom_abline(intercept = coef[1], slope = coef[2], size = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) +
    coord_cartesian(ylim = c(min(df$ord.x), max(df$ord.x))) + 
    labs(x = xlab, y = ylab) +
    theme_bw(base_size = basefont) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # conditional additions
  if(main != "")(p = p + ggtitle(main))
  
  return(p)
}

# ---- plots 4 residual plots -------------------------------------------------------

residplots = function(actual, fit, binwidth = NULL, from = NULL, to = NULL, by = NULL, histlabel.y = -10, basefont = 20)
{
  require(ggplot2)
  
  residual = actual - fit 
  DF = data.frame("actual" = actual, "fit" = fit, "residual" = residual)
  
  rvfPlot = ggplot(DF, aes(x = fit, y = residual)) + 
    geom_point(na.rm = TRUE) +
    stat_smooth(method = "loess", se = FALSE, na.rm = TRUE, color = "blue") +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    xlab("Fitted values") +
    ylab("Residuals") +
    ggtitle("Residual vs Fitted Plot") + 
    theme_bw(base_size = basefont) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  ggqq = function(x, distribution = "norm", ..., conf = 0.95, probs = c(0.25, 0.75), note = TRUE, alpha = 0.33, main = "", xlab = "\nTheoretical Quantiles", ylab = "Empirical Quantiles\n")
  {
    # compute the sample quantiles and theoretical quantiles
    q.function = eval(parse(text = paste0("q", distribution)))
    d.function = eval(parse(text = paste0("d", distribution)))
    x = na.omit(x)
    ord = order(x)
    n = length(x)
    P = ppoints(length(x))
    DF = data.frame(ord.x = x[ord], z = q.function(P, ...))
    
    # compute the quantile line
    Q.x = quantile(DF$ord.x, c(probs[1], probs[2]))
    Q.z = q.function(c(probs[1], probs[2]), ...)
    b = diff(Q.x) / diff(Q.z)
    coef = c(Q.x[1] - (b * Q.z[1]), b)
    
    # compute the confidence interval band
    zz = qnorm(1 - (1 - conf) / 2)
    SE = (coef[2] / d.function(DF$z, ...)) * sqrt(P * (1 - P) / n)
    fit.value = coef[1] + (coef[2] * DF$z)
    DF$upper = fit.value + (zz * SE)
    DF$lower = fit.value - (zz * SE)
    
    # plot the qqplot
    p = ggplot(DF, aes(x = z, y = ord.x)) + 
      geom_point(color = "black", alpha = alpha) +
      geom_abline(intercept = coef[1], slope = coef[2], size = 1, color = "blue") +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15) +
      coord_cartesian(ylim = c(min(DF$ord.x), max(DF$ord.x))) + 
      labs(x = xlab, y = ylab) +
      theme_bw(base_size = basefont) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    # conditional additions
    if(main != "")(p = p + ggtitle(main))
    
    return(p)
  }
  
  qqPlot = ggqq(residual, 
                alpha = 1,				  
                main = "Normal Q-Q Plot", 
                xlab = "Theoretical Quantiles", 
                ylab = "Residuals")
  
  rvtPlot = ggplot(data.frame("x" = 1:length(DF$residual), "y" = DF$residual), aes(x = x, y = y)) + 
    geom_line(na.rm = TRUE) +
    geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
    xlab("Obs. Number") +
    ylab("Residuals") +
    ggtitle("Residual Time Series") + 
    theme_bw(base_size = basefont) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  test = t.test(DF$residual)
  
  CI = data.frame("x" = test$estimate, 
                  "LCB" = test$conf.int[1], 
                  "UCB" = test$conf.int[2], 
                  row.names = 1)
  
  histPlot = ggplot(DF, aes(x = residual)) +
    geom_histogram(color = "white", fill = "black", binwidth = binwidth) +
    geom_segment(data = CI, aes(x = LCB, xend = LCB, y = 0, yend = Inf), color = "blue") +
    geom_segment(data = CI, aes(x = UCB, xend = UCB, y = 0, yend = Inf), color = "blue") +
    annotate("text", x = CI$x, y = histlabel.y, 
             label = "T-Test C.I.", size = 5, 
             color = "blue", fontface = 2) + 
    ggtitle("Residual Histogram") +
    labs(x = "Residuals", y = "Frequency") +
    theme_bw(base_size = basefont) +
    theme(legend.key.size = unit(.25, "in"),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  if(class(from) != "NULL" & class(to) != "NULL" & class(by) != "NULL") (histPlot = histPlot + scale_x_continuous(breaks = seq(from = from, to = to, by = by)))
  
  return(list("rvfPlot" = rvfPlot, 
              "qqPlot" = qqPlot, 
              "rvtPlot" = rvtPlot,  
              "histPlot" = histPlot))
}

# ---- builds a square confusion matrix ---------------------------------------------

confusion = function(ytrue, ypred)
{
  require(gtools)
  
  # make predicted and actual vectors into factors, if they aren't already
  if(class(ytrue) != "factor") ytrue = factor(ytrue)
  if(class(ypred) != "factor") ypred = factor(ypred)
  
  # combine their levels into one unique set of levels
  common.levels = mixedsort(unique(c(levels(ytrue), levels(ypred))))
  
  # give each vector the same levels
  ytrue = factor(ytrue, levels = common.levels)
  ypred = factor(ypred, levels = common.levels)
  
  # build the confusion matrix
  output = table("Actual" = ytrue, "Predicted" = ypred)
  
  # return the confusion matrix
  return(output)
}

# ---- generates a logarithmically spaced sequence ----------------------------------

lseq = function(from, to, length.out)
{
  return(exp(seq(log(from), log(to), length.out = length.out)))
}

}

# -----------------------------------------------------------------------------------
# ---- Prepare the Data -------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# set the work directory
setwd(mywd)

# import the data
sol = fread("Solution Data - Long.csv")
risk = fread("Solution Data with Risk - Short.csv")
ant = fread("bundle-antigen-mapping-52B.csv")

# get market level data from risk
mark = data.table(risk[,.(Markets, MarketID, Country_Risk, GNIpc, Births, Birth_Mortality)])

# aggregate Country_Risk, GNIpc, Births, Birth_Mortality by market
mark = mark[,.(Country_Risk = mean(Country_Risk, na.rm = TRUE), 
               GNIpc = mean(GNIpc, na.rm = TRUE), 
               Birth_Mortality = mean(Birth_Mortality, na.rm = TRUE)), 
            by = .(Markets, MarketID)]

# get ABP output from risk
tss = data.table(risk[,.(fileID, TSS, TCS, TPF, CMV, PMV)])

# remove duplicates from tss
tss = tss[!duplicated(tss)]

# give sol an id column to preserve row order
sol[, id := 1:nrow(sol)]

# join mark onto sol
setkey(mark, Markets, MarketID)
setkey(sol, Markets, MarketID)
sol = mark[sol]

# join tss onto sol
setkey(tss, fileID)
setkey(sol, fileID)
sol = tss[sol]

# order sol by id and then remove it
sol = sol[order(id)]
sol[, id := NULL]

# create a SolutionID variable
sol[, SolutionID := paste(fileID, MarketID, sep = "_")]
sol[, SolutionID := as.numeric(factor(SolutionID, levels = unique(SolutionID)))]

# create a subset of sol to prepare solution schedules
dat = data.table(sol[,.(SolutionID, Bundle, Selling_Qty, Selling_Price_Low, Selling_Price_High)])

# convert Selling_Qty, Selling_Price_Low, Selling_Price_High into a single column of values
dat = melt(dat, id.vars = c("SolutionID", "Bundle"))

# update variable
dat[, variable := paste(variable, Bundle, sep = "_")]
dat[, variable := factor(variable, levels = unique(variable))]

# remove Bundle from dat
dat[, Bundle := NULL]

# convert dat into wide format where the rows are SolutionID, the columns are variable, and the values in the matrix are value
dat = dcast(dat, formula = SolutionID ~ variable, value.var = "value")

# remove SolutionID as a column from dat
dat[, SolutionID := NULL]

# remove unneeded objects
rm(mark, risk, tss)

# free up RAM
gc()

}

# -----------------------------------------------------------------------------------
# ---- Model the Data ---------------------------------------------------------------
# -----------------------------------------------------------------------------------

{

# are the models already built?
models.built = TRUE

if(models.built)
{
  # choose the number of workers for parallel processing
  workers = getDTthreads() - 2
  
  # initialize the h2o instance
  h2o.init(nthreads = workers, max_mem_size = "8g")
  
  # remove any objects in the h2o instance
  h2o.removeAll()
  
  # remove the progress bar when model building
  h2o.no_progress()
  
  # read in the data
  train = fread("train.csv")
  valid = fread("valid.csv")
  test = fread("test.csv")
  
  # load the models
  glm.mod = h2o.loadModel(path = paste0(getwd(), "/glm.mod"))
  rf.mod = h2o.loadModel(path = paste0(getwd(), "/rf.mod"))
  gb.mod = h2o.loadModel(path = paste0(getwd(), "/gb.mod"))
  nnet.mod = h2o.loadModel(path = paste0(getwd(), "/nnet.mod"))
  stack.mod = h2o.loadModel(path = paste0(getwd(), "/stack.mod"))
  
  # rename the columns of train, valid, and test
  feature.names = glm.mod@model$coefficients_table$names[-1]
  setnames(train, c("Cluster", feature.names))
  setnames(valid, c("Cluster", feature.names))
  setnames(test, c("Cluster", feature.names))
  
  # convert train, valid, and test into h2o objects
  train.YX.h2o = as.h2o(train)
  valid.YX.h2o = as.h2o(valid)
  test.YX.h2o = as.h2o(test)
  
  # ---------------
  # ---- GLM ------
  # ---------------
  
  {
    # make predictions on each data set
    ynew.train = as.matrix(predict(glm.mod, newdata = train.YX.h2o)[,-1])
    ynew.valid = as.matrix(predict(glm.mod, newdata = valid.YX.h2o)[,-1])
    ynew.test = as.matrix(predict(glm.mod, newdata = test.YX.h2o)[,-1])
    
    # get the true values from each data set
    ytrue.train = as.data.frame(train.YX.h2o[,1])
    ytrue.valid = as.data.frame(valid.YX.h2o[,1])
    ytrue.test = as.data.frame(test.YX.h2o[,1])
    
    # convert the true values to factor data types
    target.levels = gsub("p", "", colnames(ynew.train))
    ytrue.train[,1] = factor(ytrue.train[,1], levels = target.levels)
    ytrue.valid[,1] = factor(ytrue.valid[,1], levels = target.levels)
    ytrue.test[,1] = factor(ytrue.test[,1], levels = target.levels)
    
    # ---- compute multi log loss ----
    
    # build a matrix indicating the true class values for each data set
    ytrue.train.mat = model.matrix(~., data = ytrue.train)[,-1]
    ytrue.train.mat = cbind(1 - rowSums(ytrue.train.mat), ytrue.train.mat)
    
    ytrue.valid.mat = model.matrix(~., data = ytrue.valid)[,-1]
    ytrue.valid.mat = cbind(1 - rowSums(ytrue.valid.mat), ytrue.valid.mat)
    
    ytrue.test.mat = model.matrix(~., data = ytrue.test)[,-1]
    ytrue.test.mat = cbind(1 - rowSums(ytrue.test.mat), ytrue.test.mat)
    
    # compute the multi-class logarithmic loss for each data set
    mll.train = MultiLogLoss(y_pred = ynew.train, y_true = ytrue.train.mat)
    mll.valid = MultiLogLoss(y_pred = ynew.valid, y_true = ytrue.valid.mat)
    mll.test = MultiLogLoss(y_pred = ynew.test, y_true = ytrue.test.mat)
    
    # free up RAM
    gc()
    
    # ---- compute kappa ----
    
    # get the predicted classes and actual classes for each data set
    ynew.train.code = apply(ynew.train, 1, which.max) - 1
    ytrue.train.code = factor(as.numeric(ytrue.train[,1]) - 1, levels = 0:(ncol(ytrue.train.mat) - 1))
    
    ynew.valid.code = apply(ynew.valid, 1, which.max) - 1
    ytrue.valid.code = factor(as.numeric(ytrue.valid[,1]) - 1, levels = 0:(ncol(ytrue.valid.mat) - 1))
    
    ynew.test.code = apply(ynew.test, 1, which.max) - 1
    ytrue.test.code = factor(as.numeric(ytrue.test[,1]) - 1, levels = 0:(ncol(ytrue.test.mat) - 1))
    
    # build a square confusion matrix for each data set
    conf.train = confusion(ytrue = ytrue.train.code, ypred = ynew.train.code)
    conf.valid = confusion(ytrue = ytrue.valid.code, ypred = ynew.valid.code)
    conf.test = confusion(ytrue = ytrue.test.code, ypred = ynew.test.code)
    
    # get the total number of observations for each data set
    n.train = sum(conf.train)
    n.valid = sum(conf.valid)
    n.test = sum(conf.test)
    
    # get the vector of correct predictions for each data set 
    dia.train = diag(conf.train)
    dia.valid = diag(conf.valid)
    dia.test = diag(conf.test)
    
    # get the vector of the number of observations per class for each data set
    rsum.train = rowSums(conf.train)
    rsum.valid = rowSums(conf.valid)
    rsum.test = rowSums(conf.test)
    
    # get the vector of the number of predictions per class for each data set
    csum.train = colSums(conf.train)
    csum.valid = colSums(conf.valid)
    csum.test = colSums(conf.test)
    
    # get the proportion of observations per class for each data set
    p.train = rsum.train / n.train
    p.valid = rsum.valid / n.valid
    p.test = rsum.test / n.test
    
    # get the proportion of predcitions per class for each data set
    q.train = csum.train / n.train
    q.valid = csum.valid / n.valid
    q.test = csum.test / n.test
    
    # compute accuracy for each data set
    acc.train = sum(dia.train) / n.train
    acc.valid = sum(dia.valid) / n.valid
    acc.test = sum(dia.test) / n.test
    
    # compute expected accuracy for each data set
    exp.acc.train = sum(p.train * q.train)
    exp.acc.valid = sum(p.valid * q.valid)
    exp.acc.test = sum(p.test * q.test)
    
    # compute kappa for each data set
    kap.train = (acc.train - exp.acc.train) / (1 - exp.acc.train)
    kap.valid = (acc.valid - exp.acc.valid) / (1 - exp.acc.valid)
    kap.test = (acc.test - exp.acc.test) / (1 - exp.acc.test)
    
    # ---- compute one-vs-all metrics ----
    
    # compute a binary confusion matrix for each class, for each data set
    one.v.all.train = lapply(1:nrow(conf.train), function(i)
    {
      # get the four entries of a binary confusion matrix
      v = c(conf.train[i,i], 
            rsum.train[i] - conf.train[i,i], 
            csum.train[i] - conf.train[i,i], 
            n.train - rsum.train[i] - csum.train[i] + conf.train[i,i]);
      
      # build the confusion matrix
      return(matrix(v, nrow = 2, byrow = TRUE))
    })
    
    one.v.all.valid = lapply(1:nrow(conf.valid), function(i)
    {
      # get the four entries of a binary confusion matrix
      v = c(conf.valid[i,i], 
            rsum.valid[i] - conf.valid[i,i], 
            csum.valid[i] - conf.valid[i,i], 
            n.valid - rsum.valid[i] - csum.valid[i] + conf.valid[i,i]);
      
      # build the confusion matrix
      return(matrix(v, nrow = 2, byrow = TRUE))
    })
    
    one.v.all.test = lapply(1:nrow(conf.test), function(i)
    {
      # get the four entries of a binary confusion matrix
      v = c(conf.test[i,i], 
            rsum.test[i] - conf.test[i,i], 
            csum.test[i] - conf.test[i,i], 
            n.test - rsum.test[i] - csum.test[i] + conf.test[i,i]);
      
      # build the confusion matrix
      return(matrix(v, nrow = 2, byrow = TRUE))
    })
    
    # sum up all of the matrices for each data set
    one.v.all.train = Reduce('+', one.v.all.train)
    one.v.all.valid = Reduce('+', one.v.all.valid)
    one.v.all.test = Reduce('+', one.v.all.test)
    
    # compute the micro average accuracy for each data set
    micro.acc.train = sum(diag(one.v.all.train)) / sum(one.v.all.train)
    micro.acc.valid = sum(diag(one.v.all.valid)) / sum(one.v.all.valid)
    micro.acc.test = sum(diag(one.v.all.test)) / sum(one.v.all.test)
    
    # get the macro accuracy for each data set
    macro.acc.train = acc.train
    macro.acc.valid = acc.valid
    macro.acc.test = acc.test
    
    # ---- finalize output ----
    
    # build a final metrics table for glm.mod
    glm.mod.table = data.table(Model = "Regression",
                                 Metric = c("Log_Loss", "Kappa", "Macro_Accuracy", "Micro_Accuracy"),
                                 Train = c(mll.train, kap.train, macro.acc.train, micro.acc.train),
                                 Valid = c(mll.valid, kap.valid, macro.acc.valid, micro.acc.valid),
                                 Test = c(mll.test, kap.test, macro.acc.test, micro.acc.test))
    
    # update the row and column names of the confusion matrices
    colnames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
    rownames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
    
    colnames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
    rownames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
    
    colnames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
    rownames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
    
    # build a list of confusion matrices
    glm.confusion = list("Train" = conf.train, 
                           "Valid" = conf.valid, 
                           "Test" = conf.test)
    
    # evaluate model performance
    glm.mod.table
    glm.confusion
    
    # free up RAM
    gc()
  }
  
  # ---------------
  # ---- RF -------
  # ---------------
  
  {
    # make predictions on each data set
    ynew.train = as.matrix(predict(rf.mod, newdata = train.YX.h2o)[,-1])
    ynew.valid = as.matrix(predict(rf.mod, newdata = valid.YX.h2o)[,-1])
    ynew.test = as.matrix(predict(rf.mod, newdata = test.YX.h2o)[,-1])
    
    # get the true values from each data set
    ytrue.train = as.data.frame(train.YX.h2o[,1])
    ytrue.valid = as.data.frame(valid.YX.h2o[,1])
    ytrue.test = as.data.frame(test.YX.h2o[,1])
    
    # convert the true values to factor data types
    target.levels = gsub("p", "", colnames(ynew.train))
    ytrue.train[,1] = factor(ytrue.train[,1], levels = target.levels)
    ytrue.valid[,1] = factor(ytrue.valid[,1], levels = target.levels)
    ytrue.test[,1] = factor(ytrue.test[,1], levels = target.levels)
    
    # ---- compute multi log loss ----
    
    # build a matrix indicating the true class values for each data set
    ytrue.train.mat = model.matrix(~., data = ytrue.train)[,-1]
    ytrue.train.mat = cbind(1 - rowSums(ytrue.train.mat), ytrue.train.mat)
    
    ytrue.valid.mat = model.matrix(~., data = ytrue.valid)[,-1]
    ytrue.valid.mat = cbind(1 - rowSums(ytrue.valid.mat), ytrue.valid.mat)
    
    ytrue.test.mat = model.matrix(~., data = ytrue.test)[,-1]
    ytrue.test.mat = cbind(1 - rowSums(ytrue.test.mat), ytrue.test.mat)
    
    # compute the multi-class logarithmic loss for each data set
    mll.train = MultiLogLoss(y_pred = ynew.train, y_true = ytrue.train.mat)
    mll.valid = MultiLogLoss(y_pred = ynew.valid, y_true = ytrue.valid.mat)
    mll.test = MultiLogLoss(y_pred = ynew.test, y_true = ytrue.test.mat)
    
    # free up RAM
    gc()
    
    # ---- compute kappa ----
    
    # get the predicted classes and actual classes for each data set
    ynew.train.code = apply(ynew.train, 1, which.max) - 1
    ytrue.train.code = factor(as.numeric(ytrue.train[,1]) - 1, levels = 0:(ncol(ytrue.train.mat) - 1))
    
    ynew.valid.code = apply(ynew.valid, 1, which.max) - 1
    ytrue.valid.code = factor(as.numeric(ytrue.valid[,1]) - 1, levels = 0:(ncol(ytrue.valid.mat) - 1))
    
    ynew.test.code = apply(ynew.test, 1, which.max) - 1
    ytrue.test.code = factor(as.numeric(ytrue.test[,1]) - 1, levels = 0:(ncol(ytrue.test.mat) - 1))
    
    # build a square confusion matrix for each data set
    conf.train = confusion(ytrue = ytrue.train.code, ypred = ynew.train.code)
    conf.valid = confusion(ytrue = ytrue.valid.code, ypred = ynew.valid.code)
    conf.test = confusion(ytrue = ytrue.test.code, ypred = ynew.test.code)
    
    # get the total number of observations for each data set
    n.train = sum(conf.train)
    n.valid = sum(conf.valid)
    n.test = sum(conf.test)
    
    # get the vector of correct predictions for each data set 
    dia.train = diag(conf.train)
    dia.valid = diag(conf.valid)
    dia.test = diag(conf.test)
    
    # get the vector of the number of observations per class for each data set
    rsum.train = rowSums(conf.train)
    rsum.valid = rowSums(conf.valid)
    rsum.test = rowSums(conf.test)
    
    # get the vector of the number of predictions per class for each data set
    csum.train = colSums(conf.train)
    csum.valid = colSums(conf.valid)
    csum.test = colSums(conf.test)
    
    # get the proportion of observations per class for each data set
    p.train = rsum.train / n.train
    p.valid = rsum.valid / n.valid
    p.test = rsum.test / n.test
    
    # get the proportion of predcitions per class for each data set
    q.train = csum.train / n.train
    q.valid = csum.valid / n.valid
    q.test = csum.test / n.test
    
    # compute accuracy for each data set
    acc.train = sum(dia.train) / n.train
    acc.valid = sum(dia.valid) / n.valid
    acc.test = sum(dia.test) / n.test
    
    # compute expected accuracy for each data set
    exp.acc.train = sum(p.train * q.train)
    exp.acc.valid = sum(p.valid * q.valid)
    exp.acc.test = sum(p.test * q.test)
    
    # compute kappa for each data set
    kap.train = (acc.train - exp.acc.train) / (1 - exp.acc.train)
    kap.valid = (acc.valid - exp.acc.valid) / (1 - exp.acc.valid)
    kap.test = (acc.test - exp.acc.test) / (1 - exp.acc.test)
    
    # ---- compute one-vs-all metrics ----
    
    # compute a binary confusion matrix for each class, for each data set
    one.v.all.train = lapply(1:nrow(conf.train), function(i)
    {
      # get the four entries of a binary confusion matrix
      v = c(conf.train[i,i], 
            rsum.train[i] - conf.train[i,i], 
            csum.train[i] - conf.train[i,i], 
            n.train - rsum.train[i] - csum.train[i] + conf.train[i,i]);
      
      # build the confusion matrix
      return(matrix(v, nrow = 2, byrow = TRUE))
    })
    
    one.v.all.valid = lapply(1:nrow(conf.valid), function(i)
    {
      # get the four entries of a binary confusion matrix
      v = c(conf.valid[i,i], 
            rsum.valid[i] - conf.valid[i,i], 
            csum.valid[i] - conf.valid[i,i], 
            n.valid - rsum.valid[i] - csum.valid[i] + conf.valid[i,i]);
      
      # build the confusion matrix
      return(matrix(v, nrow = 2, byrow = TRUE))
    })
    
    one.v.all.test = lapply(1:nrow(conf.test), function(i)
    {
      # get the four entries of a binary confusion matrix
      v = c(conf.test[i,i], 
            rsum.test[i] - conf.test[i,i], 
            csum.test[i] - conf.test[i,i], 
            n.test - rsum.test[i] - csum.test[i] + conf.test[i,i]);
      
      # build the confusion matrix
      return(matrix(v, nrow = 2, byrow = TRUE))
    })
    
    # sum up all of the matrices for each data set
    one.v.all.train = Reduce('+', one.v.all.train)
    one.v.all.valid = Reduce('+', one.v.all.valid)
    one.v.all.test = Reduce('+', one.v.all.test)
    
    # compute the micro average accuracy for each data set
    micro.acc.train = sum(diag(one.v.all.train)) / sum(one.v.all.train)
    micro.acc.valid = sum(diag(one.v.all.valid)) / sum(one.v.all.valid)
    micro.acc.test = sum(diag(one.v.all.test)) / sum(one.v.all.test)
    
    # get the macro accuracy for each data set
    macro.acc.train = acc.train
    macro.acc.valid = acc.valid
    macro.acc.test = acc.test
    
    # ---- finalize output ----
    
    # build a final metrics table for rf.mod
    rf.mod.table = data.table(Model = "Random Forest",
                                 Metric = c("Log_Loss", "Kappa", "Macro_Accuracy", "Micro_Accuracy"),
                                 Train = c(mll.train, kap.train, macro.acc.train, micro.acc.train),
                                 Valid = c(mll.valid, kap.valid, macro.acc.valid, micro.acc.valid),
                                 Test = c(mll.test, kap.test, macro.acc.test, micro.acc.test))
    
    # update the row and column names of the confusion matrices
    colnames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
    rownames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
    
    colnames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
    rownames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
    
    colnames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
    rownames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
    
    # build a list of confusion matrices
    rf.confusion = list("Train" = conf.train, 
                           "Valid" = conf.valid, 
                           "Test" = conf.test)
    
    # evaluate model performance
    rf.mod.table
    rf.confusion
    
    # free up RAM
    gc()
  }
  
  # ---------------
  # ---- GB -------
  # ---------------
  
  {
    # make predictions on each data set
    ynew.train = as.matrix(predict(gb.mod, newdata = train.YX.h2o)[,-1])
    ynew.valid = as.matrix(predict(gb.mod, newdata = valid.YX.h2o)[,-1])
    ynew.test = as.matrix(predict(gb.mod, newdata = test.YX.h2o)[,-1])
    
    # get the true values from each data set
    ytrue.train = as.data.frame(train.YX.h2o[,1])
    ytrue.valid = as.data.frame(valid.YX.h2o[,1])
    ytrue.test = as.data.frame(test.YX.h2o[,1])
    
    # convert the true values to factor data types
    target.levels = gsub("p", "", colnames(ynew.train))
    ytrue.train[,1] = factor(ytrue.train[,1], levels = target.levels)
    ytrue.valid[,1] = factor(ytrue.valid[,1], levels = target.levels)
    ytrue.test[,1] = factor(ytrue.test[,1], levels = target.levels)
    
    # ---- compute multi log loss ----
    
    # build a matrix indicating the true class values for each data set
    ytrue.train.mat = model.matrix(~., data = ytrue.train)[,-1]
    ytrue.train.mat = cbind(1 - rowSums(ytrue.train.mat), ytrue.train.mat)
    
    ytrue.valid.mat = model.matrix(~., data = ytrue.valid)[,-1]
    ytrue.valid.mat = cbind(1 - rowSums(ytrue.valid.mat), ytrue.valid.mat)
    
    ytrue.test.mat = model.matrix(~., data = ytrue.test)[,-1]
    ytrue.test.mat = cbind(1 - rowSums(ytrue.test.mat), ytrue.test.mat)
    
    # compute the multi-class logarithmic loss for each data set
    mll.train = MultiLogLoss(y_pred = ynew.train, y_true = ytrue.train.mat)
    mll.valid = MultiLogLoss(y_pred = ynew.valid, y_true = ytrue.valid.mat)
    mll.test = MultiLogLoss(y_pred = ynew.test, y_true = ytrue.test.mat)
    
    # free up RAM
    gc()
    
    # ---- compute kappa ----
    
    # get the predicted classes and actual classes for each data set
    ynew.train.code = apply(ynew.train, 1, which.max) - 1
    ytrue.train.code = factor(as.numeric(ytrue.train[,1]) - 1, levels = 0:(ncol(ytrue.train.mat) - 1))
    
    ynew.valid.code = apply(ynew.valid, 1, which.max) - 1
    ytrue.valid.code = factor(as.numeric(ytrue.valid[,1]) - 1, levels = 0:(ncol(ytrue.valid.mat) - 1))
    
    ynew.test.code = apply(ynew.test, 1, which.max) - 1
    ytrue.test.code = factor(as.numeric(ytrue.test[,1]) - 1, levels = 0:(ncol(ytrue.test.mat) - 1))
    
    # build a square confusion matrix for each data set
    conf.train = confusion(ytrue = ytrue.train.code, ypred = ynew.train.code)
    conf.valid = confusion(ytrue = ytrue.valid.code, ypred = ynew.valid.code)
    conf.test = confusion(ytrue = ytrue.test.code, ypred = ynew.test.code)
    
    # get the total number of observations for each data set
    n.train = sum(conf.train)
    n.valid = sum(conf.valid)
    n.test = sum(conf.test)
    
    # get the vector of correct predictions for each data set 
    dia.train = diag(conf.train)
    dia.valid = diag(conf.valid)
    dia.test = diag(conf.test)
    
    # get the vector of the number of observations per class for each data set
    rsum.train = rowSums(conf.train)
    rsum.valid = rowSums(conf.valid)
    rsum.test = rowSums(conf.test)
    
    # get the vector of the number of predictions per class for each data set
    csum.train = colSums(conf.train)
    csum.valid = colSums(conf.valid)
    csum.test = colSums(conf.test)
    
    # get the proportion of observations per class for each data set
    p.train = rsum.train / n.train
    p.valid = rsum.valid / n.valid
    p.test = rsum.test / n.test
    
    # get the proportion of predcitions per class for each data set
    q.train = csum.train / n.train
    q.valid = csum.valid / n.valid
    q.test = csum.test / n.test
    
    # compute accuracy for each data set
    acc.train = sum(dia.train) / n.train
    acc.valid = sum(dia.valid) / n.valid
    acc.test = sum(dia.test) / n.test
    
    # compute expected accuracy for each data set
    exp.acc.train = sum(p.train * q.train)
    exp.acc.valid = sum(p.valid * q.valid)
    exp.acc.test = sum(p.test * q.test)
    
    # compute kappa for each data set
    kap.train = (acc.train - exp.acc.train) / (1 - exp.acc.train)
    kap.valid = (acc.valid - exp.acc.valid) / (1 - exp.acc.valid)
    kap.test = (acc.test - exp.acc.test) / (1 - exp.acc.test)
    
    # ---- compute one-vs-all metrics ----
    
    # compute a binary confusion matrix for each class, for each data set
    one.v.all.train = lapply(1:nrow(conf.train), function(i)
    {
      # get the four entries of a binary confusion matrix
      v = c(conf.train[i,i], 
            rsum.train[i] - conf.train[i,i], 
            csum.train[i] - conf.train[i,i], 
            n.train - rsum.train[i] - csum.train[i] + conf.train[i,i]);
      
      # build the confusion matrix
      return(matrix(v, nrow = 2, byrow = TRUE))
    })
    
    one.v.all.valid = lapply(1:nrow(conf.valid), function(i)
    {
      # get the four entries of a binary confusion matrix
      v = c(conf.valid[i,i], 
            rsum.valid[i] - conf.valid[i,i], 
            csum.valid[i] - conf.valid[i,i], 
            n.valid - rsum.valid[i] - csum.valid[i] + conf.valid[i,i]);
      
      # build the confusion matrix
      return(matrix(v, nrow = 2, byrow = TRUE))
    })
    
    one.v.all.test = lapply(1:nrow(conf.test), function(i)
    {
      # get the four entries of a binary confusion matrix
      v = c(conf.test[i,i], 
            rsum.test[i] - conf.test[i,i], 
            csum.test[i] - conf.test[i,i], 
            n.test - rsum.test[i] - csum.test[i] + conf.test[i,i]);
      
      # build the confusion matrix
      return(matrix(v, nrow = 2, byrow = TRUE))
    })
    
    # sum up all of the matrices for each data set
    one.v.all.train = Reduce('+', one.v.all.train)
    one.v.all.valid = Reduce('+', one.v.all.valid)
    one.v.all.test = Reduce('+', one.v.all.test)
    
    # compute the micro average accuracy for each data set
    micro.acc.train = sum(diag(one.v.all.train)) / sum(one.v.all.train)
    micro.acc.valid = sum(diag(one.v.all.valid)) / sum(one.v.all.valid)
    micro.acc.test = sum(diag(one.v.all.test)) / sum(one.v.all.test)
    
    # get the macro accuracy for each data set
    macro.acc.train = acc.train
    macro.acc.valid = acc.valid
    macro.acc.test = acc.test
    
    # ---- finalize output ----
    
    # build a final metrics table for gb.mod
    gb.mod.table = data.table(Model = "Gradient Boosting",
                                 Metric = c("Log_Loss", "Kappa", "Macro_Accuracy", "Micro_Accuracy"),
                                 Train = c(mll.train, kap.train, macro.acc.train, micro.acc.train),
                                 Valid = c(mll.valid, kap.valid, macro.acc.valid, micro.acc.valid),
                                 Test = c(mll.test, kap.test, macro.acc.test, micro.acc.test))
    
    # update the row and column names of the confusion matrices
    colnames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
    rownames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
    
    colnames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
    rownames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
    
    colnames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
    rownames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
    
    # build a list of confusion matrices
    gb.confusion = list("Train" = conf.train, 
                           "Valid" = conf.valid, 
                           "Test" = conf.test)
    
    # evaluate model performance
    gb.mod.table
    gb.confusion
    
    # free up RAM
    gc()
  }
  
  # ---------------
  # ---- NNET -----
  # ---------------
  
  {
    # make predictions on each data set
    ynew.train = as.matrix(predict(nnet.mod, newdata = train.YX.h2o)[,-1])
    ynew.valid = as.matrix(predict(nnet.mod, newdata = valid.YX.h2o)[,-1])
    ynew.test = as.matrix(predict(nnet.mod, newdata = test.YX.h2o)[,-1])
    
    # get the true values from each data set
    ytrue.train = as.data.frame(train.YX.h2o[,1])
    ytrue.valid = as.data.frame(valid.YX.h2o[,1])
    ytrue.test = as.data.frame(test.YX.h2o[,1])
    
    # convert the true values to factor data types
    target.levels = gsub("p", "", colnames(ynew.train))
    ytrue.train[,1] = factor(ytrue.train[,1], levels = target.levels)
    ytrue.valid[,1] = factor(ytrue.valid[,1], levels = target.levels)
    ytrue.test[,1] = factor(ytrue.test[,1], levels = target.levels)
    
    # ---- compute multi log loss ----
    
    # build a matrix indicating the true class values for each data set
    ytrue.train.mat = model.matrix(~., data = ytrue.train)[,-1]
    ytrue.train.mat = cbind(1 - rowSums(ytrue.train.mat), ytrue.train.mat)
    
    ytrue.valid.mat = model.matrix(~., data = ytrue.valid)[,-1]
    ytrue.valid.mat = cbind(1 - rowSums(ytrue.valid.mat), ytrue.valid.mat)
    
    ytrue.test.mat = model.matrix(~., data = ytrue.test)[,-1]
    ytrue.test.mat = cbind(1 - rowSums(ytrue.test.mat), ytrue.test.mat)
    
    # compute the multi-class logarithmic loss for each data set
    mll.train = MultiLogLoss(y_pred = ynew.train, y_true = ytrue.train.mat)
    mll.valid = MultiLogLoss(y_pred = ynew.valid, y_true = ytrue.valid.mat)
    mll.test = MultiLogLoss(y_pred = ynew.test, y_true = ytrue.test.mat)
    
    # free up RAM
    gc()
    
    # ---- compute kappa ----
    
    # get the predicted classes and actual classes for each data set
    ynew.train.code = apply(ynew.train, 1, which.max) - 1
    ytrue.train.code = factor(as.numeric(ytrue.train[,1]) - 1, levels = 0:(ncol(ytrue.train.mat) - 1))
    
    ynew.valid.code = apply(ynew.valid, 1, which.max) - 1
    ytrue.valid.code = factor(as.numeric(ytrue.valid[,1]) - 1, levels = 0:(ncol(ytrue.valid.mat) - 1))
    
    ynew.test.code = apply(ynew.test, 1, which.max) - 1
    ytrue.test.code = factor(as.numeric(ytrue.test[,1]) - 1, levels = 0:(ncol(ytrue.test.mat) - 1))
    
    # build a square confusion matrix for each data set
    conf.train = confusion(ytrue = ytrue.train.code, ypred = ynew.train.code)
    conf.valid = confusion(ytrue = ytrue.valid.code, ypred = ynew.valid.code)
    conf.test = confusion(ytrue = ytrue.test.code, ypred = ynew.test.code)
    
    # get the total number of observations for each data set
    n.train = sum(conf.train)
    n.valid = sum(conf.valid)
    n.test = sum(conf.test)
    
    # get the vector of correct predictions for each data set 
    dia.train = diag(conf.train)
    dia.valid = diag(conf.valid)
    dia.test = diag(conf.test)
    
    # get the vector of the number of observations per class for each data set
    rsum.train = rowSums(conf.train)
    rsum.valid = rowSums(conf.valid)
    rsum.test = rowSums(conf.test)
    
    # get the vector of the number of predictions per class for each data set
    csum.train = colSums(conf.train)
    csum.valid = colSums(conf.valid)
    csum.test = colSums(conf.test)
    
    # get the proportion of observations per class for each data set
    p.train = rsum.train / n.train
    p.valid = rsum.valid / n.valid
    p.test = rsum.test / n.test
    
    # get the proportion of predcitions per class for each data set
    q.train = csum.train / n.train
    q.valid = csum.valid / n.valid
    q.test = csum.test / n.test
    
    # compute accuracy for each data set
    acc.train = sum(dia.train) / n.train
    acc.valid = sum(dia.valid) / n.valid
    acc.test = sum(dia.test) / n.test
    
    # compute expected accuracy for each data set
    exp.acc.train = sum(p.train * q.train)
    exp.acc.valid = sum(p.valid * q.valid)
    exp.acc.test = sum(p.test * q.test)
    
    # compute kappa for each data set
    kap.train = (acc.train - exp.acc.train) / (1 - exp.acc.train)
    kap.valid = (acc.valid - exp.acc.valid) / (1 - exp.acc.valid)
    kap.test = (acc.test - exp.acc.test) / (1 - exp.acc.test)
    
    # ---- compute one-vs-all metrics ----
    
    # compute a binary confusion matrix for each class, for each data set
    one.v.all.train = lapply(1:nrow(conf.train), function(i)
    {
      # get the four entries of a binary confusion matrix
      v = c(conf.train[i,i], 
            rsum.train[i] - conf.train[i,i], 
            csum.train[i] - conf.train[i,i], 
            n.train - rsum.train[i] - csum.train[i] + conf.train[i,i]);
      
      # build the confusion matrix
      return(matrix(v, nrow = 2, byrow = TRUE))
    })
    
    one.v.all.valid = lapply(1:nrow(conf.valid), function(i)
    {
      # get the four entries of a binary confusion matrix
      v = c(conf.valid[i,i], 
            rsum.valid[i] - conf.valid[i,i], 
            csum.valid[i] - conf.valid[i,i], 
            n.valid - rsum.valid[i] - csum.valid[i] + conf.valid[i,i]);
      
      # build the confusion matrix
      return(matrix(v, nrow = 2, byrow = TRUE))
    })
    
    one.v.all.test = lapply(1:nrow(conf.test), function(i)
    {
      # get the four entries of a binary confusion matrix
      v = c(conf.test[i,i], 
            rsum.test[i] - conf.test[i,i], 
            csum.test[i] - conf.test[i,i], 
            n.test - rsum.test[i] - csum.test[i] + conf.test[i,i]);
      
      # build the confusion matrix
      return(matrix(v, nrow = 2, byrow = TRUE))
    })
    
    # sum up all of the matrices for each data set
    one.v.all.train = Reduce('+', one.v.all.train)
    one.v.all.valid = Reduce('+', one.v.all.valid)
    one.v.all.test = Reduce('+', one.v.all.test)
    
    # compute the micro average accuracy for each data set
    micro.acc.train = sum(diag(one.v.all.train)) / sum(one.v.all.train)
    micro.acc.valid = sum(diag(one.v.all.valid)) / sum(one.v.all.valid)
    micro.acc.test = sum(diag(one.v.all.test)) / sum(one.v.all.test)
    
    # get the macro accuracy for each data set
    macro.acc.train = acc.train
    macro.acc.valid = acc.valid
    macro.acc.test = acc.test
    
    # ---- finalize output ----
    
    # build a final metrics table for nnet.mod
    nnet.mod.table = data.table(Model = "Neural Network",
                                 Metric = c("Log_Loss", "Kappa", "Macro_Accuracy", "Micro_Accuracy"),
                                 Train = c(mll.train, kap.train, macro.acc.train, micro.acc.train),
                                 Valid = c(mll.valid, kap.valid, macro.acc.valid, micro.acc.valid),
                                 Test = c(mll.test, kap.test, macro.acc.test, micro.acc.test))
    
    # update the row and column names of the confusion matrices
    colnames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
    rownames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
    
    colnames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
    rownames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
    
    colnames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
    rownames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
    
    # build a list of confusion matrices
    nnet.confusion = list("Train" = conf.train, 
                           "Valid" = conf.valid, 
                           "Test" = conf.test)
    
    # evaluate model performance
    nnet.mod.table
    nnet.confusion
    
    # free up RAM
    gc()
  }
  
  # ---------------
  # ---- STACK ----
  # ---------------
  
  {
    # make predictions on each data set
    ynew.train = as.matrix(predict(stack.mod, newdata = train.YX.h2o)[,-1])
    ynew.valid = as.matrix(predict(stack.mod, newdata = valid.YX.h2o)[,-1])
    ynew.test = as.matrix(predict(stack.mod, newdata = test.YX.h2o)[,-1])
    
    # get the true values from each data set
    ytrue.train = as.data.frame(train.YX.h2o[,1])
    ytrue.valid = as.data.frame(valid.YX.h2o[,1])
    ytrue.test = as.data.frame(test.YX.h2o[,1])
    
    # convert the true values to factor data types
    target.levels = gsub("p", "", colnames(ynew.train))
    ytrue.train[,1] = factor(ytrue.train[,1], levels = target.levels)
    ytrue.valid[,1] = factor(ytrue.valid[,1], levels = target.levels)
    ytrue.test[,1] = factor(ytrue.test[,1], levels = target.levels)
    
    # ---- compute multi log loss ----
    
    # build a matrix indicating the true class values for each data set
    ytrue.train.mat = model.matrix(~., data = ytrue.train)[,-1]
    ytrue.train.mat = cbind(1 - rowSums(ytrue.train.mat), ytrue.train.mat)
    
    ytrue.valid.mat = model.matrix(~., data = ytrue.valid)[,-1]
    ytrue.valid.mat = cbind(1 - rowSums(ytrue.valid.mat), ytrue.valid.mat)
    
    ytrue.test.mat = model.matrix(~., data = ytrue.test)[,-1]
    ytrue.test.mat = cbind(1 - rowSums(ytrue.test.mat), ytrue.test.mat)
    
    # compute the multi-class logarithmic loss for each data set
    mll.train = MultiLogLoss(y_pred = ynew.train, y_true = ytrue.train.mat)
    mll.valid = MultiLogLoss(y_pred = ynew.valid, y_true = ytrue.valid.mat)
    mll.test = MultiLogLoss(y_pred = ynew.test, y_true = ytrue.test.mat)
    
    # free up RAM
    gc()
    
    # ---- compute kappa ----
    
    # get the predicted classes and actual classes for each data set
    ynew.train.code = apply(ynew.train, 1, which.max) - 1
    ytrue.train.code = factor(as.numeric(ytrue.train[,1]) - 1, levels = 0:(ncol(ytrue.train.mat) - 1))
    
    ynew.valid.code = apply(ynew.valid, 1, which.max) - 1
    ytrue.valid.code = factor(as.numeric(ytrue.valid[,1]) - 1, levels = 0:(ncol(ytrue.valid.mat) - 1))
    
    ynew.test.code = apply(ynew.test, 1, which.max) - 1
    ytrue.test.code = factor(as.numeric(ytrue.test[,1]) - 1, levels = 0:(ncol(ytrue.test.mat) - 1))
    
    # build a square confusion matrix for each data set
    conf.train = confusion(ytrue = ytrue.train.code, ypred = ynew.train.code)
    conf.valid = confusion(ytrue = ytrue.valid.code, ypred = ynew.valid.code)
    conf.test = confusion(ytrue = ytrue.test.code, ypred = ynew.test.code)
    
    # get the total number of observations for each data set
    n.train = sum(conf.train)
    n.valid = sum(conf.valid)
    n.test = sum(conf.test)
    
    # get the vector of correct predictions for each data set 
    dia.train = diag(conf.train)
    dia.valid = diag(conf.valid)
    dia.test = diag(conf.test)
    
    # get the vector of the number of observations per class for each data set
    rsum.train = rowSums(conf.train)
    rsum.valid = rowSums(conf.valid)
    rsum.test = rowSums(conf.test)
    
    # get the vector of the number of predictions per class for each data set
    csum.train = colSums(conf.train)
    csum.valid = colSums(conf.valid)
    csum.test = colSums(conf.test)
    
    # get the proportion of observations per class for each data set
    p.train = rsum.train / n.train
    p.valid = rsum.valid / n.valid
    p.test = rsum.test / n.test
    
    # get the proportion of predcitions per class for each data set
    q.train = csum.train / n.train
    q.valid = csum.valid / n.valid
    q.test = csum.test / n.test
    
    # compute accuracy for each data set
    acc.train = sum(dia.train) / n.train
    acc.valid = sum(dia.valid) / n.valid
    acc.test = sum(dia.test) / n.test
    
    # compute expected accuracy for each data set
    exp.acc.train = sum(p.train * q.train)
    exp.acc.valid = sum(p.valid * q.valid)
    exp.acc.test = sum(p.test * q.test)
    
    # compute kappa for each data set
    kap.train = (acc.train - exp.acc.train) / (1 - exp.acc.train)
    kap.valid = (acc.valid - exp.acc.valid) / (1 - exp.acc.valid)
    kap.test = (acc.test - exp.acc.test) / (1 - exp.acc.test)
    
    # ---- compute one-vs-all metrics ----
    
    # compute a binary confusion matrix for each class, for each data set
    one.v.all.train = lapply(1:nrow(conf.train), function(i)
    {
      # get the four entries of a binary confusion matrix
      v = c(conf.train[i,i], 
            rsum.train[i] - conf.train[i,i], 
            csum.train[i] - conf.train[i,i], 
            n.train - rsum.train[i] - csum.train[i] + conf.train[i,i]);
      
      # build the confusion matrix
      return(matrix(v, nrow = 2, byrow = TRUE))
    })
    
    one.v.all.valid = lapply(1:nrow(conf.valid), function(i)
    {
      # get the four entries of a binary confusion matrix
      v = c(conf.valid[i,i], 
            rsum.valid[i] - conf.valid[i,i], 
            csum.valid[i] - conf.valid[i,i], 
            n.valid - rsum.valid[i] - csum.valid[i] + conf.valid[i,i]);
      
      # build the confusion matrix
      return(matrix(v, nrow = 2, byrow = TRUE))
    })
    
    one.v.all.test = lapply(1:nrow(conf.test), function(i)
    {
      # get the four entries of a binary confusion matrix
      v = c(conf.test[i,i], 
            rsum.test[i] - conf.test[i,i], 
            csum.test[i] - conf.test[i,i], 
            n.test - rsum.test[i] - csum.test[i] + conf.test[i,i]);
      
      # build the confusion matrix
      return(matrix(v, nrow = 2, byrow = TRUE))
    })
    
    # sum up all of the matrices for each data set
    one.v.all.train = Reduce('+', one.v.all.train)
    one.v.all.valid = Reduce('+', one.v.all.valid)
    one.v.all.test = Reduce('+', one.v.all.test)
    
    # compute the micro average accuracy for each data set
    micro.acc.train = sum(diag(one.v.all.train)) / sum(one.v.all.train)
    micro.acc.valid = sum(diag(one.v.all.valid)) / sum(one.v.all.valid)
    micro.acc.test = sum(diag(one.v.all.test)) / sum(one.v.all.test)
    
    # get the macro accuracy for each data set
    macro.acc.train = acc.train
    macro.acc.valid = acc.valid
    macro.acc.test = acc.test
    
    # ---- finalize output ----
    
    # build a final metrics table for stack.mod
    stack.mod.table = data.table(Model = "Super Learner",
                                 Metric = c("Log_Loss", "Kappa", "Macro_Accuracy", "Micro_Accuracy"),
                                 Train = c(mll.train, kap.train, macro.acc.train, micro.acc.train),
                                 Valid = c(mll.valid, kap.valid, macro.acc.valid, micro.acc.valid),
                                 Test = c(mll.test, kap.test, macro.acc.test, micro.acc.test))
    
    # update the row and column names of the confusion matrices
    colnames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
    rownames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
    
    colnames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
    rownames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
    
    colnames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
    rownames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
    
    # build a list of confusion matrices
    stack.confusion = list("Train" = conf.train, 
                           "Valid" = conf.valid, 
                           "Test" = conf.test)
    
    # evaluate model pestackormance
    stack.mod.table
    stack.confusion
    
    # free up RAM
    gc()
  }
  
} else
{
    
  # ---------------------------------------------------------------------------------
  # ---- Cluster the Data -----------------------------------------------------------
  # ---------------------------------------------------------------------------------
  
  {
    
    # choose the number of workers for parallel processing
    workers = getDTthreads()
    
    # initialize the h2o instance
    h2o.init(nthreads = workers, max_mem_size = "9g")
    
    # remove any objects in the h2o instance
    h2o.removeAll()
    
    # remove the progress bar when model building
    h2o.no_progress()
    
    # create an h2o object of dat
    dat.h2o = as.h2o(dat)
    
    # get the column names of dat
    x = names(dat)
    
    # determine a set of clusters to try
    min.k = 4
    max.k = 36
    tot.k = 12
    set.k = unique(round(lseq(from = min.k, to = max.k, length.out = tot.k), 0))
    
    # set up hyperparameters of interest
    km.hyper.params = list(k = set.k, max_iterations = 10)
    
    # lets use a random grid search and specify a time limit and/or model limit
    minutes = 20
    km.search.criteria = list(strategy = "RandomDiscrete", 
                              max_runtime_secs = minutes * 60, 
                              # max_models = 100, 
                              seed = 42)
    
    # lets run a grid search for a good model, without drop out ratios
    h2o.rm("km.random.grid")
    km.random.grid = h2o.grid(algorithm = "kmeans",
                              grid_id = "km.random.grid",
                              x = x,
                              training_frame = dat.h2o,
                              nfolds = 5,
                              fold_assignment = "Modulo",
                              seed = 21,
                              hyper_params = km.hyper.params,
                              search_criteria = km.search.criteria)
    
    # free up RAM
    gc()
    
    # rank each model in the random grids
    km.grid = h2o.getGrid("km.random.grid", sort_by  = "tot_withinss", decreasing = FALSE)
    
    # get the summary table of the grid search
    DT.km.grid = data.table(km.grid@summary_table)
    
    # convert k, max_iterations, tot_withinss to numeric data types
    DT.km.grid[, k := as.numeric(k)]
    DT.km.grid[, max_iterations := as.numeric(max_iterations)]
    DT.km.grid[, tot_withinss := as.numeric(tot_withinss)]
    
    # plot tot_withinss v. k
    twss.plot = ggplot(data = DT.km.grid, aes(x = k, y = tot_withinss, color = tot_withinss)) + 
      geom_smooth(size = 1.5, method = 'loess', color = "black", linetype = "dashed", fill = NA) + 
      geom_point(size = 7) + 
      scale_color_continuous(low = "royalblue", high = "orangered") + 
      scale_y_continuous(label = comma) + 
      ggtitle("Cluster Analysis") + 
      labs(x = "No. Clusters", y = "Total Within Sum of Squares", color = "TWSS") + 
      theme_bw(base_size = 25) +
      theme(legend.position = "none", 
            legend.key.size = unit(.25, "in"),
            plot.title = element_text(hjust = 0.5),
            # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
    
    twss.plot
    
    # pick the top model from all grid searches
    p = 3
    km.mod = h2o.getModel(km.grid@model_ids[[p]])
    
    # get clusters
    clus = predict(km.mod, newdata = dat.h2o)
    colnames(clus) = "Cluster"
    
    # pick the number of PC to retain
    r = 5
    
    # build PC
    pca.mod = h2o.prcomp(training_frame = dat.h2o, x = x, seed = 21, 
                         k = r, pca_method = "GramSVD",
                         impute_missing = TRUE, transform = "NORMALIZE")
    
    
    # get eigenvalues of pca.mod
    eigens = pca.mod@model$model_summary[1,]^2
    row.names(eigens) = "Eigenvalue"
    
    # show importance of PC
    pca.imp = rbind(pca.mod@model$model_summary, eigens)
    pca.imp
    
    # get PC
    pca = predict(pca.mod, newdata = dat.h2o)
    
    # add clus to pca
    pca = as.data.table(h2o.cbind(pca, clus))
    
    # make Cluster into a factor
    pca[, Cluster := factor(Cluster, levels = sort(unique(pca$Cluster)))]
    
    # build a simple plot to grab a legend from
    legend.plot = ggplot(pca, aes(x = PC1, y = PC2, color = Cluster)) + 
      geom_point() + 
      theme(legend.position = "top")
    
    # plot pca v. clusters
    pca.plot = ggpairs(pca, columns = names(pca)[-ncol(pca)], 
                       mapping = aes(color = Cluster), # axisLabels = "internal",
                       legend = grab_legend(legend.plot) + theme_bw(25),
                       upper = list(continuous = wrap("points", alpha = 1/3)),
                       # upper = list(continuous = wrap(ggally_cor, size = 10)),
                       lower = list(continuous = wrap("points", alpha = 1/3)),
                       # diag = list(continuous = wrap("diagAxis", labelSize = 15)),
                       diag = list(continuous = wrap("densityDiag", fill = "black", alpha = 1/3)),
                       title = "Cluster Analysis with PCA") + 
      theme_bw(25) +
      theme(plot.title = element_text(hjust = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    # pca.plot
    
    # pick x, y, and z variables
    x = pca$PC1
    y = pca$PC2
    z = pca$PC3
    
    # make Cluster into a numeric
    pca[, Cluster := as.numeric(as.character(Cluster))]
    
    # create Cluster.factor and Cluster.number columns
    pca[, Cluster.factor := factor(Cluster, levels = sort(unique(pca$Cluster)))]
    pca[, Cluster.number := as.numeric(Cluster.factor)]
    
    # plot PC v. clusters
    # https://rpubs.com/yoshio/95844
    scatter3D(x = x, y = -y, z = -z, # xlim = c(0, 11), zlim = c(0, 90),
              main = "Cluster Analysis with PCA",
              colvar = pca$Cluster.number, 
              col = jet.col(n = length(unique(pca$Cluster.number)), alpha = 2/3),
              # col = ramp.col(col = c("cornflowerblue", "forestgreen", "orangered", "mediumorchid", "chocolate"), n = length(unique(pca$Cluster)), alpha = 1/3),
              colkey = list(at = sort(unique(pca$Cluster.number)), labels = levels(pca$Cluster.factor), dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0), 
              # colkey = list(dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0),
              xlab = "\nPC1\n22.4%", ylab = "\n\nPC2\n13.5%", zlab = "\n\nPC3\n7.3%", clab = "Cluster", 
              # surf = list(x = x.pred, y = y.pred, z = z.pred, col = alpha("black", 1/4), fit = fit.values),
              theta = 20, phi = 40,
              pch = 16, ticktype = "detailed", type = "p", bty = "b2",
              # bty = "u", col.panel = "gray33", col.axis = "white", col.grid = "white", 
              cex = 1, cex.main = 2.5, cex.axis = 2, cex.lab = 2.5)
    
    x = pca$PC1
    y = pca$PC2
    z = pca$PC4
    
    scatter3D(x = x, y = -y, z = z, # xlim = c(0, 11), zlim = c(0, 90),
              main = "Cluster Analysis with PCA",
              colvar = pca$Cluster.number, 
              col = jet.col(n = length(unique(pca$Cluster.number)), alpha = 2/3),
              # col = ramp.col(col = c("cornflowerblue", "forestgreen", "orangered", "mediumorchid", "chocolate"), n = length(unique(pca$Cluster)), alpha = 1/3),
              colkey = list(at = sort(unique(pca$Cluster.number)), labels = levels(pca$Cluster.factor), dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0), 
              # colkey = list(dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0),
              xlab = "\nPC1\n22.4%", ylab = "\n\nPC2\n13.5%", zlab = "\n\nPC4\n5.5%", clab = "Cluster", 
              # surf = list(x = x.pred, y = y.pred, z = z.pred, col = alpha("black", 1/4), fit = fit.values),
              theta = 20, phi = 40,
              pch = 16, ticktype = "detailed", type = "p", bty = "b2",
              # bty = "u", col.panel = "gray33", col.axis = "white", col.grid = "white", 
              cex = 1, cex.main = 2.5, cex.axis = 2, cex.lab = 2.5)
    
    # add SolutionID to clus
    clus.dat = as.data.table(clus)
    clus.dat[, SolutionID := 1:nrow(clus.dat)]
    
    # give sol an id variable
    sol[, id := 1:nrow(sol)]
    
    # join clus.dat onto sol
    setkey(clus.dat, SolutionID)
    setkey(sol, SolutionID)
    sol = clus.dat[sol]
    
    # order sol by id an remove it
    sol = sol[order(id)]
    sol[, id := NULL]
    
    # get cluster related metrics
    clus.dat = data.table(sol[,.(fileID, SolutionID, Cluster, Markets, MarketID, Bundle, Bundle_Name, Birth_Cohort, Bundle_Cost = (Production_Cost / as.numeric(gsub(" Markets", "", Markets))) * Produce_Bundle, Supply_Qty = Selling_Qty / Bundle_Demand, Selling_Price_Low, Selling_Price_High, Reservation_Price, Surplus_Low, Surplus_High, Revenue_Low, Revenue_High, Cohort_Drop, MARR, MarketImpacted, Country_Risk, GNIpc, Birth_Mortality)])
    
    # create Debt, Profit, and Value indicators, and update Surplus
    clus.dat[, Debt_Low := ifelse(Surplus_Low < 0, -Surplus_Low, 0)]
    clus.dat[, Debt_High := ifelse(Surplus_High < 0, -Surplus_High, 0)]
    clus.dat[, Surplus_Low := ifelse(Surplus_Low < 0, 0, Surplus_Low)]
    clus.dat[, Surplus_High := ifelse(Surplus_High < 0, 0, Surplus_High)]
    clus.dat[, Profit_Low := Revenue_Low - Bundle_Cost]
    clus.dat[, Profit_High := Revenue_High - Bundle_Cost]
    clus.dat[, Value_Low := Profit_Low + Surplus_Low - Debt_Low]
    clus.dat[, Value_High := Profit_High + Surplus_High - Debt_High]
    
    # make Markets, Cohort_Drop, and MARR into factors
    clus.dat[, Markets := factor(Markets, levels = unique(Markets))]
    clus.dat[, Cohort_Drop := factor(Cohort_Drop, levels = unique(Cohort_Drop))]
    clus.dat[, MARR := factor(MARR, levels = unique(MARR))]
    
    # convert Markets, Cohort_Drop, and MARR into binary variables
    binary = data.table(clus.dat[,.(Markets, Cohort_Drop, MARR)])
    binary = data.table(model.matrix(~ ., data = binary, 
                                     contrasts.arg = lapply(binary, contrasts, contrasts = FALSE))[,-1])
    
    # remove the extra "Markets" from some column names
    # remove the "%" from some column names
    # replace "-" with "_" in some columns
    setnames(binary, gsub("-", "_", gsub("%", "", gsub(" Markets", "", names(binary)))))
    
    # append binary to clus.dat
    clus.dat = cbind(clus.dat, binary)
    
    # get a copy of clus.dat to study antigens
    clus.dat.ant = data.table(clus.dat)
    
    # give clus.dat.ant an id column to preserve row order
    clus.dat.ant[, id := 1:nrow(clus.dat.ant)]
    
    # join antigen labels to clus.dat.ant
    setkey(clus.dat.ant, Bundle)
    setkey(ant, Bundle)
    clus.dat.ant = ant[clus.dat.ant, allow.cartesian = TRUE]
    
    # order clus.dat.ant by id and remove it
    clus.dat.ant = clus.dat.ant[order(id)]
    clus.dat.ant[, id := NULL]
    
    # get the expected value of various metrics by SolutionID
    clus.dat = clus.dat[, .(Surplus_Low = sum(Surplus_Low),
                            Surplus_High = sum(Surplus_High),
                            Revenue_Low = sum(Revenue_Low),
                            Revenue_High = sum(Revenue_High),
                            Debt_Low = sum(Debt_Low),
                            Debt_High = sum(Debt_High),
                            Profit_Low = sum(Profit_Low),
                            Profit_High = sum(Profit_High),
                            Value_Low = sum(Value_Low),
                            Value_High = sum(Value_High),
                            Bundle_Cost = sum(Bundle_Cost),
                            Country_Risk = mean(Country_Risk), 
                            GNIpc = mean(GNIpc),
                            Birth_Mortality = mean(Birth_Mortality),
                            Birth_Cohort = mean(as.numeric(Birth_Cohort)),
                            Markets2 = mean(Markets2), 
                            Markets4 = mean(Markets4), 
                            Markets8 = mean(Markets8), 
                            Markets12 = mean(Markets12), 
                            Cohort_Drop1_12 = mean(Cohort_Drop1_12),
                            Cohort_Drop13_26 = mean(Cohort_Drop13_26), 
                            Cohort_Drop27_40 = mean(Cohort_Drop27_40), 
                            MARR5 = mean(MARR5), 
                            MARR10 = mean(MARR10), 
                            MARR15 = mean(MARR15), 
                            MARR20 = mean(MARR20)), 
                        by = .(Cluster, SolutionID)]
    
    # get the expected value of various metrics by Cluster
    clus.dat = clus.dat[, .(Size = length(unique(SolutionID)),
                            Surplus_Low = mean(Surplus_Low),
                            Surplus_High = mean(Surplus_High),
                            Revenue_Low = mean(Revenue_Low),
                            Revenue_High = mean(Revenue_High),
                            Debt_Low = mean(Debt_Low),
                            Debt_High = mean(Debt_High),
                            Profit_Low = mean(Profit_Low),
                            Profit_High = mean(Profit_High),
                            Value_Low = mean(Value_Low),
                            Value_High = mean(Value_High),
                            Bundle_Cost = mean(Bundle_Cost),
                            Country_Risk = mean(Country_Risk), 
                            GNIpc = mean(GNIpc),
                            Birth_Mortality = mean(Birth_Mortality),
                            Birth_Cohort = mean(Birth_Cohort),
                            Markets2 = sum(Markets2), 
                            Markets4 = sum(Markets4), 
                            Markets8 = sum(Markets8), 
                            Markets12 = sum(Markets12), 
                            Cohort_Drop1_12 = sum(Cohort_Drop1_12),
                            Cohort_Drop13_26 = sum(Cohort_Drop13_26), 
                            Cohort_Drop27_40 = sum(Cohort_Drop27_40), 
                            MARR5 = sum(MARR5), 
                            MARR10 = sum(MARR10), 
                            MARR15 = sum(MARR15), 
                            MARR20 = sum(MARR20)), 
                        by = .(Cluster)]
    
    # compute percentages for binary indicators
    clus.dat[, Markets2_p := (Markets2 / Size) * 100]
    clus.dat[, Markets4_p := (Markets4 / Size) * 100]
    clus.dat[, Markets8_p := (Markets8 / Size) * 100]
    clus.dat[, Markets12_p := (Markets12 / Size) * 100]
    
    clus.dat[, Markets2_d := (Markets2 / sum(Markets2)) * 100]
    clus.dat[, Markets4_d := (Markets4 / sum(Markets4)) * 100]
    clus.dat[, Markets8_d := (Markets8 / sum(Markets8)) * 100]
    clus.dat[, Markets12_d := (Markets12 / sum(Markets12)) * 100]
    
    clus.dat[, Cohort_Drop1_12_p := (Cohort_Drop1_12 / Size) * 100]
    clus.dat[, Cohort_Drop13_26_p := (Cohort_Drop13_26 / Size) * 100]
    clus.dat[, Cohort_Drop27_40_p := (Cohort_Drop27_40 / Size) * 100]
    
    clus.dat[, Cohort_Drop1_12_d := (Cohort_Drop1_12 / sum(Cohort_Drop1_12)) * 100]
    clus.dat[, Cohort_Drop13_26_d := (Cohort_Drop13_26 / sum(Cohort_Drop13_26)) * 100]
    clus.dat[, Cohort_Drop27_40_d := (Cohort_Drop27_40 / sum(Cohort_Drop27_40)) * 100]
    
    clus.dat[, MARR5_p := (MARR5 / Size) * 100]
    clus.dat[, MARR10_p := (MARR10 / Size) * 100]
    clus.dat[, MARR15_p := (MARR15 / Size) * 100]
    clus.dat[, MARR20_p := (MARR20 / Size) * 100]
    
    clus.dat[, MARR5_d := (MARR5 / sum(MARR5)) * 100]
    clus.dat[, MARR10_d := (MARR10 / sum(MARR10)) * 100]
    clus.dat[, MARR15_d := (MARR15 / sum(MARR15)) * 100]
    clus.dat[, MARR20_d := (MARR20 / sum(MARR20)) * 100]
    
    # get the expected value of various metrics by Antigen
    clus.dat.ant = clus.dat.ant[, .(Birth_Cohort = mean(as.numeric(Birth_Cohort)),
                                    Supply_Qty = sum(Supply_Qty)), 
                                by = .(Cluster, SolutionID, Antigen)]
    
    # compute coverage metrics by antigen
    clus.dat.ant = clus.dat.ant[, .(Coverage_per_Sol = 100 * mean(Supply_Qty / Birth_Cohort),
                                    Coverage_Total = 100 * (sum(Supply_Qty) / sum(Birth_Cohort))), 
                                by = .(Cluster, Antigen)]
    
    # widen out antigen into its own columns with Coverage values
    clus.dat.ant = dcast(data = clus.dat.ant, formula = Cluster ~ Antigen, value.var = "Coverage_Total")
    
    # rename clus.dat.ant
    setnames(clus.dat.ant, c("Cluster", "DTP_Coverage", "HepB_Coverage", "Hib_Coverage", "IPV_Coverage", "MMR_Coverage", "V_Coverage"))
    
    # append clus.dat.ant to clus.dat
    clus.dat.ant = clus.dat.ant[order(Cluster)]
    clus.dat = clus.dat[order(Cluster)]
    clus.dat = cbind(clus.dat, clus.dat.ant[,!"Cluster"])
    
    # compute the portion of the data in each cluster
    clus.dat[, Portion := round((Size / max(sol$SolutionID)) * 100, 4)]
    
    # only keep clusters that represent at least 1% of the data
    clus.dat = clus.dat[Portion >= 1]
    
    # this is how much of the data (%) is retained after removing clusters
    sum(clus.dat$Portion)
    
    # rank clus.dat based on some performance metrics
    rank.clus = data.table(clus.dat[,.(Cluster, Portion, Value_Low, GNIpc, Country_Risk, HepB_Coverage, Debt_Low, Surplus_Low, Profit_High)])
    
    # remove missing values to study correlation
    cor.clus = data.table(na.omit(rank.clus[, !c("Cluster", "Portion"), with = FALSE]))
    
    # plot correlations in the data
    cor.plot = corrplot(cor(cor.clus), 
                        diag = FALSE, tl.offset = 0.5, tl.srt = 5, tl.cex = 1, # tl.pos = "n",
                        # p.mat = cors.sig$p, insig = "blank",
                        type = "lower", method = "square", order = "FPC", # addgrid.col = "transparent", 
                        col = c("red", "blue"), bg = "white", cl.cex = 1.25)
    
    # lets use Country_Risk, HepB_Coverage, Debt_Low, and Value_Low to rank clusters
    rank.clus = rank.clus[,.(Cluster, Country_Risk, HepB_Coverage, Debt_Low, Value_Low)]
    
    # put all variables on a 0-100 scale representing bad-good
    rank.clus[, Country_Risk := rescale(-Country_Risk, to = c(0, 100))]
    rank.clus[, HepB_Coverage := rescale(HepB_Coverage, to = c(0, 100))]
    rank.clus[, Debt_Low := rescale(-Debt_Low, to = c(0, 100))]
    rank.clus[, Value_Low := rescale(Value_Low, to = c(0, 100))]
    
    # average the rankings
    rank.clus[, Rank := (Country_Risk + HepB_Coverage + Debt_Low + Value_Low) / 4]
    
    # order by Rank
    rank.clus = rank.clus[order(-Rank)]
    
    # only keep Cluster and Rank
    rank.clus = rank.clus[,.(Cluster, Rank)]
    
    # join rank.clus onto clus.dat
    setkey(rank.clus, Cluster)
    setkey(clus.dat, Cluster)
    clus.dat = rank.clus[clus.dat]
    
    # order by Rank
    clus.dat = clus.dat[order(-Rank)]
    
    # plot some of the data
    x = clus.dat$GNIpc
    y = clus.dat$Birth_Cohort
    z = clus.dat$Rank
    
    # create Cluster.factor and Cluster.number columns
    clus.dat[, Cluster.factor := factor(Cluster, levels = sort(unique(clus.dat$Cluster)))]
    clus.dat[, Cluster.number := as.numeric(Cluster.factor)]
    
    # plot characteristics of the clusters
    scatter3D(x = x, y = y, z = z, xlim = c(0, 120000), # zlim = c(0, 90),
              main = "Characterising Clusters",
              colvar = clus.dat$Cluster.number, 
              col = ggcolor(n = length(unique(clus.dat$Cluster.number)), a = 3/4),
              # col = ramp.col(col = c("cornflowerblue", "forestgreen", "orangered", "mediumorchid", "chocolate"), n = length(unique(pca$Cluster)), alpha = 1/3),
              colkey = list(at = sort(unique(clus.dat$Cluster.number)), labels = levels(clus.dat$Cluster.factor), dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0), 
              # colkey = list(dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0),
              xlab = "\n\nGNIpc\n(USD)", ylab = "\n\nAnnual\nBirths", zlab = "\nCluster Rank", clab = "Cluster", 
              # surf = list(x = x.pred, y = y.pred, z = z.pred, col = alpha("black", 1/4), fit = fit.values),
              theta = 42, phi = 20,
              pch = 16, ticktype = "detailed", type = "h", bty = "g",
              # bty = "u", col.panel = "gray33", col.axis = "white", col.grid = "white", 
              cex = 3, cex.main = 2.5, cex.axis = 1.5, cex.lab = 2.5)
    
    # plot characteristics of the clusters
    x = clus.dat$Debt_Low
    y = clus.dat$Value_Low
    z = clus.dat$Rank
    scatter3D(x = x, y = y, z = z, # xlim = c(0, 120000), # zlim = c(0, 90),
              main = "Characterising Clusters",
              colvar = clus.dat$Cluster.number, 
              col = ggcolor(n = length(unique(clus.dat$Cluster.number)), a = 3/4),
              # col = ramp.col(col = c("cornflowerblue", "forestgreen", "orangered", "mediumorchid", "chocolate"), n = length(unique(pca$Cluster)), alpha = 1/3),
              colkey = list(at = sort(unique(clus.dat$Cluster.number)), labels = levels(clus.dat$Cluster.factor), dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0), 
              # colkey = list(dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0),
              xlab = "\n\nDebt\n(USD)", ylab = "\n\n\nMarket Value\n(USD)", zlab = "\nCluster Rank", clab = "Cluster", 
              # surf = list(x = x.pred, y = y.pred, z = z.pred, col = alpha("black", 1/4), fit = fit.values),
              theta = 42, phi = 20,
              pch = 16, ticktype = "detailed", type = "h", bty = "g",
              # bty = "u", col.panel = "gray33", col.axis = "white", col.grid = "white", 
              cex = 3, cex.main = 2.5, cex.axis = 1.5, cex.lab = 2.5)
    
    # plot characteristics of the clusters
    x = clus.dat$HepB_Coverage
    y = clus.dat$Surplus_Low
    z = clus.dat$Rank
    scatter3D(x = x, y = y, z = z, # xlim = c(0, 120000), # zlim = c(0, 90),
              main = "Characterising Clusters",
              colvar = clus.dat$Cluster.number, 
              col = ggcolor(n = length(unique(clus.dat$Cluster.number)), a = 3/4),
              # col = ramp.col(col = c("cornflowerblue", "forestgreen", "orangered", "mediumorchid", "chocolate"), n = length(unique(pca$Cluster)), alpha = 1/3),
              colkey = list(at = sort(unique(clus.dat$Cluster.number)), labels = levels(clus.dat$Cluster.factor), dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0), 
              # colkey = list(dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0),
              xlab = "\n\nHepB\nCoverage\n(%)", ylab = "\n\n\nSurplus\n(USD)", zlab = "\nCluster Rank", clab = "Cluster", 
              # surf = list(x = x.pred, y = y.pred, z = z.pred, col = alpha("black", 1/4), fit = fit.values),
              theta = 42, phi = 20,
              pch = 16, ticktype = "detailed", type = "h", bty = "g",
              # bty = "u", col.panel = "gray33", col.axis = "white", col.grid = "white", 
              cex = 3, cex.main = 2.5, cex.axis = 1.5, cex.lab = 2.5)
    
    # plot characteristics of the clusters
    x = clus.dat$Country_Risk
    y = clus.dat$Profit_High
    z = clus.dat$Rank
    scatter3D(x = x, y = y, z = z, # xlim = c(0, 120000), # zlim = c(0, 90),
              main = "Characterising Clusters",
              colvar = clus.dat$Cluster.number, 
              col = ggcolor(n = length(unique(clus.dat$Cluster.number)), a = 3/4),
              # col = ramp.col(col = c("cornflowerblue", "forestgreen", "orangered", "mediumorchid", "chocolate"), n = length(unique(pca$Cluster)), alpha = 1/3),
              colkey = list(at = sort(unique(clus.dat$Cluster.number)), labels = levels(clus.dat$Cluster.factor), dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0), 
              # colkey = list(dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0),
              xlab = "\nRisk", ylab = "\n\n\nProfit\n(USD)", zlab = "\nCluster Rank", clab = "Cluster", 
              # surf = list(x = x.pred, y = y.pred, z = z.pred, col = alpha("black", 1/4), fit = fit.values),
              theta = 42, phi = 20,
              pch = 16, ticktype = "detailed", type = "h", bty = "g",
              # bty = "u", col.panel = "gray33", col.axis = "white", col.grid = "white", 
              cex = 3, cex.main = 2.5, cex.axis = 1.5, cex.lab = 2.5)
    
    # export clus.dat
    fwrite(clus.dat, "Clustering Solutions - Birth Cohort Experiment with Risk.csv")
    
  }
  
  # ---------------------------------------------------------------------------------
  # ---- Predict the Clusters -------------------------------------------------------
  # ---------------------------------------------------------------------------------
  
  {
    
    # ---- Prepare Features ----
    
    {
      
      # set up sol for modeling 
      DT = data.table(sol[,.(Cluster, SolutionID, MARR, Markets, GNIpc, Birth_Cohort, Birth_Mortality, Country_Risk, Reservation_Price)])
      
      # make MARR and Markets into factors
      DT[, MARR := factor(MARR, levels = unique(MARR))]
      DT[, Markets := factor(Markets, levels = unique(Markets))]
      
      # convert factor variables into binary variables
      binary = data.table(DT[,.(MARR, Markets)])
      binary = data.table(model.matrix(~ ., data = binary, 
                                       contrasts.arg = lapply(binary, contrasts, contrasts = FALSE))[,-1])
      
      # add binary varaibles back to DT
      DT = cbind(DT[,!c("MARR", "Markets"), with = FALSE], 
                 binary)
      
      # update column names
      setnames(DT, gsub(" Markets", "", gsub("%", "", names(DT))))
      
      # average the data to the solution level
      DT = DT[, .(GNIpc = mean(GNIpc),
                  Birth_Cohort = mean(as.numeric(Birth_Cohort)),
                  Birth_Mortality = mean(Birth_Mortality),
                  Country_Risk = mean(Country_Risk), 
                  Reservation_Price = mean(Reservation_Price), 
                  Markets2 = mean(Markets2),
                  Markets4 = mean(Markets4), 
                  Markets8 = mean(Markets8), 
                  Markets12 = mean(Markets12), 
                  MARR5 = mean(MARR5),
                  MARR10 = mean(MARR10), 
                  MARR15 = mean(MARR15), 
                  MARR20 = mean(MARR20)), 
              by = .(Cluster, SolutionID)]
      
      # remove SolutionID
      DT[, SolutionID := NULL]
      
      # only keep clusters in clus.dat
      DT = DT[Cluster %in% clus.dat$Cluster]
      
      # convert Cluster into a factor
      DT[, Cluster := factor(Cluster, levels = sort(unique(Cluster)))]
      
      # identify predictors (x) and response (y)
      y = "Cluster"
      x = names(DT[, !"Cluster", with = FALSE])
      
      # get the predictors
      X.h2o = as.h2o(DT[, x, with = FALSE])
      
      # set up hyperparameters of interest for the autoencoder
      auto.hyper.params = list(hidden = list(c(round((2/3) * length(x), 0), round((2/3)^2 * length(x), 0), round((2/3) * length(x), 0))),
                               epochs = c(30),
                               activation = "Tanh",
                               l1 = 1e-5,
                               l2 = 0,
                               rho = 0.95,
                               epsilon = 1e-8,
                               adaptive_rate = TRUE)
      
      # lets use a random grid search and specify a time limit and/or model limit
      minutes = 30
      auto.search.criteria = list(strategy = "RandomDiscrete", 
                                  max_runtime_secs = minutes * 60, 
                                  # max_models = 100, 
                                  seed = 42)
      
      # run a random grid search for a good model
      h2o.rm("auto.random.grid")
      auto.random.grid = h2o.grid(algorithm = "deeplearning",
                                  grid_id = "auto.random.grid",
                                  autoencoder = TRUE,
                                  x = x,
                                  training_frame = X.h2o,
                                  # score_each_iteration = TRUE,
                                  seed = 3,
                                  hyper_params = auto.hyper.params,
                                  search_criteria = auto.search.criteria)
      
      # rank each model in the random grid
      auto.grid = h2o.getGrid("auto.random.grid", sort_by = "rmse", decreasing = FALSE)
      
      # check out model ranking
      auto.grid
      
      # get the best model from our grid search
      auto.mod = h2o.getModel(auto.grid@model_ids[[1]])
      
      # get features from auto.mod
      features = lapply(1:3, function(i) 
        h2o.deepfeatures(object = auto.mod, data = X.h2o, layer = i))
      
      # combine features into a single table
      features = do.call("h2o.cbind", features)
      
      # add features to DT
      DT = cbind(DT, as.data.table(features))
      
      # identify predictors (x) and response (y)
      y = "Cluster"
      x = names(DT[, !"Cluster", with = FALSE])
      
      # make an h2o version of DT
      DT.h2o = as.h2o(DT)
      
      # compute the max class weight
      max.class.weight = table(unname(unlist(DT[, y, with = FALSE])))
      max.class.weight = as.numeric(max(max(max.class.weight) / max.class.weight))
      
      # set up hyperparameters of interest
      rf.hyper.params = list(ntrees = 250,
                             # min_rows = c(1, 11, 25),
                             # max_depth = c(20, 40, 60),
                             min_rows = 11,
                             max_depth = 40,
                             stopping_metric = "mean_per_class_error")
      
      # lets use a random grid search and specify a time limit and/or model limit
      minutes = 30
      rf.search.criteria = list(strategy = "RandomDiscrete", 
                                max_runtime_secs = minutes * 60, 
                                # max_models = 100, 
                                seed = 42)
      
      # lets run a grid search for a good model, without drop out ratios
      h2o.rm("rf.random.grid")
      rf.random.grid = h2o.grid(algorithm = "randomForest",
                                grid_id = "rf.random.grid",
                                y = y,
                                x = x,
                                training_frame = DT.h2o,
                                # stopping_rounds = 3,
                                # histogram_type = "RoundRobin",
                                # nfolds = 3,
                                # fold_assignment = "Stratified",
                                seed = 3,
                                # score_each_iteration = TRUE,
                                balance_classes = TRUE,
                                max_after_balance_size = max.class.weight,
                                hyper_params = rf.hyper.params,
                                search_criteria = rf.search.criteria)
      
      # free up RAM
      gc()
      
      # rank each model in the random grids
      rf.grid = h2o.getGrid("rf.random.grid", sort_by = "mean_per_class_error", decreasing = FALSE)
      
      # pick the top model from all grid searches
      imp.rf = h2o.getModel(rf.grid@model_ids[[1]])
      
      # extract variable importance
      imp = data.table(imp.rf@model$variable_importances)
      
      # make variable into a factor
      imp[, variable := factor(variable, levels = unique(variable))]
      
      # pick a cut off value
      # this value should show where importance drops the most (look for the center of the "knee")
      cutoff = 0.25
      
      # plot a barplot of variable importance
      imp.plot1 = ggplot(imp, aes(x = variable, y = scaled_importance, fill = scaled_importance, color = scaled_importance)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_hline(yintercept = cutoff, color = "blue", linetype = "dashed", size = 1.1) + 
        ggtitle("GINI Importance\nRandom Forests") +
        labs(x = "Variable", y = "Scaled Importance") +
        scale_y_continuous(labels = percent) +
        scale_fill_gradient(low = "yellow", high = "red") +
        scale_color_gradient(low = "yellow", high = "red") +
        theme_dark(25) +
        theme(legend.position = "none", 
              plot.title = element_text(hjust = 0.5), 
              axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              panel.grid.major.x = element_blank())
      
      # plot a density plot of variable importance
      imp.plot2 = ggplot(imp, aes(x = scaled_importance)) +
        geom_density(fill = "cornflowerblue", alpha = 2/3) +
        geom_vline(xintercept = cutoff, color = "red", linetype = "dashed", size = 1.1) + 
        ggtitle("GINI Importance\nRandom Forests") +
        labs(x = "Scaled Importance", y = "Density") +
        scale_x_continuous(labels = percent) +
        coord_flip() + 
        theme_bw(25) +
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      
      # show the importance plots
      grid.arrange(imp.plot1, imp.plot2, nrow = 1)
      
      # check out imp
      imp
      
      # reverse the order of the levels of variable
      imp[, variable := factor(variable, levels = rev(levels(variable)))]
      
      # plot imp
      cutoff = 0.25
      imp.plot3 = ggplot(imp, aes(x = variable, y = scaled_importance, fill = scaled_importance, color = scaled_importance)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_hline(yintercept = cutoff, color = "blue", linetype = "dashed", size = 1.1) + 
        ggtitle("Random Forests: GINI Importance") +
        labs(x = "Variable", y = "Scaled Importance") +
        scale_y_continuous(labels = percent) +
        scale_fill_gradient(low = "yellow", high = "red") +
        scale_color_gradient(low = "yellow", high = "red") +
        coord_flip() + 
        theme_dark(20) +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5), 
              axis.title.x = element_blank(),
              axis.title.y = element_blank(), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      
      imp.plot3
      
      # find which indicators meet the cutoff value
      keep.indicators = as.character(unname(unlist(imp[scaled_importance >= cutoff, .(variable)])))
      
      # find the top m indicators
      m = 15
      keep.indicators = as.character(imp$variable)[1:m]
      
      # update DT to only contain indicators of interest
      DT = DT[, c(y, keep.indicators), with = FALSE]
      
      # determine which features to categorize
      exclude.columns = c(y)
      features = data.table(DT[, !exclude.columns, with = FALSE])
      feature.names = names(features)
      
      # go through each column of features and categorize it based on 10 quantiles
      features = lapply(1:ncol(features), function(i)
      {
        # get column i from features
        v = unname(unlist(features[, i, with = FALSE]))
        
        # categorize v
        v = cut(x = v, 
                breaks = unique(c(min(v) - 1e-6, 
                                  as.numeric(quantile(v, probs = seq(0.1, 0.9, 0.1))),
                                  max(v) + 1e-6)), 
                ordered_result = TRUE)
        
        # make v into a column again
        v = data.table(v)
        
        # give v its name back
        setnames(v, feature.names[i])
        
        return(v)
      })
      
      # combine features into a single table
      features = do.call("cbind", features)
      
      # convert features into binary values
      features = data.table(model.matrix(~ ., data = features, 
                                         contrasts.arg = lapply(features, contrasts, contrasts = FALSE))[,-1])
      
      # keep a numeric version of DT
      DT.numeric = data.table(DT)
      
      # plot some features colored by cluster
      DT.plot = ggpairs(DT.numeric, columns = names(DT.numeric)[2:6], 
                        mapping = aes(color = Cluster),
                        upper = list(continuous = wrap("density", alpha = 1)),
                        lower = list(continuous = wrap("points", alpha = 1/3, size = 1.5)),
                        diag = list(continuous = wrap("densityDiag", alpha = 1)),
                        title = "Features for Clusters") + 
        theme_bw(20) +
        theme(plot.title = element_text(hjust = 0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      
      # DT.plot
      
      # create Cluster.factor and Cluster.number columns
      DT.numeric[, Cluster.factor := factor(Cluster, levels = sort(unique(Cluster)))]
      DT.numeric[, Cluster.number := as.numeric(Cluster.factor)]
      
      # plot some features colored by cluster
      x = DT.numeric$Birth_Cohort
      y = DT.numeric$GNIpc
      z = DT.numeric$DF.L3.C2
      
      # plot features for the clusters
      scatter3D(x = -x / 1e6, y = y, z = z, # xlim = c(0, 120000), # zlim = c(0, 90),
                main = "Features for Clusters",
                colvar = DT.numeric$Cluster.number, 
                col = ggcolor(n = length(unique(DT.numeric$Cluster.number)), a = 1/3),
                colkey = list(at = sort(unique(DT.numeric$Cluster.number)), labels = levels(DT.numeric$Cluster.factor), dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0), 
                xlab = "\n\nAnnual Births\n(Millions)", ylab = "\n\nGNIpc\n(USD)", zlab = "\n\nDeep Feature\nL3C2", clab = "Cluster", 
                theta = 50, phi = 20,
                pch = 16, ticktype = "detailed", type = "p", bty = "g",
                cex = 1.5, cex.main = 2.5, cex.axis = 1.5, cex.lab = 2.5)
      
      # plot some features colored by cluster
      x = DT.numeric$Reservation_Price
      y = DT.numeric$Country_Risk
      z = DT.numeric$DF.L1.C6
      
      # plot features for the clusters
      scatter3D(x = -x, y = y, z = z, # xlim = c(0, 120000), # zlim = c(0, 90),
                main = "Features for Clusters",
                colvar = DT.numeric$Cluster.number, 
                col = ggcolor(n = length(unique(DT.numeric$Cluster.number)), a = 1/3),
                colkey = list(at = sort(unique(DT.numeric$Cluster.number)), labels = levels(DT.numeric$Cluster.factor), dist = -0.18, length = 0.75, cex.axis = 1.5, cex.clab = 1.5, adj.clab = 0), 
                xlab = "\n\n\nReservation\nPrice\n(USD)", ylab = "\nRisk", zlab = "\n\nDeep Feature\nL1C6", clab = "Cluster", 
                theta = 50, phi = 20,
                pch = 16, ticktype = "detailed", type = "p", bty = "g",
                cex = 1.5, cex.main = 2.5, cex.axis = 1.5, cex.lab = 2.5)
      
      # update DT to only have binary features
      DT = cbind(DT[, exclude.columns, with = FALSE], features)
      
      # identify predictors (x) and response (y)
      y = "Cluster"
      x = names(DT[, !"Cluster", with = FALSE])
      
      # make an h2o version of DT
      DT.h2o = as.h2o(DT)
      
      # compute the max class weight
      max.class.weight = table(unname(unlist(DT[, y, with = FALSE])))
      max.class.weight = as.numeric(max(max(max.class.weight) / max.class.weight))
      
      # set up hyperparameters of interest
      rf.hyper.params = list(ntrees = 250,
                             # min_rows = c(1, 11, 25),
                             # max_depth = c(20, 40, 60),
                             min_rows = 11,
                             max_depth = 40,
                             stopping_metric = "mean_per_class_error")
      
      # lets use a random grid search and specify a time limit and/or model limit
      minutes = 30
      rf.search.criteria = list(strategy = "RandomDiscrete", 
                                max_runtime_secs = minutes * 60, 
                                # max_models = 100, 
                                seed = 42)
      
      # lets run a grid search for a good model, without drop out ratios
      h2o.rm("rf.random.grid")
      rf.random.grid = h2o.grid(algorithm = "randomForest",
                                grid_id = "rf.random.grid",
                                y = y,
                                x = x,
                                training_frame = DT.h2o,
                                # stopping_rounds = 3,
                                # nfolds = 3,
                                # fold_assignment = "Stratified",
                                seed = 3,
                                # score_each_iteration = TRUE,
                                
                                # control how many bins are created for splitting on a feature
                                nbins_cats = 3, 
                                nbins = 3,
                                nbins_top_level = 3,
                                
                                balance_classes = TRUE,
                                max_after_balance_size = max.class.weight,
                                hyper_params = rf.hyper.params,
                                search_criteria = rf.search.criteria)
      
      # free up RAM
      gc()
      
      # rank each model in the random grids
      rf.grid = h2o.getGrid("rf.random.grid", sort_by = "mean_per_class_error", decreasing = FALSE)
      
      # pick the top model from all grid searches
      imp.rf = h2o.getModel(rf.grid@model_ids[[1]])
      
      # extract variable importance
      imp = data.table(imp.rf@model$variable_importances)
      
      # make variable into a factor
      imp[, variable := factor(variable, levels = unique(variable))]
      
      # pick a cut off value
      # this value should show where importance drops the most (look for the center of the "knee")
      cutoff = 0.23
      
      # plot a barplot of variable importance
      imp.plot1 = ggplot(imp, aes(x = variable, y = scaled_importance, fill = scaled_importance, color = scaled_importance)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_hline(yintercept = cutoff, color = "blue", linetype = "dashed", size = 1.1) + 
        ggtitle("GINI Importance\nRandom Forests") +
        labs(x = "Variable", y = "Scaled Importance") +
        scale_y_continuous(labels = percent) +
        scale_fill_gradient(low = "yellow", high = "red") +
        scale_color_gradient(low = "yellow", high = "red") +
        theme_dark(25) +
        theme(legend.position = "none", 
              plot.title = element_text(hjust = 0.5), 
              axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              panel.grid.major.x = element_blank())
      
      # plot a density plot of variable importance
      imp.plot2 = ggplot(imp, aes(x = scaled_importance)) +
        geom_density(fill = "cornflowerblue", alpha = 2/3) +
        geom_vline(xintercept = cutoff, color = "red", linetype = "dashed", size = 1.1) + 
        ggtitle("GINI Importance\nRandom Forests") +
        labs(x = "Scaled Importance", y = "Density") +
        scale_x_continuous(labels = percent) +
        coord_flip() + 
        theme_bw(25) +
        theme(plot.title = element_text(hjust = 0.5),
              axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      
      # show the importance plots
      grid.arrange(imp.plot1, imp.plot2, nrow = 1)
      
      # check out imp
      imp
      
      # reverse the order of the levels of variable
      imp[, variable := factor(variable, levels = rev(levels(variable)))]
      
      # plot imp
      cutoff = 0.17
      imp.plot3 = ggplot(imp, aes(x = variable, y = scaled_importance, fill = scaled_importance, color = scaled_importance)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_hline(yintercept = cutoff, color = "blue", linetype = "dashed", size = 1.1) + 
        ggtitle("Random Forests: GINI Importance") +
        labs(x = "Variable", y = "Scaled Importance") +
        scale_y_continuous(labels = percent) +
        scale_fill_gradient(low = "yellow", high = "red") +
        scale_color_gradient(low = "yellow", high = "red") +
        coord_flip() + 
        theme_dark(20) +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5), 
              axis.title.x = element_blank(),
              axis.title.y = element_blank(), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      
      imp.plot3
      
      # find which indicators meet the cutoff value
      keep.indicators = as.character(unname(unlist(imp[scaled_importance >= cutoff, .(variable)])))
      
      # find the top m indicators
      m = 50
      keep.indicators = as.character(imp$variable)[1:m]
      
      # update DT to only contain indicators of interest
      DT = DT[, c(y, keep.indicators), with = FALSE]
      
      # identify predictors (x) and response (y)
      y = "Cluster"
      x = names(DT[, !"Cluster", with = FALSE])
      
      # build the fold assignment
      set.seed(42)
      k.folds = 5
      folds = createFolds(y = unname(unlist(DT[, y, with = FALSE])), k = k.folds)
      
      # split up DT into train, valid, and test
      train.rows = unname(unlist(lapply(1:(k.folds - 2), function(f) folds[[f]])))
      train = data.table(DT[train.rows])
      
      valid.rows = unname(unlist(folds[[k.folds - 1]]))
      valid = data.table(DT[valid.rows])
      
      test.rows = unname(unlist(folds[[k.folds]]))
      test = data.table(DT[test.rows])
      
      # split up YX.h2o into train, valid, and test
      train.YX.h2o = as.h2o(train[, c(y, x), with = FALSE])
      valid.YX.h2o = as.h2o(valid[, c(y, x), with = FALSE])
      test.YX.h2o = as.h2o(test[, c(y, x), with = FALSE])
      
      # write out train, valid, test, and DT.numeric
      fwrite(train, "train.csv")
      fwrite(valid, "valid.csv")
      fwrite(test, "test.csv")
      fwrite(DT.numeric, "numeric.csv")
    }
    
    # ---- GLM ----
    
    {
      
      # ---- grid search for models ----
      
      # set up hyperparameters of interest
      glm.hyper.params = list(lambda = c(1, 0.5, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0),
                              alpha = c(0, 0.5, 1))
      
      # lets use a random grid search and specify a time limit and/or model limit
      minutes = 15
      glm.search.criteria = list(strategy = "RandomDiscrete", 
                                 max_runtime_secs = minutes * 60, 
                                 # max_models = 100, 
                                 seed = 42)
      
      # lets run a grid search for a good model with intercept = FALSE and standardize = FALSE
      h2o.rm("glm.random.gridA")
      glm.random.gridA = h2o.grid(algorithm = "glm",
                                  grid_id = "glm.random.gridA",
                                  y = y,
                                  x = x,
                                  training_frame = train.YX.h2o,
                                  validation_frame = valid.YX.h2o,
                                  early_stopping = TRUE,
                                  nfolds = 3,
                                  score_each_iteration = TRUE,
                                  keep_cross_validation_predictions = TRUE,
                                  fold_assignment = "Modulo",
                                  family = "multinomial",
                                  max_iterations = 100,
                                  intercept = FALSE,
                                  standardize = FALSE,
                                  seed = 21,
                                  solver = "COORDINATE_DESCENT",
                                  hyper_params = glm.hyper.params,
                                  search_criteria = glm.search.criteria)
      
      # lets run a grid search for a good model with intercept = TRUE and standardize = FALSE
      h2o.rm("glm.random.gridB")
      glm.random.gridB = h2o.grid(algorithm = "glm",
                                  grid_id = "glm.random.gridB",
                                  y = y,
                                  x = x,
                                  training_frame = train.YX.h2o,
                                  validation_frame = valid.YX.h2o,
                                  early_stopping = TRUE,
                                  nfolds = 3,
                                  score_each_iteration = TRUE,
                                  keep_cross_validation_predictions = TRUE,
                                  fold_assignment = "Modulo",
                                  family = "multinomial",
                                  max_iterations = 100,
                                  intercept = TRUE,
                                  standardize = FALSE,
                                  seed = 21,
                                  solver = "COORDINATE_DESCENT",
                                  hyper_params = glm.hyper.params,
                                  search_criteria = glm.search.criteria)
      
      # lets run a grid search for a good model with intercept = TRUE and standardize = TRUE
      h2o.rm("glm.random.gridC")
      glm.random.gridC = h2o.grid(algorithm = "glm",
                                  grid_id = "glm.random.gridC",
                                  y = y,
                                  x = x,
                                  training_frame = train.YX.h2o,
                                  validation_frame = valid.YX.h2o,
                                  early_stopping = TRUE,
                                  nfolds = 3,
                                  score_each_iteration = TRUE,
                                  keep_cross_validation_predictions = TRUE,
                                  fold_assignment = "Modulo",
                                  family = "multinomial",
                                  max_iterations = 100,
                                  intercept = TRUE,
                                  standardize = TRUE,
                                  seed = 21,
                                  solver = "COORDINATE_DESCENT",
                                  hyper_params = glm.hyper.params,
                                  search_criteria = glm.search.criteria)
      
      # lets run a grid search for a good model with intercept = FALSE and standardize = TRUE
      h2o.rm("glm.random.gridD")
      glm.random.gridD = h2o.grid(algorithm = "glm",
                                  grid_id = "glm.random.gridD",
                                  y = y,
                                  x = x,
                                  training_frame = train.YX.h2o,
                                  validation_frame = valid.YX.h2o,
                                  early_stopping = TRUE,
                                  nfolds = 3,
                                  score_each_iteration = TRUE,
                                  keep_cross_validation_predictions = TRUE,
                                  fold_assignment = "Modulo",
                                  family = "multinomial",
                                  max_iterations = 100,
                                  intercept = FALSE,
                                  standardize = TRUE,
                                  seed = 21,
                                  solver = "COORDINATE_DESCENT",
                                  hyper_params = glm.hyper.params,
                                  search_criteria = glm.search.criteria)
      
      # free up RAM
      gc()
      
      # rank each model in the random grids
      glm.gridA = h2o.getGrid("glm.random.gridA", sort_by = "mean_per_class_error", decreasing = FALSE)
      glm.gridB = h2o.getGrid("glm.random.gridB", sort_by = "mean_per_class_error", decreasing = FALSE)
      glm.gridC = h2o.getGrid("glm.random.gridC", sort_by = "mean_per_class_error", decreasing = FALSE)
      glm.gridD = h2o.getGrid("glm.random.gridD", sort_by = "mean_per_class_error", decreasing = FALSE)
      
      # combine all the grid tables into one grid table that considers the options for intercept and standardize
      glm.grid = rbind(cbind(data.table(glm.gridA@summary_table), intercept = "FALSE", standardize = "FALSE", search = 1),
                       cbind(data.table(glm.gridB@summary_table), intercept = "TRUE", standardize = "FALSE", search = 2),
                       cbind(data.table(glm.gridC@summary_table), intercept = "TRUE", standardize = "TRUE", search = 3),
                       cbind(data.table(glm.gridD@summary_table), intercept = "FALSE", standardize = "TRUE", search = 4))
      
      # combine all the grid models
      glm.grid.models = list(glm.gridA, glm.gridB, glm.gridC, glm.gridD)
      
      # order grid by mean_per_class_error
      glm.grid = glm.grid[order(as.numeric(mean_per_class_error), decreasing = FALSE)]
      
      # get the summary table of the grid search
      DT.glm.grid = data.table(glm.grid)
      DT.glm.grid
      
      # set up the data types for each column in DT.grid for plotting purposes
      DT.glm.grid = DT.glm.grid[, .(alpha = factor(as.numeric(gsub("[", "", gsub("]", "", alpha, fixed = TRUE), fixed = TRUE)), levels = c(1, 0.5, 0)),
                                    lambda = factor(as.numeric(gsub("[", "", gsub("]", "", lambda, fixed = TRUE), fixed = TRUE)), levels = c(1, 0.5, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0)),
                                    intercept = as.factor(intercept),
                                    standardize = as.factor(standardize),
                                    model = removePunctuation(gsub("[A-z]+", "", model_ids)),
                                    mean_per_class_error = as.numeric(mean_per_class_error),
                                    search = as.numeric(search))]
      
      # plot mean_per_class_error v. standardize, lambda, and alpha to see which structure is most robust
      plot.glm.grid = ggplot(DT.glm.grid, aes(x = lambda, y = mean_per_class_error, color = standardize, fill = standardize)) + 
        # geom_boxplot() + 
        geom_jitter(size = 3, alpha = 2/3) + 
        # scale_y_continuous(labels = dollar) + 
        ggtitle("Cross Validation Error") + 
        labs(x = "Strength of Regularization", y = "Log Loss", color = "Standardize", fill = "Standardize") + 
        facet_wrap(~paste("L1/L2 Distribution:", alpha), nrow = 1) +
        theme_bw(base_size = 25) +
        theme(legend.position = "top", 
              legend.key.size = unit(.25, "in"),
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
      
      plot.glm.grid
      
      # stack models?
      stack.models = FALSE
      if(stack.models)
      {
        # find the top k models from all grid searches
        k = 10
        glm.mod.list = lapply(1:k, function(i)
        {
          # get grid search (s), the model id (m), and the position of the model (p)
          s = DT.glm.grid$search[i]
          m = DT.glm.grid$model[i]
          p = which(removePunctuation(gsub("[A-z]+", "", unlist(glm.grid.models[[s]]@model_ids))) == m)
          
          # get the model of interest
          output = glm.grid.models[[s]]@model_ids[[p]]
          
          return(output)
        })
        
        # stack the top k models from all grid searches
        h2o.rm("glm.mod")
        glm.mod = h2o.stackedEnsemble(x = x,
                                      y = y,
                                      training_frame = train.YX.h2o,
                                      validation_frame = valid.YX.h2o,
                                      model_id = "glm.mod",
                                      base_models = glm.mod.list)
        
      } else
      {
        # pick the top model from all grid searches
        glm.mod = h2o.getModel(glm.grid.models[[DT.glm.grid$search[1]]]@model_ids[[1]])
      }
      
      
      
      # make predictions on each data set
      ynew.train = as.matrix(predict(glm.mod, newdata = train.YX.h2o)[,-1])
      ynew.valid = as.matrix(predict(glm.mod, newdata = valid.YX.h2o)[,-1])
      ynew.test = as.matrix(predict(glm.mod, newdata = test.YX.h2o)[,-1])
      
      # get the true values from each data set
      ytrue.train = as.data.frame(train.YX.h2o[,1])
      ytrue.valid = as.data.frame(valid.YX.h2o[,1])
      ytrue.test = as.data.frame(test.YX.h2o[,1])
      
      # ---- compute multi log loss ----
      
      # build a matrix indicating the true class values for each data set
      ytrue.train.mat = model.matrix(~., data = ytrue.train)[,-1]
      ytrue.train.mat = cbind(1 - rowSums(ytrue.train.mat), ytrue.train.mat)
      
      ytrue.valid.mat = model.matrix(~., data = ytrue.valid)[,-1]
      ytrue.valid.mat = cbind(1 - rowSums(ytrue.valid.mat), ytrue.valid.mat)
      
      ytrue.test.mat = model.matrix(~., data = ytrue.test)[,-1]
      ytrue.test.mat = cbind(1 - rowSums(ytrue.test.mat), ytrue.test.mat)
      
      # compute the multi-class logarithmic loss for each data set
      mll.train = MultiLogLoss(y_pred = ynew.train, y_true = ytrue.train.mat)
      mll.valid = MultiLogLoss(y_pred = ynew.valid, y_true = ytrue.valid.mat)
      mll.test = MultiLogLoss(y_pred = ynew.test, y_true = ytrue.test.mat)
      
      # free up RAM
      gc()
      
      # ---- compute kappa ----
      
      # get the predicted classes and actual classes for each data set
      ynew.train.code = apply(ynew.train, 1, which.max) - 1
      ytrue.train.code = factor(as.numeric(ytrue.train[,1]) - 1, levels = 0:(ncol(ytrue.train.mat) - 1))
      
      ynew.valid.code = apply(ynew.valid, 1, which.max) - 1
      ytrue.valid.code = factor(as.numeric(ytrue.valid[,1]) - 1, levels = 0:(ncol(ytrue.valid.mat) - 1))
      
      ynew.test.code = apply(ynew.test, 1, which.max) - 1
      ytrue.test.code = factor(as.numeric(ytrue.test[,1]) - 1, levels = 0:(ncol(ytrue.test.mat) - 1))
      
      # build a square confusion matrix for each data set
      conf.train = confusion(ytrue = ytrue.train.code, ypred = ynew.train.code)
      conf.valid = confusion(ytrue = ytrue.valid.code, ypred = ynew.valid.code)
      conf.test = confusion(ytrue = ytrue.test.code, ypred = ynew.test.code)
      
      # get the total number of observations for each data set
      n.train = sum(conf.train)
      n.valid = sum(conf.valid)
      n.test = sum(conf.test)
      
      # get the vector of correct predictions for each data set 
      dia.train = diag(conf.train)
      dia.valid = diag(conf.valid)
      dia.test = diag(conf.test)
      
      # get the vector of the number of observations per class for each data set
      rsum.train = rowSums(conf.train)
      rsum.valid = rowSums(conf.valid)
      rsum.test = rowSums(conf.test)
      
      # get the vector of the number of predictions per class for each data set
      csum.train = colSums(conf.train)
      csum.valid = colSums(conf.valid)
      csum.test = colSums(conf.test)
      
      # get the proportion of observations per class for each data set
      p.train = rsum.train / n.train
      p.valid = rsum.valid / n.valid
      p.test = rsum.test / n.test
      
      # get the proportion of predcitions per class for each data set
      q.train = csum.train / n.train
      q.valid = csum.valid / n.valid
      q.test = csum.test / n.test
      
      # compute accuracy for each data set
      acc.train = sum(dia.train) / n.train
      acc.valid = sum(dia.valid) / n.valid
      acc.test = sum(dia.test) / n.test
      
      # compute expected accuracy for each data set
      exp.acc.train = sum(p.train * q.train)
      exp.acc.valid = sum(p.valid * q.valid)
      exp.acc.test = sum(p.test * q.test)
      
      # compute kappa for each data set
      kap.train = (acc.train - exp.acc.train) / (1 - exp.acc.train)
      kap.valid = (acc.valid - exp.acc.valid) / (1 - exp.acc.valid)
      kap.test = (acc.test - exp.acc.test) / (1 - exp.acc.test)
      
      # ---- compute one-vs-all metrics ----
      
      # compute a binary confusion matrix for each class, for each data set
      one.v.all.train = lapply(1:nrow(conf.train), function(i)
      {
        # get the four entries of a binary confusion matrix
        v = c(conf.train[i,i], 
              rsum.train[i] - conf.train[i,i], 
              csum.train[i] - conf.train[i,i], 
              n.train - rsum.train[i] - csum.train[i] + conf.train[i,i]);
        
        # build the confusion matrix
        return(matrix(v, nrow = 2, byrow = TRUE))
      })
      
      one.v.all.valid = lapply(1:nrow(conf.valid), function(i)
      {
        # get the four entries of a binary confusion matrix
        v = c(conf.valid[i,i], 
              rsum.valid[i] - conf.valid[i,i], 
              csum.valid[i] - conf.valid[i,i], 
              n.valid - rsum.valid[i] - csum.valid[i] + conf.valid[i,i]);
        
        # build the confusion matrix
        return(matrix(v, nrow = 2, byrow = TRUE))
      })
      
      one.v.all.test = lapply(1:nrow(conf.test), function(i)
      {
        # get the four entries of a binary confusion matrix
        v = c(conf.test[i,i], 
              rsum.test[i] - conf.test[i,i], 
              csum.test[i] - conf.test[i,i], 
              n.test - rsum.test[i] - csum.test[i] + conf.test[i,i]);
        
        # build the confusion matrix
        return(matrix(v, nrow = 2, byrow = TRUE))
      })
      
      # sum up all of the matrices for each data set
      one.v.all.train = Reduce('+', one.v.all.train)
      one.v.all.valid = Reduce('+', one.v.all.valid)
      one.v.all.test = Reduce('+', one.v.all.test)
      
      # compute the micro average accuracy for each data set
      micro.acc.train = sum(diag(one.v.all.train)) / sum(one.v.all.train)
      micro.acc.valid = sum(diag(one.v.all.valid)) / sum(one.v.all.valid)
      micro.acc.test = sum(diag(one.v.all.test)) / sum(one.v.all.test)
      
      # get the macro accuracy for each data set
      macro.acc.train = acc.train
      macro.acc.valid = acc.valid
      macro.acc.test = acc.test
      
      # ---- finalize output ----
      
      # build a final metrics table for glm.mod
      glm.mod.table = data.table(Model = "Regression",
                                 Metric = c("Log_Loss", "Kappa", "Macro_Accuracy", "Micro_Accuracy"),
                                 Train = c(mll.train, kap.train, macro.acc.train, micro.acc.train),
                                 Valid = c(mll.valid, kap.valid, macro.acc.valid, micro.acc.valid),
                                 Test = c(mll.test, kap.test, macro.acc.test, micro.acc.test))
      
      # build a final grid search table
      glm.grid.search = data.table(cbind(Model = rep("Regression", nrow(DT.glm.grid)), 
                                         DT.glm.grid))
      
      # update the row and column names of the confusion matrices
      colnames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
      rownames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
      
      colnames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
      rownames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
      
      colnames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
      rownames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
      
      # build a list of confusion matrices
      glm.confusion = list("Train" = conf.train, 
                           "Valid" = conf.valid, 
                           "Test" = conf.test)
      
      # evaluate model performance
      glm.mod.table
      glm.confusion
      
      # save the model
      h2o.saveModel(object = glm.mod, path = getwd())
      
      # load the model
      # glm.mod = h2o.loadModel(path = paste0(getwd(), "/glm.mod"))
      
      # remove some objects
      rm(glm.grid.models, glm.gridA, glm.gridB, glm.gridC, glm.gridD, glm.hyper.params, 
         dia.train, ynew.train.code, ytrue.train.code, ytrue.train.mat, 
         dia.valid, ynew.valid.code, ytrue.valid.code, ytrue.valid.mat, 
         dia.test, ynew.test.code, ytrue.test.code, ytrue.test.mat, 
         p.train, q.train, acc.train, exp.acc.train,
         p.valid, q.valid, acc.valid, exp.acc.valid,
         p.test, q.test, acc.test, exp.acc.test,
         conf.train, rsum.train, csum.train, n.train, 
         conf.valid, rsum.valid, csum.valid, n.valid, 
         conf.test, rsum.test, csum.test, n.test, 
         one.v.all.train, one.v.all.valid, one.v.all.test,
         mll.train, kap.train, macro.acc.train, micro.acc.train,
         mll.valid, kap.valid, macro.acc.valid, micro.acc.valid,
         mll.test, kap.test, macro.acc.test, micro.acc.test,
         glm.random.gridA, glm.random.gridB, glm.random.gridC, glm.random.gridD, glm.search.criteria,
         plot.glm.grid, ynew.test, ynew.train, ynew.valid, ytrue.test, ytrue.train, ytrue.valid)
      
      # free up RAM
      gc()
      
      
    }
    
    # ---- RF ----
    
    {
      
      # ---- grid search for models ----
      
      # identify predictors (x) and response (y)
      y = "Cluster"
      x = names(DT[, !"Cluster", with = FALSE])
      
      # compute the max class weight
      max.class.weight = table(unname(unlist(DT[, y, with = FALSE])))
      max.class.weight = as.numeric(max(max(max.class.weight) / max.class.weight))
      
      # set up hyperparameters of interest, with drop out ratios
      rf.hyper.params = list(
        
        # -- initial tuning --
        # ntrees = c(50, 250, 500),
        min_rows = c(1, 5, 11),
        # max_depth = c(10, 20, 40),
        min_split_improvement = c(0, 1e-5),
        
        # -- re-fine tuning --
        ntrees = 50,
        # min_rows = 11,
        max_depth = 40,
        # min_split_improvement = 0,
        
        # -- scoring function -- 
        stopping_metric = "mean_per_class_error")
      
      # lets use a random grid search and specify a time limit and/or model limit
      minutes = 30
      rf.search.criteria = list(strategy = "RandomDiscrete", 
                                max_runtime_secs = minutes * 60, 
                                # max_models = 100, 
                                seed = 42)
      
      # lets run a grid search for a good model, without drop out ratios
      h2o.rm("rf.random.grid")
      rf.random.grid = h2o.grid(algorithm = "randomForest",
                                grid_id = "rf.random.grid",
                                y = y,
                                x = x,
                                training_frame = train.YX.h2o,
                                validation_frame = valid.YX.h2o,
                                # stopping_rounds = 3,
                                nfolds = 3,
                                
                                # control how many bins are created for splitting on a feature
                                nbins_cats = 3, 
                                nbins = 3,
                                nbins_top_level = 3,
                                
                                score_each_iteration = TRUE,
                                keep_cross_validation_predictions = TRUE,
                                fold_assignment = "Modulo",
                                seed = 21,
                                balance_classes = TRUE,
                                max_after_balance_size = max.class.weight,
                                hyper_params = rf.hyper.params,
                                search_criteria = rf.search.criteria)
      
      # free up RAM
      gc()
      
      # rank each model in the random grids
      rf.grid = h2o.getGrid("rf.random.grid", sort_by = "mean_per_class_error", decreasing = FALSE)
      
      # get the summary table of the grid search
      DT.rf.grid = data.table(rf.grid@summary_table)
      DT.rf.grid
      
      # set up the data types for each column in DT.grid for plotting purposes
      DT.rf.grid = DT.rf.grid[, .(ntrees = as.numeric(ntrees),
                                  min_rows = as.factor(min_rows),
                                  max_depth = as.factor(max_depth),
                                  min_split_improvement = as.factor(min_split_improvement),
                                  model = removePunctuation(gsub("[A-z]+", "", model_ids)),
                                  mean_per_class_error = as.numeric(mean_per_class_error))]
      
      # plot mean_per_class_error v. max_depth and ntrees to see which structure is most robust
      plot.rf.grid = ggplot(DT.rf.grid, aes(x = max_depth, y = mean_per_class_error, color = ntrees, fill = ntrees)) + 
        # geom_boxplot() + 
        geom_jitter(size = 3, alpha = 2/3) + 
        # scale_y_continuous(labels = dollar) + 
        ggtitle("Cross Validation Error") + 
        labs(x = "Max Depth", y = "Mean Per Class Error", color = "Trees", fill = "Trees") + 
        facet_wrap(~paste("Min Rows:", min_rows), nrow = 2) +
        theme_bw(base_size = 25) +
        theme(legend.position = "top", 
              legend.key.size = unit(.25, "in"),
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE),
               fill = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
      
      plot.rf.grid
      
      # stack models?
      stack.models = FALSE
      if(stack.models)
      {
        # find the top k models from all grid searches
        k = 5
        rf.mod.list = lapply(1:k, function(i) rf.grid@model_ids[[i]])
        
        # stack the top k models from all grid searches
        h2o.rm("rf.mod")
        rf.mod = h2o.stackedEnsemble(x = x,
                                     y = y,
                                     training_frame = train.YX.h2o,
                                     validation_frame = valid.YX.h2o,
                                     model_id = "rf.mod",
                                     base_models = rf.mod.list)
        
      } else
      {
        # pick the top model from all grid searches
        rf.mod = h2o.getModel(rf.grid@model_ids[[1]])
      }
      
      # plot variable importance
      # h2o.varimp_plot(rf.mod)
      
      
      
      # make predictions on each data set
      ynew.train = as.matrix(predict(rf.mod, newdata = train.YX.h2o)[,-1])
      ynew.valid = as.matrix(predict(rf.mod, newdata = valid.YX.h2o)[,-1])
      ynew.test = as.matrix(predict(rf.mod, newdata = test.YX.h2o)[,-1])
      
      # get the true values from each data set
      ytrue.train = as.data.frame(train.YX.h2o[,1])
      ytrue.valid = as.data.frame(valid.YX.h2o[,1])
      ytrue.test = as.data.frame(test.YX.h2o[,1])
      
      # ---- compute multi log loss ----
      
      # build a matrix indicating the true class values for each data set
      ytrue.train.mat = model.matrix(~., data = ytrue.train)[,-1]
      ytrue.train.mat = cbind(1 - rowSums(ytrue.train.mat), ytrue.train.mat)
      
      ytrue.valid.mat = model.matrix(~., data = ytrue.valid)[,-1]
      ytrue.valid.mat = cbind(1 - rowSums(ytrue.valid.mat), ytrue.valid.mat)
      
      ytrue.test.mat = model.matrix(~., data = ytrue.test)[,-1]
      ytrue.test.mat = cbind(1 - rowSums(ytrue.test.mat), ytrue.test.mat)
      
      # compute the multi-class logarithmic loss for each data set
      mll.train = MultiLogLoss(y_pred = ynew.train, y_true = ytrue.train.mat)
      mll.valid = MultiLogLoss(y_pred = ynew.valid, y_true = ytrue.valid.mat)
      mll.test = MultiLogLoss(y_pred = ynew.test, y_true = ytrue.test.mat)
      
      # free up RAM
      gc()
      
      # ---- compute kappa ----
      
      # get the predicted classes and actual classes for each data set
      ynew.train.code = apply(ynew.train, 1, which.max) - 1
      ytrue.train.code = factor(as.numeric(ytrue.train[,1]) - 1, levels = 0:(ncol(ytrue.train.mat) - 1))
      
      ynew.valid.code = apply(ynew.valid, 1, which.max) - 1
      ytrue.valid.code = factor(as.numeric(ytrue.valid[,1]) - 1, levels = 0:(ncol(ytrue.valid.mat) - 1))
      
      ynew.test.code = apply(ynew.test, 1, which.max) - 1
      ytrue.test.code = factor(as.numeric(ytrue.test[,1]) - 1, levels = 0:(ncol(ytrue.test.mat) - 1))
      
      # build a square confusion matrix for each data set
      conf.train = confusion(ytrue = ytrue.train.code, ypred = ynew.train.code)
      conf.valid = confusion(ytrue = ytrue.valid.code, ypred = ynew.valid.code)
      conf.test = confusion(ytrue = ytrue.test.code, ypred = ynew.test.code)
      
      # get the total number of observations for each data set
      n.train = sum(conf.train)
      n.valid = sum(conf.valid)
      n.test = sum(conf.test)
      
      # get the vector of correct predictions for each data set 
      dia.train = diag(conf.train)
      dia.valid = diag(conf.valid)
      dia.test = diag(conf.test)
      
      # get the vector of the number of observations per class for each data set
      rsum.train = rowSums(conf.train)
      rsum.valid = rowSums(conf.valid)
      rsum.test = rowSums(conf.test)
      
      # get the vector of the number of predictions per class for each data set
      csum.train = colSums(conf.train)
      csum.valid = colSums(conf.valid)
      csum.test = colSums(conf.test)
      
      # get the proportion of observations per class for each data set
      p.train = rsum.train / n.train
      p.valid = rsum.valid / n.valid
      p.test = rsum.test / n.test
      
      # get the proportion of predcitions per class for each data set
      q.train = csum.train / n.train
      q.valid = csum.valid / n.valid
      q.test = csum.test / n.test
      
      # compute accuracy for each data set
      acc.train = sum(dia.train) / n.train
      acc.valid = sum(dia.valid) / n.valid
      acc.test = sum(dia.test) / n.test
      
      # compute expected accuracy for each data set
      exp.acc.train = sum(p.train * q.train)
      exp.acc.valid = sum(p.valid * q.valid)
      exp.acc.test = sum(p.test * q.test)
      
      # compute kappa for each data set
      kap.train = (acc.train - exp.acc.train) / (1 - exp.acc.train)
      kap.valid = (acc.valid - exp.acc.valid) / (1 - exp.acc.valid)
      kap.test = (acc.test - exp.acc.test) / (1 - exp.acc.test)
      
      # ---- compute one-vs-all metrics ----
      
      # compute a binary confusion matrix for each class, for each data set
      one.v.all.train = lapply(1:nrow(conf.train), function(i)
      {
        # get the four entries of a binary confusion matrix
        v = c(conf.train[i,i], 
              rsum.train[i] - conf.train[i,i], 
              csum.train[i] - conf.train[i,i], 
              n.train - rsum.train[i] - csum.train[i] + conf.train[i,i]);
        
        # build the confusion matrix
        return(matrix(v, nrow = 2, byrow = TRUE))
      })
      
      one.v.all.valid = lapply(1:nrow(conf.valid), function(i)
      {
        # get the four entries of a binary confusion matrix
        v = c(conf.valid[i,i], 
              rsum.valid[i] - conf.valid[i,i], 
              csum.valid[i] - conf.valid[i,i], 
              n.valid - rsum.valid[i] - csum.valid[i] + conf.valid[i,i]);
        
        # build the confusion matrix
        return(matrix(v, nrow = 2, byrow = TRUE))
      })
      
      one.v.all.test = lapply(1:nrow(conf.test), function(i)
      {
        # get the four entries of a binary confusion matrix
        v = c(conf.test[i,i], 
              rsum.test[i] - conf.test[i,i], 
              csum.test[i] - conf.test[i,i], 
              n.test - rsum.test[i] - csum.test[i] + conf.test[i,i]);
        
        # build the confusion matrix
        return(matrix(v, nrow = 2, byrow = TRUE))
      })
      
      # sum up all of the matrices for each data set
      one.v.all.train = Reduce('+', one.v.all.train)
      one.v.all.valid = Reduce('+', one.v.all.valid)
      one.v.all.test = Reduce('+', one.v.all.test)
      
      # compute the micro average accuracy for each data set
      micro.acc.train = sum(diag(one.v.all.train)) / sum(one.v.all.train)
      micro.acc.valid = sum(diag(one.v.all.valid)) / sum(one.v.all.valid)
      micro.acc.test = sum(diag(one.v.all.test)) / sum(one.v.all.test)
      
      # get the macro accuracy for each data set
      macro.acc.train = acc.train
      macro.acc.valid = acc.valid
      macro.acc.test = acc.test
      
      # ---- finalize output ----
      
      # build a final metrics table for rf.mod
      rf.mod.table = data.table(Model = "Random Forest",
                                Metric = c("Log_Loss", "Kappa", "Macro_Accuracy", "Micro_Accuracy"),
                                Train = c(mll.train, kap.train, macro.acc.train, micro.acc.train),
                                Valid = c(mll.valid, kap.valid, macro.acc.valid, micro.acc.valid),
                                Test = c(mll.test, kap.test, macro.acc.test, micro.acc.test))
      
      # build a final grid search table
      rf.grid.search = data.table(cbind(Model = rep("Random Forest", nrow(DT.rf.grid)), 
                                        DT.rf.grid))
      
      # update the row and column names of the confusion matrices
      colnames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
      rownames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
      
      colnames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
      rownames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
      
      colnames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
      rownames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
      
      # build a list of confusion matrices
      rf.confusion = list("Train" = conf.train, 
                          "Valid" = conf.valid, 
                          "Test" = conf.test)
      
      # evaluate model performance
      rf.mod.table
      rf.confusion
      
      # save the model
      h2o.saveModel(object = rf.mod, path = getwd())
      
      # load the model
      # rf.mod = h2o.loadModel(path = paste0(getwd(), "/rf.mod"))
      
      # remove some objects
      rm(rf.hyper.params, 
         dia.train, ynew.train.code, ytrue.train.code, ytrue.train.mat, 
         dia.valid, ynew.valid.code, ytrue.valid.code, ytrue.valid.mat, 
         dia.test, ynew.test.code, ytrue.test.code, ytrue.test.mat, 
         p.train, q.train, acc.train, exp.acc.train,
         p.valid, q.valid, acc.valid, exp.acc.valid,
         p.test, q.test, acc.test, exp.acc.test,
         conf.train, rsum.train, csum.train, n.train, 
         conf.valid, rsum.valid, csum.valid, n.valid, 
         conf.test, rsum.test, csum.test, n.test,
         one.v.all.train, one.v.all.valid, one.v.all.test,
         mll.train, kap.train, macro.acc.train, micro.acc.train,
         mll.valid, kap.valid, macro.acc.valid, micro.acc.valid,
         mll.test, kap.test, macro.acc.test, micro.acc.test,
         rf.search.criteria, rf.random.grid,
         plot.rf.grid, minutes, ynew.test, ynew.train, ynew.valid, ytrue.test, ytrue.train, ytrue.valid)
      
      # free up RAM
      gc()
      
    }
    
    # ---- GB ----
    
    {
      
      # ---- grid search for models ----
      
      # compute the max class weight
      max.class.weight = table(unname(unlist(DT[, y, with = FALSE])))
      max.class.weight = as.numeric(max(max(max.class.weight) / max.class.weight))
      
      # set up hyperparameters of interest
      gb.hyper.params = list(
        
        # -- initial tuning --
        # ntrees = c(50, 250, 500),
        learn_rate = c(0.025, 0.05, 0.1),
        max_depth = c(5, 10, 20), 
        # min_rows = c(1, 5, 11),
        # sample_rate = c(0.7, 1),
        # col_sample_rate = c(0.7, 1),
        # min_split_improvement = c(0, 1e-5),
        
        # -- re-fine tuning --
        ntrees = 50,
        # learn_rate = 0.1,
        # max_depth = 5, 
        min_rows = 11,
        sample_rate = 0.7,
        col_sample_rate = 0.7,
        min_split_improvement = 0,
        
        # -- scoring function -- 
        stopping_metric = "mean_per_class_error")
      
      # lets use a random grid search and specify a time limit and/or model limit
      minutes = 30
      gb.search.criteria = list(strategy = "RandomDiscrete", 
                                max_runtime_secs = minutes * 60, 
                                # max_models = 100, 
                                seed = 42)
      
      # lets run a grid search for a good model, without drop out ratios
      h2o.rm("gb.random.grid")
      gb.random.grid = h2o.grid(algorithm = "gbm",
                                grid_id = "gb.random.grid",
                                distribution = "multinomial",
                                y = y,
                                x = x,
                                training_frame = train.YX.h2o,
                                validation_frame = valid.YX.h2o,
                                # stopping_rounds = 3,
                                nfolds = 3,
                                
                                # control how many bins are created for splitting on a feature
                                nbins_cats = 3, 
                                nbins = 3,
                                nbins_top_level = 3,
                                
                                score_each_iteration = TRUE,
                                keep_cross_validation_predictions = TRUE,
                                fold_assignment = "Modulo",
                                seed = 21,
                                balance_classes = TRUE,
                                max_after_balance_size = max.class.weight,
                                hyper_params = gb.hyper.params,
                                search_criteria = gb.search.criteria)
      
      # free up RAM
      gc()
      
      # rank each model in the random grids
      gb.grid = h2o.getGrid("gb.random.grid", sort_by = "mean_per_class_error", decreasing = FALSE)
      
      # get the summary table of the grid search
      DT.gb.grid = data.table(gb.grid@summary_table)
      DT.gb.grid
      
      # set up the data types for each column in DT.grid for plotting purposes
      DT.gb.grid = DT.gb.grid[, .(ntrees = as.numeric(ntrees),
                                  learn_rate = as.numeric(learn_rate),
                                  max_depth = as.factor(max_depth),
                                  min_rows = as.factor(min_rows),
                                  sample_rate = as.numeric(sample_rate),
                                  col_sample_rate = as.numeric(col_sample_rate),
                                  min_split_improvement = as.numeric(min_split_improvement),
                                  model = removePunctuation(gsub("[A-z]+", "", model_ids)),
                                  mean_per_class_error = as.numeric(mean_per_class_error))]
      
      # plot mean_per_class_error v. min_rows and activation to see which structure is most robust
      plot.gb.grid = ggplot(DT.gb.grid, aes(x = min_rows, y = mean_per_class_error, color = ntrees, fill = ntrees)) + 
        # geom_boxplot() + 
        geom_jitter(size = 3, alpha = 2/3) + 
        # scale_y_continuous(labels = dollar) + 
        ggtitle("Cross Validation Error") + 
        labs(x = "Minimum Child Weight", y = "Mean Per Class Error", color = "Number of Trees", fill = "Number of Trees") + 
        facet_wrap(~paste("Max Depth:", max_depth), nrow = 2) +
        theme_bw(base_size = 25) +
        theme(legend.position = "top", 
              legend.key.size = unit(.25, "in"),
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
      
      plot.gb.grid
      
      # stack models?
      stack.models = FALSE
      if(stack.models)
      {
        # find the top k models from all grid searches
        k = 10
        gb.mod.list = lapply(1:k, function(i) gb.grid@model_ids[[i]])
        
        # stack the top k models from all grid searches
        h2o.rm("gb.mod")
        gb.mod = h2o.stackedEnsemble(x = x,
                                     y = y,
                                     training_frame = train.YX.h2o,
                                     validation_frame = valid.YX.h2o,
                                     model_id = "gb.mod",
                                     base_models = gb.mod.list)
        
      } else
      {
        # pick the top model from all grid searches
        gb.mod = h2o.getModel(gb.grid@model_ids[[1]])
      }
      
      # plot variable importance
      # h2o.varimp_plot(gb.mod)
      
      
      
      # make predictions on each data set
      ynew.train = as.matrix(predict(gb.mod, newdata = train.YX.h2o)[,-1])
      ynew.valid = as.matrix(predict(gb.mod, newdata = valid.YX.h2o)[,-1])
      ynew.test = as.matrix(predict(gb.mod, newdata = test.YX.h2o)[,-1])
      
      # get the true values from each data set
      ytrue.train = as.data.frame(train.YX.h2o[,1])
      ytrue.valid = as.data.frame(valid.YX.h2o[,1])
      ytrue.test = as.data.frame(test.YX.h2o[,1])
      
      # ---- compute multi log loss ----
      
      # build a matrix indicating the true class values for each data set
      ytrue.train.mat = model.matrix(~., data = ytrue.train)[,-1]
      ytrue.train.mat = cbind(1 - rowSums(ytrue.train.mat), ytrue.train.mat)
      
      ytrue.valid.mat = model.matrix(~., data = ytrue.valid)[,-1]
      ytrue.valid.mat = cbind(1 - rowSums(ytrue.valid.mat), ytrue.valid.mat)
      
      ytrue.test.mat = model.matrix(~., data = ytrue.test)[,-1]
      ytrue.test.mat = cbind(1 - rowSums(ytrue.test.mat), ytrue.test.mat)
      
      # compute the multi-class logarithmic loss for each data set
      mll.train = MultiLogLoss(y_pred = ynew.train, y_true = ytrue.train.mat)
      mll.valid = MultiLogLoss(y_pred = ynew.valid, y_true = ytrue.valid.mat)
      mll.test = MultiLogLoss(y_pred = ynew.test, y_true = ytrue.test.mat)
      
      # free up RAM
      gc()
      
      # ---- compute kappa ----
      
      # get the predicted classes and actual classes for each data set
      ynew.train.code = apply(ynew.train, 1, which.max) - 1
      ytrue.train.code = factor(as.numeric(ytrue.train[,1]) - 1, levels = 0:(ncol(ytrue.train.mat) - 1))
      
      ynew.valid.code = apply(ynew.valid, 1, which.max) - 1
      ytrue.valid.code = factor(as.numeric(ytrue.valid[,1]) - 1, levels = 0:(ncol(ytrue.valid.mat) - 1))
      
      ynew.test.code = apply(ynew.test, 1, which.max) - 1
      ytrue.test.code = factor(as.numeric(ytrue.test[,1]) - 1, levels = 0:(ncol(ytrue.test.mat) - 1))
      
      # build a square confusion matrix for each data set
      conf.train = confusion(ytrue = ytrue.train.code, ypred = ynew.train.code)
      conf.valid = confusion(ytrue = ytrue.valid.code, ypred = ynew.valid.code)
      conf.test = confusion(ytrue = ytrue.test.code, ypred = ynew.test.code)
      
      # get the total number of observations for each data set
      n.train = sum(conf.train)
      n.valid = sum(conf.valid)
      n.test = sum(conf.test)
      
      # get the vector of correct predictions for each data set 
      dia.train = diag(conf.train)
      dia.valid = diag(conf.valid)
      dia.test = diag(conf.test)
      
      # get the vector of the number of observations per class for each data set
      rsum.train = rowSums(conf.train)
      rsum.valid = rowSums(conf.valid)
      rsum.test = rowSums(conf.test)
      
      # get the vector of the number of predictions per class for each data set
      csum.train = colSums(conf.train)
      csum.valid = colSums(conf.valid)
      csum.test = colSums(conf.test)
      
      # get the proportion of observations per class for each data set
      p.train = rsum.train / n.train
      p.valid = rsum.valid / n.valid
      p.test = rsum.test / n.test
      
      # get the proportion of predcitions per class for each data set
      q.train = csum.train / n.train
      q.valid = csum.valid / n.valid
      q.test = csum.test / n.test
      
      # compute accuracy for each data set
      acc.train = sum(dia.train) / n.train
      acc.valid = sum(dia.valid) / n.valid
      acc.test = sum(dia.test) / n.test
      
      # compute expected accuracy for each data set
      exp.acc.train = sum(p.train * q.train)
      exp.acc.valid = sum(p.valid * q.valid)
      exp.acc.test = sum(p.test * q.test)
      
      # compute kappa for each data set
      kap.train = (acc.train - exp.acc.train) / (1 - exp.acc.train)
      kap.valid = (acc.valid - exp.acc.valid) / (1 - exp.acc.valid)
      kap.test = (acc.test - exp.acc.test) / (1 - exp.acc.test)
      
      # ---- compute one-vs-all metrics ----
      
      # compute a binary confusion matrix for each class, for each data set
      one.v.all.train = lapply(1:nrow(conf.train), function(i)
      {
        # get the four entries of a binary confusion matrix
        v = c(conf.train[i,i], 
              rsum.train[i] - conf.train[i,i], 
              csum.train[i] - conf.train[i,i], 
              n.train - rsum.train[i] - csum.train[i] + conf.train[i,i]);
        
        # build the confusion matrix
        return(matrix(v, nrow = 2, byrow = TRUE))
      })
      
      one.v.all.valid = lapply(1:nrow(conf.valid), function(i)
      {
        # get the four entries of a binary confusion matrix
        v = c(conf.valid[i,i], 
              rsum.valid[i] - conf.valid[i,i], 
              csum.valid[i] - conf.valid[i,i], 
              n.valid - rsum.valid[i] - csum.valid[i] + conf.valid[i,i]);
        
        # build the confusion matrix
        return(matrix(v, nrow = 2, byrow = TRUE))
      })
      
      one.v.all.test = lapply(1:nrow(conf.test), function(i)
      {
        # get the four entries of a binary confusion matrix
        v = c(conf.test[i,i], 
              rsum.test[i] - conf.test[i,i], 
              csum.test[i] - conf.test[i,i], 
              n.test - rsum.test[i] - csum.test[i] + conf.test[i,i]);
        
        # build the confusion matrix
        return(matrix(v, nrow = 2, byrow = TRUE))
      })
      
      # sum up all of the matrices for each data set
      one.v.all.train = Reduce('+', one.v.all.train)
      one.v.all.valid = Reduce('+', one.v.all.valid)
      one.v.all.test = Reduce('+', one.v.all.test)
      
      # compute the micro average accuracy for each data set
      micro.acc.train = sum(diag(one.v.all.train)) / sum(one.v.all.train)
      micro.acc.valid = sum(diag(one.v.all.valid)) / sum(one.v.all.valid)
      micro.acc.test = sum(diag(one.v.all.test)) / sum(one.v.all.test)
      
      # get the macro accuracy for each data set
      macro.acc.train = acc.train
      macro.acc.valid = acc.valid
      macro.acc.test = acc.test
      
      # ---- finalize output ----
      
      # build a final metrics table for gb.mod
      gb.mod.table = data.table(Model = "Gradient Boosting",
                                Metric = c("Log_Loss", "Kappa", "Macro_Accuracy", "Micro_Accuracy"),
                                Train = c(mll.train, kap.train, macro.acc.train, micro.acc.train),
                                Valid = c(mll.valid, kap.valid, macro.acc.valid, micro.acc.valid),
                                Test = c(mll.test, kap.test, macro.acc.test, micro.acc.test))
      
      # build a final grid search table
      gb.grid.search = data.table(cbind(Model = rep("Gradient Boosting", nrow(DT.gb.grid)), 
                                        DT.gb.grid))
      
      # update the row and column names of the confusion matrices
      colnames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
      rownames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
      
      colnames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
      rownames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
      
      colnames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
      rownames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
      
      # build a list of confusion matrices
      gb.confusion = list("Train" = conf.train, 
                          "Valid" = conf.valid, 
                          "Test" = conf.test)
      
      # evaluate model pegbormance
      gb.mod.table
      gb.confusion
      
      # save the model
      h2o.saveModel(object = gb.mod, path = getwd())
      
      # load the model
      # gb.mod = h2o.loadModel(path = paste0(getwd(), "/gb.mod"))
      
      # remove some objects
      rm(gb.hyper.params, 
         dia.train, ynew.train.code, ytrue.train.code, ytrue.train.mat, 
         dia.valid, ynew.valid.code, ytrue.valid.code, ytrue.valid.mat, 
         dia.test, ynew.test.code, ytrue.test.code, ytrue.test.mat, 
         p.train, q.train, acc.train, exp.acc.train,
         p.valid, q.valid, acc.valid, exp.acc.valid,
         p.test, q.test, acc.test, exp.acc.test,
         conf.train, rsum.train, csum.train, n.train, 
         conf.valid, rsum.valid, csum.valid, n.valid, 
         conf.test, rsum.test, csum.test, n.test,
         one.v.all.train, one.v.all.valid, one.v.all.test,
         mll.train, kap.train, macro.acc.train, micro.acc.train,
         mll.valid, kap.valid, macro.acc.valid, micro.acc.valid,
         mll.test, kap.test, macro.acc.test, micro.acc.test,
         gb.search.criteria, gb.random.grid,
         plot.gb.grid, minutes, ynew.test, ynew.train, ynew.valid, ytrue.test, ytrue.train, ytrue.valid)
      
      # free up RAM
      gc()
      
      
      
    }
    
    # ---- NNET ----
    
    {
      
      # ---- grid search for models ----
      
      # compute the max class weight
      max.class.weight = table(unname(unlist(DT[, y, with = FALSE])))
      max.class.weight = as.numeric(max(max(max.class.weight) / max.class.weight))
      
      # compute small, medium, and large hidden layer sizes relative to the input layer
      small = round(length(x) * (1 + 0.5), 0)
      medium = round(length(x) * (1 + 0.5)^2, 0)
      large = round(length(x) * (1 + 0.5)^3, 0)
      
      # tuning order:
      # 1. hidden & rho
      # 2. epsilon & activation (only use Rectifier if you have alot of data)
      # 3. L1/L2 & ratios
      # 4. epochs
      
      # set up hyperparameters of interest
      nnet.hyper.params = list(adaptive_rate = TRUE,
                               
                               # -- initial tuning --
                               # epochs = c(10, 100, 1000),
                               # hidden = list(small, medium, large, c(small, small), c(medium, medium), c(large, large), 
                               #               c(small, small, small), c(medium, medium, medium), c(large, large, large)),
                               # activation = c("RectifierWithDropout", "TanhWithDropout"),
                               # input_dropout_ratio = c(0, 0.15),
                               # l1 = c(0, 1e-5),
                               # l2 = c(0, 1e-5),
                               # rho = c(0.9, 0.95, 0.99, 0.999),
                               # epsilon = c(1e-10, 1e-8),
                               
                               # -- re-fine tuning --
                               epochs = 100,
                               hidden = list(c(small, small, small), c(medium, medium, medium), c(large, large, large)),
                               activation = "RectifierWithDropout",
                               input_dropout_ratio = 0,
                               l1 = 0,
                               l2 = 0,
                               rho = 0.925,
                               epsilon = 1e-8,
                               
                               # -- scoring function -- 
                               stopping_metric = "mean_per_class_error")
      
      # lets use a random grid search and specify a time limit and/or model limit
      minutes = 30
      nnet.search.criteria = list(strategy = "RandomDiscrete", 
                                  max_runtime_secs = minutes * 60, 
                                  # max_models = 100, 
                                  seed = 42)
      
      # lets run a grid search for a good model
      h2o.rm("nnet.random.grid")
      nnet.random.grid = h2o.grid(algorithm = "deeplearning",
                                  grid_id = "nnet.random.grid",
                                  y = y,
                                  x = x,
                                  training_frame = train.YX.h2o,
                                  validation_frame = valid.YX.h2o,
                                  nfolds = 3,
                                  score_each_iteration = TRUE,
                                  keep_cross_validation_predictions = TRUE,
                                  score_validation_sampling = "Stratified",
                                  fold_assignment = "Modulo",
                                  seed = 21,
                                  variable_importances = FALSE,
                                  balance_classes = TRUE,
                                  max_after_balance_size = max.class.weight,
                                  hyper_params = nnet.hyper.params,
                                  search_criteria = nnet.search.criteria)
      
      # free up RAM
      gc()
      
      # rank each model in the random grids
      nnet.grid = h2o.getGrid("nnet.random.grid", sort_by = "mean_per_class_error", decreasing = FALSE)
      
      # get the summary table of the grid search
      DT.nnet.grid = data.table(nnet.grid@summary_table)
      DT.nnet.grid
      
      # set up the data types for each column in DT.grid for plotting purposes
      DT.nnet.grid = DT.nnet.grid[, .(activation = as.factor(activation),
                                      epochs = as.numeric(epochs),
                                      hidden = as.factor(hidden),
                                      input_dropout_ratio = as.factor(input_dropout_ratio),
                                      rho = as.factor(rho),
                                      epsilon = as.factor(epsilon),
                                      l1 = as.factor(l1),
                                      l2 = as.factor(l2),
                                      model = removePunctuation(gsub("[A-z]+", "", model_ids)),
                                      mean_per_class_error = as.numeric(mean_per_class_error))]
      
      # plot mean_per_class_error v. hidden and activation to see which structure is most robust
      plot.nnet.grid = ggplot(DT.nnet.grid, aes(x = hidden, y = mean_per_class_error, color = epsilon, fill = epsilon)) + 
        # geom_boxplot() + 
        geom_jitter(size = 3, alpha = 2/3) + 
        scale_y_continuous(labels = percent) + 
        ggtitle("Cross Validation Error") + 
        labs(x = "Hidden Layer Structure", y = "Mean Per Class Error", color = "Epsilon", fill = "Epsilon") + 
        facet_wrap(~paste("Learning Rate:", rho), nrow = 1) +
        theme_bw(base_size = 25) +
        theme(legend.position = "top", 
              legend.key.size = unit(.25, "in"),
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        guides(color = guide_legend(override.aes = list(size = 10, linetype = 1, alpha = 1), nrow = 1, byrow = TRUE))
      
      plot.nnet.grid
      
      # stack models?
      stack.models = FALSE
      if(stack.models)
      {
        # find the top k models from all grid searches
        k = 5
        nnet.mod.list = lapply(1:k, function(i) nnet.grid@model_ids[[i]])
        
        # stack the top k models from all grid searches
        h2o.rm("nnet.mod")
        nnet.mod = h2o.stackedEnsemble(x = x,
                                       y = y,
                                       training_frame = train.YX.h2o,
                                       validation_frame = valid.YX.h2o,
                                       model_id = "nnet.mod",
                                       base_models = nnet.mod.list)
        
      } else
      {
        # pick the top model from all grid searches
        nnet.mod = h2o.getModel(nnet.grid@model_ids[[1]])
      }
      
      # plot variable importance
      # h2o.varimp_plot(nnet.mod)
      
      
      
      # make predictions on each data set
      ynew.train = as.matrix(predict(nnet.mod, newdata = train.YX.h2o)[,-1])
      ynew.valid = as.matrix(predict(nnet.mod, newdata = valid.YX.h2o)[,-1])
      ynew.test = as.matrix(predict(nnet.mod, newdata = test.YX.h2o)[,-1])
      
      # get the true values from each data set
      ytrue.train = as.data.frame(train.YX.h2o[,1])
      ytrue.valid = as.data.frame(valid.YX.h2o[,1])
      ytrue.test = as.data.frame(test.YX.h2o[,1])
      
      # ---- compute multi log loss ----
      
      # build a matrix indicating the true class values for each data set
      ytrue.train.mat = model.matrix(~., data = ytrue.train)[,-1]
      ytrue.train.mat = cbind(1 - rowSums(ytrue.train.mat), ytrue.train.mat)
      
      ytrue.valid.mat = model.matrix(~., data = ytrue.valid)[,-1]
      ytrue.valid.mat = cbind(1 - rowSums(ytrue.valid.mat), ytrue.valid.mat)
      
      ytrue.test.mat = model.matrix(~., data = ytrue.test)[,-1]
      ytrue.test.mat = cbind(1 - rowSums(ytrue.test.mat), ytrue.test.mat)
      
      # compute the multi-class logarithmic loss for each data set
      mll.train = MultiLogLoss(y_pred = ynew.train, y_true = ytrue.train.mat)
      mll.valid = MultiLogLoss(y_pred = ynew.valid, y_true = ytrue.valid.mat)
      mll.test = MultiLogLoss(y_pred = ynew.test, y_true = ytrue.test.mat)
      
      # free up RAM
      gc()
      
      # ---- compute kappa ----
      
      # get the predicted classes and actual classes for each data set
      ynew.train.code = apply(ynew.train, 1, which.max) - 1
      ytrue.train.code = factor(as.numeric(ytrue.train[,1]) - 1, levels = 0:(ncol(ytrue.train.mat) - 1))
      
      ynew.valid.code = apply(ynew.valid, 1, which.max) - 1
      ytrue.valid.code = factor(as.numeric(ytrue.valid[,1]) - 1, levels = 0:(ncol(ytrue.valid.mat) - 1))
      
      ynew.test.code = apply(ynew.test, 1, which.max) - 1
      ytrue.test.code = factor(as.numeric(ytrue.test[,1]) - 1, levels = 0:(ncol(ytrue.test.mat) - 1))
      
      # build a square confusion matrix for each data set
      conf.train = confusion(ytrue = ytrue.train.code, ypred = ynew.train.code)
      conf.valid = confusion(ytrue = ytrue.valid.code, ypred = ynew.valid.code)
      conf.test = confusion(ytrue = ytrue.test.code, ypred = ynew.test.code)
      
      # get the total number of observations for each data set
      n.train = sum(conf.train)
      n.valid = sum(conf.valid)
      n.test = sum(conf.test)
      
      # get the vector of correct predictions for each data set 
      dia.train = diag(conf.train)
      dia.valid = diag(conf.valid)
      dia.test = diag(conf.test)
      
      # get the vector of the number of observations per class for each data set
      rsum.train = rowSums(conf.train)
      rsum.valid = rowSums(conf.valid)
      rsum.test = rowSums(conf.test)
      
      # get the vector of the number of predictions per class for each data set
      csum.train = colSums(conf.train)
      csum.valid = colSums(conf.valid)
      csum.test = colSums(conf.test)
      
      # get the proportion of observations per class for each data set
      p.train = rsum.train / n.train
      p.valid = rsum.valid / n.valid
      p.test = rsum.test / n.test
      
      # get the proportion of predcitions per class for each data set
      q.train = csum.train / n.train
      q.valid = csum.valid / n.valid
      q.test = csum.test / n.test
      
      # compute accuracy for each data set
      acc.train = sum(dia.train) / n.train
      acc.valid = sum(dia.valid) / n.valid
      acc.test = sum(dia.test) / n.test
      
      # compute expected accuracy for each data set
      exp.acc.train = sum(p.train * q.train)
      exp.acc.valid = sum(p.valid * q.valid)
      exp.acc.test = sum(p.test * q.test)
      
      # compute kappa for each data set
      kap.train = (acc.train - exp.acc.train) / (1 - exp.acc.train)
      kap.valid = (acc.valid - exp.acc.valid) / (1 - exp.acc.valid)
      kap.test = (acc.test - exp.acc.test) / (1 - exp.acc.test)
      
      # ---- compute one-vs-all metrics ----
      
      # compute a binary confusion matrix for each class, for each data set
      one.v.all.train = lapply(1:nrow(conf.train), function(i)
      {
        # get the four entries of a binary confusion matrix
        v = c(conf.train[i,i], 
              rsum.train[i] - conf.train[i,i], 
              csum.train[i] - conf.train[i,i], 
              n.train - rsum.train[i] - csum.train[i] + conf.train[i,i]);
        
        # build the confusion matrix
        return(matrix(v, nrow = 2, byrow = TRUE))
      })
      
      one.v.all.valid = lapply(1:nrow(conf.valid), function(i)
      {
        # get the four entries of a binary confusion matrix
        v = c(conf.valid[i,i], 
              rsum.valid[i] - conf.valid[i,i], 
              csum.valid[i] - conf.valid[i,i], 
              n.valid - rsum.valid[i] - csum.valid[i] + conf.valid[i,i]);
        
        # build the confusion matrix
        return(matrix(v, nrow = 2, byrow = TRUE))
      })
      
      one.v.all.test = lapply(1:nrow(conf.test), function(i)
      {
        # get the four entries of a binary confusion matrix
        v = c(conf.test[i,i], 
              rsum.test[i] - conf.test[i,i], 
              csum.test[i] - conf.test[i,i], 
              n.test - rsum.test[i] - csum.test[i] + conf.test[i,i]);
        
        # build the confusion matrix
        return(matrix(v, nrow = 2, byrow = TRUE))
      })
      
      # sum up all of the matrices for each data set
      one.v.all.train = Reduce('+', one.v.all.train)
      one.v.all.valid = Reduce('+', one.v.all.valid)
      one.v.all.test = Reduce('+', one.v.all.test)
      
      # compute the micro average accuracy for each data set
      micro.acc.train = sum(diag(one.v.all.train)) / sum(one.v.all.train)
      micro.acc.valid = sum(diag(one.v.all.valid)) / sum(one.v.all.valid)
      micro.acc.test = sum(diag(one.v.all.test)) / sum(one.v.all.test)
      
      # get the macro accuracy for each data set
      macro.acc.train = acc.train
      macro.acc.valid = acc.valid
      macro.acc.test = acc.test
      
      # ---- finalize output ----
      
      # build a final metrics table for nnet.mod
      nnet.mod.table = data.table(Model = "Neural Network",
                                  Metric = c("Log_Loss", "Kappa", "Macro_Accuracy", "Micro_Accuracy"),
                                  Train = c(mll.train, kap.train, macro.acc.train, micro.acc.train),
                                  Valid = c(mll.valid, kap.valid, macro.acc.valid, micro.acc.valid),
                                  Test = c(mll.test, kap.test, macro.acc.test, micro.acc.test))
      
      # build a final grid search table
      nnet.grid.search = data.table(cbind(Model = rep("Neural Network", nrow(DT.nnet.grid)), 
                                          DT.nnet.grid))
      
      # update the row and column names of the confusion matrices
      colnames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
      rownames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
      
      colnames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
      rownames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
      
      colnames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
      rownames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
      
      # build a list of confusion matrices
      nnet.confusion = list("Train" = conf.train, 
                            "Valid" = conf.valid, 
                            "Test" = conf.test)
      
      # evaluate model performance
      nnet.mod.table
      nnet.confusion
      
      # save the model
      h2o.saveModel(object = nnet.mod, path = getwd())
      
      # load the model
      # nnet.mod = h2o.loadModel(path = paste0(getwd(), "/nnet.mod"))
      
      # remove some objects
      rm(nnet.hyper.params, 
         dia.train, ynew.train.code, ytrue.train.code, ytrue.train.mat, 
         dia.valid, ynew.valid.code, ytrue.valid.code, ytrue.valid.mat, 
         dia.test, ynew.test.code, ytrue.test.code, ytrue.test.mat, 
         p.train, q.train, acc.train, exp.acc.train,
         p.valid, q.valid, acc.valid, exp.acc.valid,
         p.test, q.test, acc.test, exp.acc.test, 
         conf.train, rsum.train, csum.train, n.train, 
         conf.valid, rsum.valid, csum.valid, n.valid, 
         conf.test, rsum.test, csum.test, n.test,
         one.v.all.train, one.v.all.valid, one.v.all.test,
         mll.train, kap.train, macro.acc.train, micro.acc.train,
         mll.valid, kap.valid, macro.acc.valid, micro.acc.valid,
         mll.test, kap.test, macro.acc.test, micro.acc.test,
         nnet.search.criteria, nnet.random.grid,
         plot.nnet.grid, minutes, ynew.test, ynew.train, ynew.valid, ytrue.test, ytrue.train, ytrue.valid)
      
      # free up RAM
      gc()
      
    }
    
    # ---- STACK ----
    
    {
      
      # ---- stack models ----
      
      # pick models
      my.models = list(glm.mod, rf.mod, gb.mod, nnet.mod)
      
      # stack models
      h2o.rm("stack.mod")
      stack.mod = h2o.stackedEnsemble(x = x,
                                      y = y,
                                      training_frame = train.YX.h2o,
                                      validation_frame = valid.YX.h2o,
                                      model_id = "stack.mod",
                                      base_models = my.models)
      
      # make predictions on each data set
      ynew.train = as.matrix(predict(stack.mod, newdata = train.YX.h2o)[,-1])
      ynew.valid = as.matrix(predict(stack.mod, newdata = valid.YX.h2o)[,-1])
      ynew.test = as.matrix(predict(stack.mod, newdata = test.YX.h2o)[,-1])
      
      # get the true values from each data set
      ytrue.train = as.data.frame(train.YX.h2o[,1])
      ytrue.valid = as.data.frame(valid.YX.h2o[,1])
      ytrue.test = as.data.frame(test.YX.h2o[,1])
      
      # ---- compute multi log loss ----
      
      # build a matrix indicating the true class values for each data set
      ytrue.train.mat = model.matrix(~., data = ytrue.train)[,-1]
      ytrue.train.mat = cbind(1 - rowSums(ytrue.train.mat), ytrue.train.mat)
      
      ytrue.valid.mat = model.matrix(~., data = ytrue.valid)[,-1]
      ytrue.valid.mat = cbind(1 - rowSums(ytrue.valid.mat), ytrue.valid.mat)
      
      ytrue.test.mat = model.matrix(~., data = ytrue.test)[,-1]
      ytrue.test.mat = cbind(1 - rowSums(ytrue.test.mat), ytrue.test.mat)
      
      # compute the multi-class logarithmic loss for each data set
      mll.train = MultiLogLoss(y_pred = ynew.train, y_true = ytrue.train.mat)
      mll.valid = MultiLogLoss(y_pred = ynew.valid, y_true = ytrue.valid.mat)
      mll.test = MultiLogLoss(y_pred = ynew.test, y_true = ytrue.test.mat)
      
      # free up RAM
      gc()
      
      # ---- compute kappa ----
      
      # get the predicted classes and actual classes for each data set
      ynew.train.code = apply(ynew.train, 1, which.max) - 1
      ytrue.train.code = factor(as.numeric(ytrue.train[,1]) - 1, levels = 0:(ncol(ytrue.train.mat) - 1))
      
      ynew.valid.code = apply(ynew.valid, 1, which.max) - 1
      ytrue.valid.code = factor(as.numeric(ytrue.valid[,1]) - 1, levels = 0:(ncol(ytrue.valid.mat) - 1))
      
      ynew.test.code = apply(ynew.test, 1, which.max) - 1
      ytrue.test.code = factor(as.numeric(ytrue.test[,1]) - 1, levels = 0:(ncol(ytrue.test.mat) - 1))
      
      # build a square confusion matrix for each data set
      conf.train = confusion(ytrue = ytrue.train.code, ypred = ynew.train.code)
      conf.valid = confusion(ytrue = ytrue.valid.code, ypred = ynew.valid.code)
      conf.test = confusion(ytrue = ytrue.test.code, ypred = ynew.test.code)
      
      # get the total number of observations for each data set
      n.train = sum(conf.train)
      n.valid = sum(conf.valid)
      n.test = sum(conf.test)
      
      # get the vector of correct predictions for each data set 
      dia.train = diag(conf.train)
      dia.valid = diag(conf.valid)
      dia.test = diag(conf.test)
      
      # get the vector of the number of observations per class for each data set
      rsum.train = rowSums(conf.train)
      rsum.valid = rowSums(conf.valid)
      rsum.test = rowSums(conf.test)
      
      # get the vector of the number of predictions per class for each data set
      csum.train = colSums(conf.train)
      csum.valid = colSums(conf.valid)
      csum.test = colSums(conf.test)
      
      # get the proportion of observations per class for each data set
      p.train = rsum.train / n.train
      p.valid = rsum.valid / n.valid
      p.test = rsum.test / n.test
      
      # get the proportion of predcitions per class for each data set
      q.train = csum.train / n.train
      q.valid = csum.valid / n.valid
      q.test = csum.test / n.test
      
      # compute accuracy for each data set
      acc.train = sum(dia.train) / n.train
      acc.valid = sum(dia.valid) / n.valid
      acc.test = sum(dia.test) / n.test
      
      # compute expected accuracy for each data set
      exp.acc.train = sum(p.train * q.train)
      exp.acc.valid = sum(p.valid * q.valid)
      exp.acc.test = sum(p.test * q.test)
      
      # compute kappa for each data set
      kap.train = (acc.train - exp.acc.train) / (1 - exp.acc.train)
      kap.valid = (acc.valid - exp.acc.valid) / (1 - exp.acc.valid)
      kap.test = (acc.test - exp.acc.test) / (1 - exp.acc.test)
      
      # ---- compute one-vs-all metrics ----
      
      # compute a binary confusion matrix for each class, for each data set
      one.v.all.train = lapply(1:nrow(conf.train), function(i)
      {
        # get the four entries of a binary confusion matrix
        v = c(conf.train[i,i], 
              rsum.train[i] - conf.train[i,i], 
              csum.train[i] - conf.train[i,i], 
              n.train - rsum.train[i] - csum.train[i] + conf.train[i,i]);
        
        # build the confusion matrix
        return(matrix(v, nrow = 2, byrow = TRUE))
      })
      
      one.v.all.valid = lapply(1:nrow(conf.valid), function(i)
      {
        # get the four entries of a binary confusion matrix
        v = c(conf.valid[i,i], 
              rsum.valid[i] - conf.valid[i,i], 
              csum.valid[i] - conf.valid[i,i], 
              n.valid - rsum.valid[i] - csum.valid[i] + conf.valid[i,i]);
        
        # build the confusion matrix
        return(matrix(v, nrow = 2, byrow = TRUE))
      })
      
      one.v.all.test = lapply(1:nrow(conf.test), function(i)
      {
        # get the four entries of a binary confusion matrix
        v = c(conf.test[i,i], 
              rsum.test[i] - conf.test[i,i], 
              csum.test[i] - conf.test[i,i], 
              n.test - rsum.test[i] - csum.test[i] + conf.test[i,i]);
        
        # build the confusion matrix
        return(matrix(v, nrow = 2, byrow = TRUE))
      })
      
      # sum up all of the matrices for each data set
      one.v.all.train = Reduce('+', one.v.all.train)
      one.v.all.valid = Reduce('+', one.v.all.valid)
      one.v.all.test = Reduce('+', one.v.all.test)
      
      # compute the micro average accuracy for each data set
      micro.acc.train = sum(diag(one.v.all.train)) / sum(one.v.all.train)
      micro.acc.valid = sum(diag(one.v.all.valid)) / sum(one.v.all.valid)
      micro.acc.test = sum(diag(one.v.all.test)) / sum(one.v.all.test)
      
      # get the macro accuracy for each data set
      macro.acc.train = acc.train
      macro.acc.valid = acc.valid
      macro.acc.test = acc.test
      
      # ---- finalize output ----
      
      # build a final metrics table for stack.mod
      stack.mod.table = data.table(Model = "Super Learner",
                                   Metric = c("Log_Loss", "Kappa", "Macro_Accuracy", "Micro_Accuracy"),
                                   Train = c(mll.train, kap.train, macro.acc.train, micro.acc.train),
                                   Valid = c(mll.valid, kap.valid, macro.acc.valid, micro.acc.valid),
                                   Test = c(mll.test, kap.test, macro.acc.test, micro.acc.test))
      
      # update the row and column names of the confusion matrices
      colnames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
      rownames(conf.train) = levels(as.data.frame(ytrue.train)[,1])
      
      colnames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
      rownames(conf.valid) = levels(as.data.frame(ytrue.valid)[,1])
      
      colnames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
      rownames(conf.test) = levels(as.data.frame(ytrue.test)[,1])
      
      # build a list of confusion matrices
      stack.confusion = list("Train" = conf.train, 
                             "Valid" = conf.valid, 
                             "Test" = conf.test)
      
      # evaluate model pestackormance
      stack.mod.table
      stack.confusion
      
      # save the model
      h2o.saveModel(object = stack.mod, path = getwd())
      
      # load the model
      # stack.mod = h2o.loadModel(path = paste0(getwd(), "/stack.mod"))
      
      # remove some objects
      rm(dia.train, ynew.train.code, ytrue.train.code, ytrue.train.mat, 
         dia.valid, ynew.valid.code, ytrue.valid.code, ytrue.valid.mat, 
         dia.test, ynew.test.code, ytrue.test.code, ytrue.test.mat, 
         p.train, q.train, acc.train, exp.acc.train,
         p.valid, q.valid, acc.valid, exp.acc.valid,
         p.test, q.test, acc.test, exp.acc.test,
         conf.train, rsum.train, csum.train, n.train, 
         conf.valid, rsum.valid, csum.valid, n.valid, 
         conf.test, rsum.test, csum.test, n.test,
         one.v.all.train, one.v.all.valid, one.v.all.test,
         mll.train, kap.train, macro.acc.train, micro.acc.train,
         mll.valid, kap.valid, macro.acc.valid, micro.acc.valid,
         mll.test, kap.test, macro.acc.test, micro.acc.test,
         ynew.test, ynew.train, ynew.valid, ytrue.test, ytrue.train, ytrue.valid)
      
      # free up RAM
      gc()
      
    }
  }
}

# clean up the data in the h2o cluster
h2o.removeAll()

# shut down the h2o cluster
h2o.shutdown(prompt = FALSE)

# free up RAM
gc()

}









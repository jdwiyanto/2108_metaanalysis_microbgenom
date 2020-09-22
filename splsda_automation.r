#automation of PLSDA and SPLSDA analysis

########################
#set variables manually
########################

doc.title <- 'classification based on physeq.genus and mainland-immigrant model'
ds <-  physeq.genus # choose correct phyloseq object 
var <- 'eth' # variable of interest based on colnames of ds object
no.var <- sample_data(ds)[[var]] %>% factor %>% levels %>% length # set number of variables of var 
no.comp <- 20 # set number of components to test in plsda and splsda model
no.repeat <- 50 # set number of cross validation repeats in plsda and splsda model

########################
# set manually training and validation dataset below
########################

path <- 'E:/ubuntu/rev/'
setwd(path)

package <- c('phyloseq', 'ggplot2', 'dplyr', 'plyr', 
             'ape', 'propr', 'ggsci','ggpubr',
             'mixOmics', 'cowplot', 'rlang')

lapply(package, library, character.only = T)

# create directory

dir.create(paste0(path, 'splsda_', var, '/'))
setwd(paste0(path, 'splsda_', var, '/'))

# create otu table 

x <- otu_table(ds) %>% as.data.frame() %>% t() # create otu table from selected dataset
keep <- colSums(x) > 100 # exclude otu count lower than 100
x2 <- x[,keep] 
x3 <- x2 + 1 # apply pseudo-count of 1

###############################################################
#manually set how to divide the training and validation dataset
###############################################################

# option 1 - split ds 70:30 for training and validation

#set.seed(5)
#train <- sample(nrow(sample_data(ds)), 0.7*nrow(sample_data(ds)), replace = F)
#trainset <- sample_data(ds)[train,]
#validset <- sample_data(ds)[-train,]

# option 2 - split ds into immigrant and mainland communities

set.seed(5)
trainset <- subset_samples(ds, country2 != 'malaysia') %>% sample_data
validset <- subset_samples(ds, country2 == 'malaysia') %>% sample_data


# define classification variable

yt <- trainset[[var]]
yv <- validset[[var]]

# define otu table for training and validation dataset

xt <- subset(x3, rownames(x3) %in% rownames(trainset))
xv <- subset(x3, rownames(x3) %in% rownames(validset))

# plsda 

xt.clr <- logratio.transfo(xt, logratio = 'CLR') # transform abundance table with CLR

plsda.t <- plsda(xt.clr, yt, ncomp = no.comp) # create plsda object
plsda.t.perf <- perf(plsda.t, validation = 'Mfold', folds = 5, nrepeat = no.repeat, auc = T, progressBar = T) # cross-validation

## plot balanced error rate based on plsda model

er <- plsda.t.perf$error.rate$BER[,1]
component <- 1:no.comp
er.df <- data.frame(component, er)

fig.plsda.er <- ggplot(er.df, aes(x = component, y = er)) + 
  geom_point() +
  geom_smooth(linetype = 'dashed', color = 'black', size = 0.1, se = F) +
  scale_x_discrete(limits = component) +
  theme_bw() +
  xlab('No. components') +
  ylab('error rate')

## plot AUC curve for plsda model

best.comp <- plsda.t.perf$choice.ncomp[2,1] # take best no of components based on BER and max distance

plsda.auc <- auroc(plsda.t, roc.comp = best.comp)
fig.plsda.auc <- paste0('plsda.auc$graph.Comp', best.comp) %>% parse_expr %>% eval + theme_bw() + theme(legend.position = 'bottom') + guides(color = guide_legend(nrow = 6)) #AUC plot

## plot plsda ordination for the first 2 axes

plsda.plot <- plotIndiv(plsda.t, comp = c(1,2),
                        style = 'ggplot2', 
                        group = yt, 
                        ind.names = F, 
                        ellipse = T, 
                        legend = T,
                        col.per.group = color.mixo(1:no.var),
                        pch = 20,
                        title = 'PLS-DA')

fig.plsda.plot <- plsda.plot$graph + theme_bw() + theme(legend.position = 'bottom')

# validation test via confusion matrix

xv.clr <- logratio.transfo(xv, logratio = 'CLR') # transform abundance table with CLR

plsda.test <- predict(plsda.t, newdata = xv.clr, method = 'max.dist')
plsda.test2 <- plsda.test$class$max.dist %>% as.data.frame
table(yv, plsda.test2[,best.comp]) %>% prop.table(margin = 1)

# print error rate to document

sink('output.txt', append = F)

doc.title
ds
var
no.var
no.comp
no.repeat

#error rate for plsda model
plsda.t.perf$error.rate

#best number of component based on t-test
best.comp

#plsda confusion matrix based on best no of component
table(yv, plsda.test2[,best.comp])

#plsda confusion matrix in proportion
table(yv, plsda.test2[,best.comp]) %>% prop.table(margin = 1)

sink()

### splsda

splsda.t <- tune.splsda(xt.clr, 						# tune splsda
            Y = yt, 
            ncomp = no.comp, 
            multilevel = NULL, 
            test.keepX = c(seq(5,150, 5)), 
            validation = c('Mfold'), 
            folds = 5, 
            dist = 'max.dist', 
            nrepeat = no.repeat,
            cpus = 8,
            auc = T,
            progressBar = T)
			
best.comp2 <- splsda.t$choice.ncomp$ncomp				# see no. optimal component based on t test
choice.keepx <- splsda.t$choice.keepX[1:best.comp2]		# use the optimal number of components based on best.comp2

splsda.t2 <- splsda(X = xt.clr, Y = yt, ncomp = ifelse(best.comp2 < 2, 2, best.comp2), keepX = choice.keepx) # run splsda

splsda.auc <- auroc(splsda.t2, roc.comp = best.comp2)
fig.splsda.auc <- paste0('splsda.auc$graph.Comp', best.comp2) %>% parse_expr %>% eval + theme_bw() + theme(legend.position = 'bottom') + guides(color = guide_legend(nrow = 6)) #AUC plot

splsda.plot <- plotIndiv(splsda.t2, 
          ind.names = F, 
          col.per.group = color.mixo(1:no.var), 
          comp = c(1,2),
          pch = 20, 
          ellipse = TRUE,
          legend = TRUE,
          star = F,
          centroid = F,
          title = 'sPLS-DA')
		  
fig.splsda.plot <- splsda.plot$graph + theme_bw() + theme(legend.position = 'bottom')

# splsda cross validation

splsda.cv <- perf(splsda.t2, validation = 'Mfold', folds = 5, progressBar = T, nrepeat = no.repeat)

splsda.er <- splsda.cv$error.rate$BER[,1] # splsda error rate based on BER
splsda.comp <- 1:ifelse(best.comp2 < 2, 2, best.comp2)
splsda.er.df <- data.frame(splsda.comp, splsda.er)

fig.splsda.er <- ggplot(splsda.er.df, aes(x = splsda.comp, y = splsda.er)) + 
  geom_point() +
  geom_smooth(linetype = 'dashed', color = 'black', size = 0.1, se = F) +
  scale_x_discrete(limits = splsda.comp) +
  theme_bw() +
  xlab('No. components') +
  ylab('error rate')

# validation test via confusion matrix
splsda.test <- predict(splsda.t2, newdata = xv.clr, method = 'max.dist', dist)
splsda.test2 <- splsda.test$class$max.dist %>% as.data.frame
table(yv, splsda.test2[,best.comp2])

# sink splsda profile

sink('output.txt', append = T)

#error rate for splsda model
splsda.cv$error.rate

#best no of component based on splsda model
best.comp2

#confusion matrix based on splsda model
table(yv, splsda.test2[,best.comp2])

#splsda confusion matrix in proportion
table(yv, splsda.test2[,best.comp2]) %>% prop.table(margin = 1)

sink()

# compile plots
plot_grid(fig.plsda.plot, fig.splsda.plot, fig.plsda.er, fig.splsda.er, 
  fig.plsda.auc, fig.splsda.auc, labels = c('a', 'd', 'b', 'e', 'c', 'f'), nrow = 3, rel_heights = c(1.5,1,2))

ggsave('fig_all.pdf', units = 'mm', height = 250, width = 250)

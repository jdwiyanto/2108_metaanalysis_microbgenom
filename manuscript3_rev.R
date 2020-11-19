# packages used in this script
package <- c('phyloseq', 'ggplot2', 'dplyr', 'plyr', 
             'ape', 'propr', 'ggsci','ggpubr',
             'mixOmics', 'cowplot', 'tidyr')
lapply(package, library, character.only = T)

# set path

path <- 'E:/ubuntu/rev/val/'   

#################################################
##### IMPORT DATA TO CREATE PHYLOSEQ OBJECT #####
#################################################
require(stringr)

# import edge abundance table
edge <- read.csv('200718_rev.edge_tally.csv', header = T, row.names = 1)
rownames(edge) <- gsub('.exp.', '', rownames(edge))                                  # remove '.exp' from rownames


# import metadata 
meta <- read.csv('metadata.csv', header = T, row.names = 1)
meta$country2 <- meta$country2 %>% str_to_title()                                    # capitalise country
meta$eth <- meta$eth %>% str_to_title()                                              # capitalise ethnic group
meta$ce <- paste0(meta$eth, '.',meta$country2)                                       # merge ethnicity and country to new var


# import taxa table 
taxa <- read.csv('200718_rev.taxon_map.csv', header = T, row.names = 1) %>% as.matrix()


# clean data 
edge <- subset(edge, rownames(edge) %in% rownames(meta))                             # filter edge to remove jakun group (meta file contains final participant list)
edge <- edge[(order(match(rownames(edge), rownames(meta)))),]                        # sort edge data so participants in same order as metadata
colnames(edge) <- gsub('X', '', colnames(edge))                                      # remove 'X" from edge number due to transposition
edge[is.na(edge)] <- 0                                                               # convert NA to 0

ordering <- match(colnames(edge), rownames(taxa))
edge <- edge[,order(ordering)]
colnames(edge) <- paste0(colnames(edge), '_', taxa[,5])

rownames(taxa) <- paste0(rownames(taxa), '_', taxa[,5])

# create phyloseq object
physeq <- phyloseq(edge %>% otu_table(taxa_are_rows = F),
                   taxa %>% tax_table(),
                   meta %>% sample_data())


physeq.genus <- physeq %>% tax_glom(taxrank = 'genus')


# set variables

ps <- physeq.genus                          # main phyloseq object
var <- 'ce'                                 # single var to test
var2 <- c('country2', 'eth', 'ce')          # loop var to test



#############################
# ALPHA-DIVERSITY MEASUREMENT
#############################
library(vegan)

# generate rarefy curve
ps.rcurve <- subset_samples(ps, depth < 30000)                                     # remove sample with depth > 30k to beautify graph     
rarecurve(otu_table(ps.rcurve), step = 500, label = F)
dev.print(tiff, 'fig.rarefaction.tiff', units = 'mm', height = 200, width = 200, res = 300)            

ps.rare <- rarefy_even_depth(ps,                                                   # rarefy all samples to 10k depth
                             sample.size = 10000, 
                             replace = F, 
                             trimOTUs = T, 
                             verbose = T)

alpha <- estimate_richness(ps.rare, measures = c('Chao1', 'Shannon'))

for (x in var2) {

fig.alpha <- plot_richness(ps.rare, 
                           x = x,
                           measures = c('Chao1', 'Shannon')) +
                           geom_violin(fill = 'white', alpha = 0.5) +
                           geom_boxplot(aes(fill = 'blue'), alpha = 1) +
                           scale_fill_d3() +
                           theme_bw() +
                           theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                                 legend.position = 'none') +
                           xlab(element_blank()) 

assign(paste0('fig.alpha.', x), fig.alpha)

}

list.alpha <- ls(pattern = 'fig.alpha.')
ggarrange(plotlist = mget(list.alpha),
          ncol = 1,
          labels = 'auto')

ggsave('fig.alpha.pdf', units = 'mm', width = 200, height = 300)


# calculate significance across diversity

sink('alpha_sig.txt', append = F)

for (x in var2) {
  
  for (y in c('Chao1', 'Shannon')) {
  
  paste('testing variable', x, 'against', y) %>% print    
  sig.test <- kruskal.test(alpha[[y]], sample_data(ps.rare)[[x]])
  sig.test %>% print 

  }
}

sink()

require(FSA)

for (x in var2) {
  
  for (y in c('Chao1', 'Shannon')) {
    
  stat.posthoc <- dunnTest(alpha[[y]], sample_data(ps.rare)[[x]], method = 'bh') 
  stat.posthoc <- stat.posthoc$res 
  stat.posthoc$stat <-  y
  stat.posthoc$variable <- x
  
  assign(paste0('stat.posthoc.', x, '.', y), stat.posthoc)
  
  }
}

list.ph <- ls(pattern = 'stat.posthoc.')

df.posthoc <- bind_rows(mget(list.ph))

write.table(df.posthoc , file = 'df.posthoc.csv', quote = F, sep = ',', row.names = F, col.names = T)

################################
##### EXPLORATORY ANALYSIS #####
################################
require(microbiome)

ps.core <- core(ps, detection = 0.1, prevalence = 50/100, include.lowest = T)

# transform with centred-log transformation and recreate phyloseq object

edge.propr <- otu_table(ps.core)
edge.propr2 <- edge.propr %>% as.data.frame() %>% t() %>% propr()
edge.propr3 <- edge.propr2@logratio %>% t() %>% otu_table(taxa_are_rows = F)

ps.clr <- phyloseq(edge.propr3, sample_data(ps.core), tax_table(ps.core))

############
# Ordination
############

ord <- ordinate(ps.clr, method = 'PCoA', distance = 'euclidean')

for (x in var2) {

fig.pca <- plot_ordination(ps.clr, 
                type = "samples",
                ordination = ord,
                color = x) +
                geom_point(size = 2.5) +
                theme_bw() +
                theme(legend.position = 'bottom',legend.title = element_blank()) +
                scale_color_d3() +
                guides(color = guide_legend(nrow = 3))

assign(paste0('fig.pca.', x), fig.pca)

}

list.pca <- ls(pattern = 'fig.pca.')

ggarrange(plotlist = mget(list.pca),
          nrow = 1,
          labels = 'auto')

ggsave('fig.pca.pdf', units = 'mm', height = 300, width = 400)
                

###########
# permanova
###########
library(vegan)

set.seed(123)

sample_data(ps.core)$names <- sample_names(ps.core)                           
grouping <- sample_data(ps.core) %>% group_by(country2) %>% sample_n(70)            # randomly sample 70 from each country for permanova

retained <- sample_names(ps.core) %in% grouping$names
ps.perm <- subset_samples(ps.core, retained)                                        # take ps.core taxa

otu.perm <- otu_table(ps.perm) %>% data.frame
meta.perm <- sample_data(ps.perm) %>% data.frame

sink('permanova.txt', append = F)

for (x in var2) {
  
  paste('test permanova with ps.core with', x) %>% print
  
  permanova <- adonis(otu.perm ~ x %>% get,  
                    data = meta.perm,
                    permutations = 999,
                    method = "euclidean") 

  print(permanova)

}

sink()

# DRM model
require(DirichletMultinomial); require(lattice); require(tidyr)

d.otu <- otu_table(ps.core)                                          # set otu table and add pseudo-count of 1

cnts <- log10(colSums(d.otu))                                              # optional, confirm density plot is normal (neg binomial distribution)
densityplot(cnts, xlim = range(cnts))

fit <- mclapply(1:5, dmn, count = d.otu, verbose = T)
fit

lplc <- sapply(fit, laplace)
plot(lplc, type = "b", xlab = "number of Dirichlet component", ylab = "Model fit")

best <- fit[[which.min(lplc)]]
best

d.type <- mixture(best)                                                         # give cluster prediction for each sample
colnames(d.type) <- c(1:which.min(lplc))                                        # give column names to cluster prediction
cluster <- colnames(d.type)[apply(d.type, 1, which.max)]                        # assign new column with chosen best prediciton
d.type2 <- cbind(d.type, cluster) %>% as.data.frame                             # combine best cluster column with cluster prediction data frame

p0 <- fitted(fit[[1]], scale=TRUE)                                              # scale by theta
p3 <- fitted(best, scale=TRUE)
colnames(p3) <- paste("m", 1:which.min(lplc), sep="")

diff <- rowSums(abs(p3 - as.vector(p0)))                                        # generate dataframe of enterotypes
o <- order(diff, decreasing=TRUE)
cdiff <- cumsum(diff[o]) / sum(diff)
df.drm <- ( head(cbind(Mean=p0[o], p3[o,], diff=diff[o], cdiff), 30) ) %>% as.data.frame()

write.csv(df.drm, 'dirichlet.csv')

df.drm$names <- rownames(df.drm)
df.drm2 <- df.drm[1:5,c((1:which.min(lplc) +1), (which.min(lplc) + 4))] %>% 
  pivot_longer(!names, names_to = 'enterotype', values_to = 'prop')

fig.drm <- ggplot(df.drm2, aes( x = enterotype,
                                y = prop*100,
                                fill = names)) +
  geom_col() +
  theme_bw() +
  ylab('Proportion (%) explained by top 5 genera') +
  theme(legend.title = element_blank())

ggsave('fig.drm.genus.pdf', units = 'mm', width = 400, height = 300)

for (x in var2) {                                                                # loop of enterotype - variable association
  
  drm.table <- table(d.type2$cluster, meta[[x]]) %>% prop.table(margin = 1) %>% data.frame
  
  fig.drm.prev <- ggplot(drm.table, aes(x = Var1,
                                        y = Freq*100,
                                        fill = Var2)) + 
                  geom_col(position = 'stack') +
                  scale_fill_discrete(name = x) +
                  theme_bw() +
                  ylab('Proportion (%)') +
                  xlab('Enterotype') + 
                  theme(legend.position = 'bottom',legend.title = element_blank()) +
                  guides(fill = guide_legend(nrow = 3))
                  
  
  assign(paste0('fig.drm.prev.', x), fig.drm.prev)
  
}

fig.drm.list <- paste0('fig.drm.prev.', var2) 
fig.drm2 <- ggarrange(plotlist = mget(fig.drm.list), labels = c('b', 'c', 'd'))

ggarrange(fig.drm, fig.drm2, labels = c('a', ''))

ggsave('fig.drm.compiled.pdf', units = 'mm', width = 400, height = 300)


#######################
# bifidobacterium count
#######################

otu.clr <- otu_table(ps.clr) %>% data.frame
otu.core <- otu_table(ps.core) %>% data.frame

for (x in var2) {

bifido.clr <- tapply(otu.clr$X5723, sample_data(ps.clr)[[x]], mean)
bifido.raw <- tapply(otu.core$X5723, sample_data(ps.core)[[x]], mean)

bifido.all <- rbind(bifido.raw, bifido.clr) %>% as.data.frame %>% t %>% as.data.frame
bifido.all$name <- rownames(bifido.all)

bifido.all2 <- bifido.all %>% pivot_longer(-name, names_to = 'type', values_to = 'abundance')

fig.bif <- ggplot(bifido.all2, aes(x = name, y = abundance %>% log, group = type)) + 
  geom_point(aes(shape = type)) +
  geom_line(aes(linetype = type)) +
  theme_bw() +
  xlab(element_blank()) +
  ylab('log abundance')

assign(paste0('fig.bifido.', x), fig.bif)

}

list.bifido <- ls(pattern = 'fig.bifido.')
ggarrange(plotlist = mget(list.bifido), labels = 'auto', ncol = 1)

ggsave('fig.bifido.pdf', units = 'mm', height = 200, width = 300)

save.image(file = 'r_env_manuscript3.RData')

##### SPLSA #####

for (var in var2) { 
  
  for (valmodel in c('malaysia', 'all')) {                                 # determine validation model, malaysia (ONLY WITH ETH VAR)
    # or all (valid for ce, eth and country2)
    
    ds <-  ps.core                                                         # choose correct phyloseq object 
    no.comp <- 20                                                           # set number of components to test in plsda and splsda model
    no.repeat <- 50                                                         # set number of cross validation repeats in plsda and splsda model
    
    
    path2 <- paste0('splsda_', var, '_', valmodel)
    
    setwd(path)
    dir.create(path2)
    setwd(paste0(path, path2))
    
    
    
    x3 <- ( otu_table(ds) %>% as.data.frame() ) + 1                        # create otu table 
    
    set.seed(5)
    
    train <- sample(nrow(sample_data(ds)), 0.7*nrow(sample_data(ds)), replace = F) # option 1 - split ds 70:30 for training and validation
    trainset.all <- sample_data(ds)[train,]
    validset.all <- sample_data(ds)[-train,]
    
    trainset.msia <- subset_samples(ds, country2 != 'Malaysia') %>% sample_data  # option 2 - split ds into immigrant and mainland communities
    validset.msia <- subset_samples(ds, country2 == 'Malaysia') %>% sample_data
    
    trainset.call <- ifelse(valmodel == 'all', 'trainset.all', ifelse(valmodel == 'malaysia', 'trainset.msia', 'none defined'))
    validset.call <- ifelse(valmodel == 'all', 'validset.all', ifelse(valmodel == 'malaysia', 'validset.msia', 'none defined'))
    
    trainset <- trainset.call %>% get
    validset <- validset.call %>% get
    
    yt <- trainset[[var]]                 # var for training set
    yv <- validset[[var]]                 # for for validation set
    
    no.var <- yt %>% unique %>% length
    
    xt <- subset(x3, rownames(x3) %in% rownames(trainset))    # otu for training set
    xv <- subset(x3, rownames(x3) %in% rownames(validset))    # otu for validation set
    
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
    
    best.comp.pre <- plsda.t.perf$choice.ncomp[2,1]             # take best no of components based on BER and max distance
    best.comp <- ifelse(best.comp.pre < 2, 2, best.comp.pre)    # make sure 2 is the minimum component no for plotting
    
    plsda.auc <- auroc(plsda.t, roc.comp = best.comp)
    
    fig.plsda.auc <- paste0('plsda.auc$graph.Comp', best.comp) %>% 
      parse_expr %>% eval + 
      theme_bw() + 
      theme(legend.position = 'bottom') + 
      guides(color = guide_legend(nrow = 6)) #AUC plot
    
    ## plot plsda ordination for the first 2 axes
    
    plsda.plot <- plotIndiv(plsda.t, comp = c(1,2),
                            style = 'ggplot2', 
                            group = yt, 
                            ind.names = F, 
                            ellipse = F, 
                            legend = T,
                            title = 'PLS-DA')
    
    fig.plsda.plot <- plsda.plot$graph + theme_bw() + theme(legend.position = 'bottom')
    
    # validation test via confusion matrix
    
    xv.clr <- logratio.transfo(xv, logratio = 'CLR') # transform abundance table with CLR
    
    plsda.test <- predict(plsda.t, newdata = xv.clr, method = 'max.dist')
    plsda.test2 <- plsda.test$class$max.dist %>% as.data.frame
    table(yv, plsda.test2[,best.comp]) %>% prop.table(margin = 1)
    
    # print error rate to document
    
    sink('output2.txt', append = F)
    
    paste('PLSDA and SPLSDA analysis for', var) %>% print
    
    paste('physeq object used is') %>% print
    ds %>% print
    
    paste('this test used', no.comp, 'number of components, repeated for', no.repeat, 'times') %>% print
    
    paste('error rate for plsda model') %>% print
    plsda.t.perf$error.rate %>% print
    
    paste('best number of component based on t-test') %>% print
    best.comp %>% print
    
    paste('plsda confusion matrix based on best no of component') %>% print
    table(yv, plsda.test2[,best.comp]) %>% print
    
    paste('plsda confusion matrix in proportion') %>% print
    table(yv, plsda.test2[,best.comp]) %>% prop.table(margin = 1) %>% print
    
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
    
    best.comp2.pre <- splsda.t$choice.ncomp$ncomp				            # see no. optimal component based on t test
    best.comp2 <- ifelse(best.comp2.pre < 2, 2, best.comp2.pre)     # ensure best comp > 1
    choice.keepx <- splsda.t$choice.keepX[1:best.comp2]		          # use the optimal number of components based on best.comp2
    
    splsda.t2 <- splsda(X = xt.clr, Y = yt, ncomp = best.comp2, keepX = choice.keepx) # run splsda
    
    splsda.auc <- auroc(splsda.t2, roc.comp = best.comp2)
    fig.splsda.auc <- paste0('splsda.auc$graph.Comp', best.comp2) %>% parse_expr %>% eval + theme_bw() + theme(legend.position = 'bottom') + guides(color = guide_legend(nrow = 6)) #AUC plot
    
    splsda.plot <- plotIndiv(splsda.t2, 
                             ind.names = F, 
                             col.per.group = color.mixo(1:no.var), 
                             comp = c(1,2),
                             ellipse = F,
                             legend = TRUE,
                             star = F,
                             centroid = F,
                             title = 'sPLS-DA')
    
    fig.splsda.plot <- splsda.plot$graph + theme_bw() + theme(legend.position = 'bottom')
    
    # splsda cross validation
    
    splsda.cv <- perf(splsda.t2, validation = 'Mfold', folds = 5, progressBar = T, nrepeat = no.repeat)
    
    splsda.er <- splsda.cv$error.rate$BER[,1] # splsda error rate based on BER
    splsda.comp <- 1:best.comp2
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
    
    sink('output2.txt', append = T)
    
    paste('error rate for splsda model') %>% print
    splsda.cv$error.rate %>% print
    
    paste('best no of component based on splsda model') %>% print 
    best.comp2 %>% print
    
    paste('confusion matrix based on splsda model') %>% print
    table(yv, splsda.test2[,best.comp2]) %>% print
    
    paste('splsda confusion matrix in proportion') %>% print()
    table(yv, splsda.test2[,best.comp2]) %>% prop.table(margin = 1) %>% print
    
    sink()
    
    # compile plots
    
    figcomp.plsda1 <- ggarrange(fig.plsda.er, fig.plsda.auc, labels = 'auto')
    figcomp.plsda2 <- ggarrange(figcomp.plsda1, fig.plsda.plot, ncol = 1, labels = c('', 'c'))
    ggsave('fig_plsda.pdf', figcomp.plsda2, units = 'mm', height = 250, width = 250)
    
    figcomp.splsda1 <- ggarrange(fig.splsda.er, fig.splsda.auc, labels = 'auto')
    figcomp.splsda2 <- ggarrange(figcomp.splsda1, fig.splsda.plot, ncol = 1, labels = c('', 'c'))
    ggsave('fig_splsda.pdf', figcomp.splsda2, units = 'mm', height = 250, width = 250)
    
    plot_grid(fig.plsda.plot, fig.splsda.plot, fig.plsda.er, fig.splsda.er, 
              fig.plsda.auc, fig.splsda.auc, labels = c('auto'), nrow = 3, rel_heights = c(1.5,1,2))
    
    ggsave('fig_all.pdf', units = 'mm', height = 250, width = 250)
    
    
    # contribution plot -- mainly for mainland-immigrant model output to see mutual taxa differentiating both groups
    
    if (valmodel == 'malaysia') {
      if (var == 'country2') {
        
        save.image(file = paste0('splsda_', var, '_', valmodel, '.RData'))
        
        setwd(path)
        
      }
    } else {
      
    splsda.v <- tune.splsda(xv.clr, 						# tune splsda
                            Y = yv, 
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
    
    
    
    best.comp3.pre <- splsda.v$choice.ncomp$ncomp				            # see no. optimal component based on t test
    best.comp3 <- ifelse(best.comp3.pre < 2, 2, best.comp3.pre)     # ensure best comp > 1
    choice.keepx2 <- splsda.v$choice.keepX[1:best.comp3]		          # use the optimal number of components based on best.comp2
    
    splsda.v2 <- splsda(X = xv.clr, Y = yv, ncomp = best.comp3, keepX = choice.keepx2) # run splsda
    
    contib.t <- lapply(1:best.comp2, function(x){                   # generate contrib dataframe for training set
      plotLoadings(splsda.t2, 
                   comp = x, 
                   method = 'mean', 
                   contrib = 'max')})
    
    contib.t2 <- lapply(1:best.comp2, function(x) cbind(contib.t[[x]], x, contib.t[[x]] %>% rownames))
    contib.t3 <- do.call(rbind, contib.t2)
    contib.t4 <- contib.t3 %>% arrange(desc(abs(importance))) %>% group_by(x) %>% slice(1:7)
    contib.t4$genus <- contib.t4[['contib.t[[x]] %>% rownames']]
    contib.t4$group = 'training'
    
    contib.v <- lapply(1:best.comp3, function(x){                   # generate contrib dataframe for validation set
      plotLoadings(splsda.v2, 
                   comp = x, 
                   method = 'mean', 
                   contrib = 'max')})
    
    contib.v2 <- lapply(1:best.comp3, function(x) cbind(contib.v[[x]], x, contib.v[[x]] %>% rownames))
    contib.v3 <- do.call(rbind, contib.v2)
    contib.v4 <- contib.v3 %>% arrange(desc(abs(importance))) %>% group_by(x) %>% slice(1:7)
    contib.v4$genus <- contib.v4[['contib.v[[x]] %>% rownames']]
    contib.v4$group = 'validation'
    
    contib.all <- rbind(contib.t4, contib.v4)                       # combine both training and validation contrib dataframe
    
    fig.contib <- ggplot(contib.all, aes(x = genus, y = importance %>% abs, fill = GroupContrib %>% as.factor )) + 
      geom_col(position = 'dodge') + 
      coord_flip() + 
      facet_grid(~paste0(group)) +
      guides(color = guide_legend(nrow = 1)) +
      theme_bw() +
      labs(fill = 'component no.') +
      theme(legend.position = 'bottom', legend.title = element_blank()) +
      xlab(element_blank()) +
      ylab('importance') +
      scale_fill_d3()
    
    ggsave('fig_contrib.pdf', units = 'mm', height = 250, width = 250)
    
    save.image(file = paste0('splsda_', var, '_', valmodel, '.RData'))
    
    setwd(path)
    
    }
  }
}


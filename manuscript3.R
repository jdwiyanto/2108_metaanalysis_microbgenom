# packages used in this script
setwd('E:\\ubuntu\\rev')

package <- c('phyloseq', 'ggplot2', 'dplyr', 'plyr', 
             'ape', 'propr', 'ggsci','ggpubr',
             'mixOmics', 'cowplot', 'tidyr')
lapply(package, library, character.only = T)


#################################################
##### IMPORT DATA TO CREATE PHYLOSEQ OBJECT #####
#################################################

# import edge abundance table
edge <- read.csv('200718_rev.edge_tally.csv', header = T, row.names = 1)
rownames(edge) <- gsub('.exp.', '', rownames(edge)) # remove '.exp' from rownames

# import metadata 
meta <- read.csv('metadata.csv', header = T, row.names = 1)
meta$ce <- paste0(meta$eth, '.',meta$country2)

# filter edge to remove jakun group (meta file contains final participant list)
edge <- subset(edge, rownames(edge) %in% rownames(meta))

# sort edge data so participants in same order as metadata
edge <- edge[(order(match(rownames(edge), rownames(meta)))),]

# create phyloseq object for edge abundance table
edge <- t(edge) # transpose edge table to adhere with phyloseq orientation requirement
rownames(edge) <- gsub('X', '', rownames(edge)) # remove 'X" from edge number due to transposition

edge.p <- otu_table(edge, taxa_are_rows = T) # convert to phyloseq object
edge.p[is.na(edge.p)] <- 0 # convert NAs to 0

# import taxa table and transform to phyloseq object
taxa <- read.csv('200718_rev.taxon_map.csv', header = T, row.names = 1) %>% as.matrix()
taxa.p <- tax_table(taxa) # convert to phyloseq object

# transform metadata to a phyloseq object 
meta.p <- sample_data(meta)

# create phyloseq object
physeq <- phyloseq(edge.p, taxa.p, meta.p)
physeq

#############################
# ALPHA-DIVERSITY MEASUREMENT
#############################
library(vegan)

# generate rarefying curve based on dataset with depth > 100000 removed to beautify grapph
physeq.rarefaction <- subset_samples(physeq, depth < 30000)
rarecurve(otu_table(physeq.rarefaction) %>% t(), step = 500, label = F)
dev.print(tiff, 'fig.rarefaction.tiff', units = 'mm', height = 200, width = 200, res = 300)            

# rarefy to 10,000 depth reads
physeq.rare <- rarefy_even_depth(physeq, sample.size = 10000, replace = F, trimOTUs = T, verbose = T)

alpha <- estimate_richness(physeq.rare, measures = c('Chao1', 'Shannon'))

fig.alpha <- plot_richness(physeq.rare, 
                           x = 'ce',
                           measures = c('Chao1', 'Shannon')) +
  geom_violin(fill = 'white', alpha = 0.5) +
  geom_boxplot(aes(fill = ce), alpha = 1) +
  scale_fill_d3() +
  theme_bw() +
  theme(legend.position = 'none') +
  xlab(element_blank()) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

fig.alpha
fig.alpha.country
fig.alpha.eth

ggarrange(fig.alpha.country, fig.alpha.eth, fig.alpha, nrow = 3, ncol = 1, labels = c('a', 'b', 'c'))
ggsave('fig.alpha.tiff', units = 'mm', width = 200, height = 200, dpi = 300)

# calculate significance across diversity

kruskal.test(alpha$Chao1, sample_data(physeq.rare)$eth) # overall significance

library(FSA)

stat.shannon.ce <- dunnTest(alpha$Shannon, sample_data(physeq.rare)$ce, method = 'bh') # post-hoc for KW test

write.table(stat.chao.ce$res , file = 'stat.chao.ce.csv', quote = F, sep = ',', row.names = F, col.names = T)

################################
##### EXPLORATORY ANALYSIS #####
################################

# transform with centred-log transformation and recreate phyloseq object

edge.propr <- otu_table(physeq.sea)
keep <- rowSums(edge.propr) > 100 # retain taxa with minimum of 101 counts
edge.propr2 <- edge.propr[keep,]
edge.propr3 <- edge.propr2 %>% as.data.frame() %>% t() %>% propr()
edge.propr4 <- edge.propr3@logratio %>% t() %>% otu_table(taxa_are_rows = T)

physeq.sea.clr <- phyloseq(edge.propr4, sample_data(physeq.sea), tax_table(physeq.sea))

############
# Ordination
############

ord <- ordinate(physeq.clr, method = 'PCoA', distance = 'euclidean')

plot_ordination(physeq.clr, 
                type = "samples",
                ordination = ord,
                color = 'BioProject') +
                geom_point(size = 2.5) +
                theme_bw() +
                theme(legend.position = 'none',legend.title = element_blank()) +
                scale_color_d3() +
                guides(color = guide_legend(nrow = 1))
                
                
fig.pca.country
fig.pca.eth
fig.pca.ce

fig.pca.top <- ggarrange(fig.pca.country, fig.pca.eth, nrow = 1, ncol = 2, labels = c('a', 'b'))
ggarrange(fig.pca.top, fig.pca.ce, nrow = 2, labels = c('', 'c'))
ggsave('fig.pca2.tiff', units = 'mm', height = 200, width = 200, dpi = 300)

###########
# permanova
###########
library(vegan)

grouping <- sample_data(physeq) %>% group_by(country2)
grouping$name <- rownames(grouping)
grouping2 <- sample_n(grouping, 70)

grouping2$name
meta.perm <- meta[rownames(meta) %in% grouping2$name,]
physeq.perm <- subset_samples(physeq, sample_names(physeq) %in% grouping2$name)
meta.msia <- subset(meta, rownames(meta) %in% sample_names(physeq.msia))
meta.sea <- subset(meta, rownames(meta) %in% sample_names(physeq.sea))

permanova5 <- adonis(otu_table(physeq.clr) %>% t() ~ ce,
                    data = meta,
                    permutations = 999,
                    method = "euclidean") 

print(permanova2)

#permanova = physeq.clr - eth
#permanova2 = physeq.clr - country2
#Permanova3 = physeq.msia.clr - eth
#permanova4 = physeq.sea.clr ~ country2
#permanova5 = physeq.clr ~ ce

####################################
# contribution plot from splsda plot
####################################
# run this after splsda automation on mainland-immigrant model

splsda.t2 # the mainland model from the automated splsda script

# generate new splsda model for immigrant (malaysian) group

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

best.comp2.v <- splsda.v$choice.ncomp$ncomp				# see no. optimal component based on t test
choice.keepx.v <- splsda.v$choice.keepX[1:best.comp2.v]	

# the malaysian model
splsda.v2 <- splsda(X = xv.clr, Y = yv, ncomp = ifelse(best.comp2 < 2, 2, best.comp2.v), keepX = choice.keepx.v) # run splsda

# run this contib twice, once with mainland model (splsda.t2) and another with immigrant model (splsda.v2)

# splsda.t2 15 optimal components
# splsda.v2 3 optimal components

contib <- lapply(1:15, function(x){  # edit range based on no of optimal component
   plotLoadings(splsda.t2, 
             comp = x, 
             method = 'mean', 
             contrib = 'max')})

contib2 <- lapply(1:15, function(x) cbind(contib[[x]], x, contib[[x]] %>% rownames))
contib3 <- do.call(rbind, contib2)
contib4 <- contib3 %>% arrange(desc(abs(importance))) %>% group_by(x) %>% slice(1:7)
contib4$genus <- contib4$`contib[[x]] %>% rownames`

write.csv(contib4, 'contib_mainland.csv')

contib5.msia <- read.csv('contib_immigrant.csv', header = T) # after assigning genus name in excel
contib5 <- read.csv('contib_mainland.csv', header = T) # after assigning genus name in excel

contib6 <- contib5 %>% subset(contib5$GroupContrib != 'malay')
contib6.msia <- contib5.msia %>% subset(contib5.msia$GroupContrib != 'malay')

contib6$group <- 'mainland'
contib6.msia$group <- 'malaysia'

contib5$group <- 'mainland'
contib5.msia$group <- 'malaysia'

contib7 <- rbind(contib5, contib5.msia)

fig.contrib.malay <- ggplot(contib7, aes(x = genus2, y = importance %>% abs, fill = GroupContrib %>% as.factor )) + 
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
  
fig.contrib.all
fig.contrib.malay # same as fig.contrib.all but with malay

# network analysis
n <- mixOmics::network(res.splsda2, 
        comp = 1:9, 
        cutoff = 0.2,
        show.edge.labels = T)

dev.print(pdf, 'fig.network.msia.pdf')
dev.print(tiff, 'fig.network2.tiff', units = 'mm', height = 200, width = 200, res = 300)

#####################################
# DIRICHLET MULTINOMIAL MIXTURE MODEL
#####################################

library("DirichletMultinomial")
library(lattice)

# create genus-agglomerated count files from phyloseq object
physeq.msia <- subset_samples(physeq.genus, country2 == 'malaysia')
physeq.china <- subset_samples(physeq.genus, country2 == 'china')
physeq.india <- subset_samples(physeq.genus, country2 == 'india')
physeq.sea <- subset_samples(physeq.genus, country2 == 'malaysia' | country2 == 'indonesia')

d.otu <- otu_table(physeq.genus) %>% t()
d.otu <- d.otu + 1
colSums(d.otu) %>% summary()
keep <- colSums(d.otu) > 1108
d.otu <- d.otu[,keep]

cnts <- log10(colSums(d.otu))
densityplot(cnts, xlim = range(cnts))
cnts
fit <- mclapply(1:7, dmn, count = d.otu, verbose = T)
fit

lplc <- sapply(fit, laplace)
plot(lplc, type = "b", xlab = "number of Dirichlet component", ylab = "Model fit")

best <- fit[[which.min(lplc)]]
best

mixturewt(best)
head(mixture(best), 5)
splom(log(fitted(best)))

d.type <- mixture(best) # give cluster prediction for each sample
colnames(d.type) <- c("1", "2", '3', '4', '5', '6') # give column names to cluster prediction
cluster <- colnames(d.type)[apply(d.type, 1, which.max)] # assign new column with chosen best prediciton
d.type2 <- cbind(d.type, cluster) # combine best cluster column with cluster prediction data frame
d.type2 <- as.data.frame(d.type2)

p0 <- fitted(fit[[1]], scale=TRUE) # scale by theta
p3 <- fitted(best, scale=TRUE)
colnames(p3) <- paste("m", 1:6, sep="")
(meandiff <- colSums(abs(p3 - as.vector(p0))))

sum(meandiff)

diff <- rowSums(abs(p3 - as.vector(p0)))
o <- order(diff, decreasing=TRUE)
cdiff <- cumsum(diff[o]) / sum(diff)
df <- head(cbind(Mean=p0[o], p3[o,], diff=diff[o], cdiff), 30)
df

df2 <- as.data.frame(df)
write.csv(df2, "dirichlet_genus.csv")

B <- sample_data(physeq.genus)

B$enterotype <- d.type2$cluster

entero <- B$enterotype

dirichlet.prop <- table(entero, B$ce) %>% prop.table(margin = 2) %>% as.data.frame()
dirichlet.prop$entero2 <- dirichlet.prop$entero
dirichlet.prop$entero2[dirichlet.prop$entero2 == '5'] <- 1
dirichlet.prop$entero2

fig.dirichlet <- ggplot(dirichlet.prop, aes(x = Var2, 
                           y = Freq*100,
                           fill = entero2)) + 
  geom_col(position = 'stack') +
  scale_fill_d3(labels = c('1' = 'Bacteroides',
                           '2' = 'Bifidobacterium',
                           '3' = 'Prevotella',
                           '6' = 'Mixed')) + 
  theme_bw()  +
  ylab('Proportion (%)') +
  xlab(element_blank()) +
  theme(legend.title = element_blank()) 
  
  
ggsave('fig.dirichlet.tiff', units = 'mm', dpi = 'print', height = 150, width = 200)

#######################
# bifidobacterium count
#######################
bifido.clr2 <- tapply(X2.clr[,X2.clr = '5723'], Y3, mean)
bifido.raw2 <- tapply(X2[,X2 = '5723'], Y3, mean)
bifido.count2 <- rbind(bifido.raw2, bifido.clr2) %>% as.data.frame %>% t %>% as.data.frame
bifido.count2$name <- rownames(bifido.count2)
bifido.count3 <- bifido.count2 %>% pivot_longer(-name, names_to = 'type', values_to = 'abundance')

fig.bif <- ggplot(bifido.count3, aes(x = name, y = abundance %>% log, group = type)) + 
  geom_point(aes(shape = type)) +
  geom_line(aes(linetype = type)) +
  theme_bw() +
  xlab(element_blank()) +
  ylab('log abundance')

######################################
# complete figure for manuscript list
######################################

# figure 1
fig.alpha
fig.alpha.country
fig.alpha.eth

ggarrange(fig.alpha.country, fig.alpha.eth, fig.alpha, nrow = 3, ncol = 1, labels = c('a', 'b', 'c'))
ggsave('fig1.pdf', units = 'mm', width = 200, height = 200, dpi = 300)

# figure 2
plot_grid(fig.pca1$graph + guides(color = guide_legend(nrow = 6)), 
          fig.pca2$graph + guides(color = guide_legend(nrow = 6)),
          fig.pca3$graph + guides(color = guide_legend(nrow = 6)), 
          labels = c('a', 'b', 'c'), nrow =1)
fig.pca.compiled <- plot_grid(fig.pca.top, fig.pca3$graph, nrow = 2, labels = c('', 'c'))
ggsave('fig2_draft3.pdf', units = 'mm', dpi = 300, height = 200, width = 250)


# figure 3
fig.spl1 <- fig.splsda.er
fig.spl2 <- x.sp$graph.Comp13 + theme_bw() + guides(color = guide_legend(nrow = 3)) + theme(legend.position = 'bottom')
fig.spl3 <- y.sp$graph + theme_bw() + theme(legend.position = 'bottom')

fig.spl.top <- plot_grid(fig.spl1, fig.spl2, labels = c('a', 'b'))
fig.spl.compiled <- plot_grid(fig.spl.top, fig.spl3, nrow = 2, labels = c('', 'c'))
fig.spl.compiled
ggsave('fig3.pdf', units = 'mm', height = 200, width = 200, dpi = 300)

# figure 4
fig.contrib.all
ggsave('fig4.pdf', units = 'mm', height = 200, width = 200, dpi = 300)

# figure 5
fig.dirichlet 
ggsave('fig5.pdf', units = 'mm', dpi = 'print', height = 150, width = 200)

# suppl 2 figure rarefaction
physeq.rarefaction <- subset_samples(physeq, depth < 30000)
rarecurve(otu_table(physeq.rarefaction) %>% t(), step = 500, label = F)
dev.print(pdf, 'suppl_fig2.pdf')  

# suppl file 5 - PLS-DA
fig.plsda1 <- fig.plsda.er 
fig.plsda.2 <- x$graph.Comp14 + theme_bw() + theme(legend.position = 'bottom') + guides(color = guide_legend(nrow = 3))
fig.plsda3 <- y$graph + theme_bw()+ theme(legend.position = 'bottom')

plsda.bottom <- plot_grid(fig.plsda1, fig.plsda.2, labels = c('a','b')) 
fig.plsda.compiled <- plot_grid(plsda.bottom, fig.plsda3, ncol = 1, labels = c('', 'c'))
fig.plsda.compiled 
ggsave('suppl_fig5.pdf', units = 'mm', width = 200, height = 200, dpi = 300)

# suppl file 6 - sPLS-DA no malay
fig.spl2.1 <- fig.splsda.er2
fig.spl2.2 <- x.sp2$graph.Comp9 + theme_bw() + theme(legend.position = 'bottom')
fig.spl2.3 <- y.sp2$graph + theme_bw() + theme(legend.position = 'bottom')

fig.spl2.top <- plot_grid(fig.spl2.1, fig.spl2.2, labels = c('a', 'b'))
fig.spl2.compiled <- plot_grid(fig.spl2.top, fig.spl2.3, nrow = 2, labels = c('', 'c'))
fig.spl2.compiled 
ggsave('suppl_fig6.pdf', units = 'mm', height = 200, width = 200, dpi = 300)

# suppl file 7 - sPLS-DA malaysian
fig.spl.v1 <- fig.splsda.er.v
fig.spl.v2 <- x.sp.v$graph.Comp4 + theme_bw() + theme(legend.position = 'bottom')
fig.spl.v3 <- y.sp.v$graph + theme_bw() + theme(legend.position = 'bottom')

fig.spl.v.top <- plot_grid(fig.spl.v1, fig.spl.v2, labels = c('a', 'b'))
fig.spl.v.compiled <- plot_grid(fig.spl.v.top, fig.spl.v3, labels = c('', 'c'), nrow = 2)
fig.spl.v.compiled
ggsave('suppl_7.pdf', units = 'mm', height = 200, width = 200, dpi = 300)

# suppl 9 bifido
fig.bif 
ggsave('suppl_9_draft3.pdf', units = 'mm', height = 200, width = 200)

# new figure, splsda ce

fig.plsda.ce.compiled
fig.splsda.ce.compiled
ggsave('plsda_ce.pdf', units = 'mm', height = 300, width = 300)

# list of phyloseq object

physeq # raw phyloseq data
physeq.genus # physeq data agglomerated to genus level
physeq.cim # physeq.genus with only china, india and malaysia

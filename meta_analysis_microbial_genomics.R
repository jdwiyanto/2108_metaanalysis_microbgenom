setwd()
library(tidyverse)
library(phyloseq)
library(mbiome)
library(vegan)
library(ppcor)
library(corrplot)
library(LDM)

otu_table(physeq.genus) 

ps.rf = rarefy_even_depth(physeq.genus, sample.size = 10000)

ps.rf.clr = jd_clr(ps.rf)

# permanova test controlled for country
adonis(otu_table(ps.rf.clr) %>% t ~ eth, data = sample_data(ps.rf.clr) %>% data.frame, permutations = 999, method = 'euclidean')



# revised sequences

list.files()
asv = read.csv('phyloseq_otu.csv', header = T, row.names = 1)
tax = read.csv('phyloseq_tax.csv', header = T, row.names = 1) %>% as.matrix()

asv %>% rownames() %>% length

meta # from r_env_manuscript3

asv = asv[rownames(asv) %in% rownames(meta),] 
asv2 = asv[order(match(rownames(asv), rownames(meta))),]

ps.rev = phyloseq(otu_table(asv2, taxa_are_rows = F),
                  tax_table(tax),
                  sample_data(meta))

ps.rev.genus = ps.rev %>% tax_glom('Genus')

# rarefaction curve
rf.remove = (ps.rev.genus %>% otu_table() %>% data.frame %>% rowSums) < 30000 
ps.rev.genus.rare = subset_samples(ps.rev.genus, rf.remove)

pdf('fig_rarecurve.pdf', height = 7.5, width = 7.5)
rarecurve(otu_table(ps.rev.genus.rare) %>% data.frame, step = 100, label = F)
mtext(text = expression("Supplementary Figure S1: Rarefaction curve"), side = 3)
dev.off()

# import confounders metadata
confounder = read.csv('manuscript3_meta_confounders.csv', header = T)
confounder2 = confounder %>% filter(BioProject %in% sample_data(ps.rev)$BioProject) 

# trim low abundant taxa. 100 and 1000 both verified to still exhibit ethnicity significance after controlling for other factors
toremove = taxa_sums(ps.rev.genus) > 1000
ps.rev.genus.trim = prune_taxa(toremove, ps.rev.genus)

# rarefying also results in the significance of ethnicity.
ps.rev.genus.rf = ps.rev.genus %>% rarefy_even_depth(sample.size = 10000)
toremove = taxa_sums(ps.rev.genus.rf) > 1000
ps.rev.genus.trim = prune_taxa(toremove, ps.rev.genus.rf)

# ifelse to add 16sregion and kit to sample_data(ps.rf.clr)
testdf = sample_data(ps.rev.genus.trim) %>% data.frame
test_list = lapply(confounder2$BioProject, function(x) {
  testdf %>% filter(BioProject == x) %>% rownames_to_column() %>%
    mutate(preservative = confounder2[confounder2$BioProject == x,]$preservative) %>%
    mutate(vregion = confounder2[confounder2$BioProject == x,]$X16s) %>%
    mutate(kit = confounder2[confounder2$BioProject == x,]$kit)
})

testdf2 = test_list %>% bind_rows

testdf3 = testdf2 %>% column_to_rownames()
testdf4 = testdf3[order(match(rownames(testdf3), sample_names(ps.rev.genus.trim))),]
(rownames(testdf4) == sample_names(ps.rev.genus.trim)) %>% table

testdf4$vregion = testdf4$vregion %>% as.factor

# kit metadata
testdf4$kit = testdf4$kit %>% as.factor

#replace na in preservatives to none
testdf4$preservative = testdf4 %>% 
  pull(preservative) %>% 
  str_replace(pattern = '^na$', replacement = 'none') %>% 
  as.factor

testdf4$preservative_binary = ifelse(testdf4$preservative == 'none', 'none', 'yes') %>% as.factor


ps.rev.genus.trim.clr = jd_clr(ps.rev.genus.trim)
permanova_model = adonis(otu_table(ps.rev.genus.trim.clr) %>% t ~ preservative_binary + vregion + kit + ce, data = testdf4, permutations = 999, method = 'euclidean')
permanova_model

# individual PCA per included sample
sample_data(ps.rev.genus.trim.clr)$BioProject = sample_data(ps.rev.genus.trim.clr)$BioProject %>% as.factor 

# ordination function
jd_plot_ord2 = function (ps, var, method = "PCoA", distance = "euclidean") {
    sampledata = phyloseq::sample_data(ps)
    sampledata = subset(sampledata, !is.na(sampledata[[var]]))
    ps2 = phyloseq::phyloseq(sampledata, otu_table(ps), tax_table(ps))
    ord <- phyloseq::ordinate(ps2, method = method, distance = distance)
    fig <- phyloseq::plot_ordination(ps2, type = "samples", 
                                     ordination = ord, color = var) + ggplot2::theme_bw() + 
      ggsci::scale_color_igv(palette = "default")
    return(fig)
 
}

# figlist_pca = lapply(sample_data(ps.rev.genus.trim.clr)$BioProject %>% unique, function(var) {
#   df = sample_data(ps.rev.genus.trim.clr) %>% data.frame %>% filter(BioProject == var)
#   ps = prune_samples(rownames(df), ps.rev.genus.trim) %>% jd_clr()
#   jd_plot_ord2(ps, var = 'BioProject')
#   })
# library(ggpubr)
# ggarrange(plotlist = figlist_pca)

ps.rev.genus.trim.clr2 = merge_phyloseq(sample_data(testdf4), ps.rev.genus.trim.clr)

library(ggsci)

# plot each confounding variable, facet wrapped for BioProject
library(ggpubr)

confounding = c('country2', 'eth', 'kit', 'vregion', 'preservative')
figlist_pca = lapply(confounding, function(var) {
  jd_plot_ord2(ps.rev.genus.trim.clr2, var = var) + scale_color_d3() + facet_wrap(~BioProject)
})
ggarrange(plotlist = figlist_pca, nrow = 5, ncol = 1, labels = 'auto') %>% annotate_figure(top = 'Supplementary Figure S3: Ordination plot of possible confounders')
ggsave('fig_pca_confounder.pdf', height = 25, width = 10)

# complex heatmap of metadata

#tutorial data
nr1 = 4; nr2 = 8; nr3 = 6; nr = nr1 + nr2 + nr3
nc1 = 6; nc2 = 8; nc3 = 10; nc = nc1 + nc2 + nc3
mat = cbind(rbind(matrix(rnorm(nr1*nc1, mean = 1,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc1, mean = 0,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc1, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc2, mean = 0,   sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc2, mean = 1,   sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc2, mean = 0,   sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1*nc3, mean = 0.5, sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2*nc3, mean = 0.5, sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3*nc3, mean = 1,   sd = 0.5), nr = nr3))
)
mat = mat[sample(nr, nr), sample(nc, nc)] # random shuffle rows and columns
rownames(mat) = paste0("row", seq_len(nr))
colnames(mat) = paste0("column", seq_len(nc))

discrete_mat = matrix(sample(letters[1:4], 100, replace = TRUE), 10, 10)

Heatmap(discrete_mat)
Heatmap(testdf4$BioProject, name = 'test')

heatdf = confounder2 %>% select(c(BioProject, country, preservative, kit, X16s)) %>% 
  column_to_rownames(var = 'BioProject') %>% arrange(country, preservative, kit, X16s) %>%
  as.matrix()

hm_left = rowAnnotation(Preservative = heatdf %>% data.frame %>% pull(preservative),
                        Region = heatdf %>% data.frame %>% pull(X16s),
                        Kit = heatdf %>% data.frame %>% pull(kit))

pdf(file = 'fig_heatmap.pdf', width = 7.5, height = 7.5)
Heatmap(heatdf[,1], 
        row_split = heatdf[,1] %>% data.frame %>% pull(.),
        column_title = 'Supplementary Figure S2: Heatmap of metadata',
        row_title = NULL,
        row_gap = unit(5, 'mm'),
        column_title_gp = gpar(fontsize = 10),
        name = 'Country',
        border = F,
        heatmap_width = unit(52, 'mm'),
        left_annotation = hm_left,
        heatmap_legend_param = list(title = 'Country',
                                    labels = c('China', 'India', 'Indonesia', 'Malaysia')))

dev.off()  

# differential abundance analysis

da_eth = sample_data(ps.rev.genus.trim.clr2)$eth != 'Malay'
ps.ci = prune_samples(da_eth, ps.rev.genus.trim)

library(ALDEx2)

da_df = testdf4 %>% dplyr::select(c(country2, eth, preservative, vregion, kit)) %>% rownames_to_column() # %>% filter(da_eth)
da_df$rowname == sample_names(ps.rev.genus.trim)

da_asv = otu_table(ps.rev.genus.trim) %>% t()

mm = model.matrix(~ country2 + preservative + vregion + kit + eth, da_df)
da_clr = aldex.clr(da_asv, mm, mc.samples = 256, denom = 'all')
da_glm = ALDEx2::aldex.glm(da_clr)  
da_glm[da_glm[,74] < 0.05,] # col 74 refers to Indian


save.image()

##################
# sPLS-DA analysis
##################
library(mixOmics)

# control phyloseq object = 
ps.rev.genus.trim # for Chinese Indian Malay

ps.rev.genus.rf.ci = ps.rev.genus.rf %>% subset_samples(eth == 'Chinese' | eth == 'Indian') # for chinese and indian only
tokeep = (ps.rev.genus.rf.ci %>% taxa_sums()) > 1000
ps.rev.genus.trim.ci = prune_taxa(tokeep, ps.rev.genus.rf.ci)



# set training and validation set
trainset = ps.rev.genus.trim %>% subset_samples(country2 != 'Malaysia')
validset = ps.rev.genus.trim %>% subset_samples(country2 == 'Malaysia')

spl_asv = trainset %>% otu_table %>% data.frame %>% + 0.1 %>% logratio.transfo(logratio = 'CLR')
spl_meta = trainset %>% sample_data %>% data.frame %>% pull(eth)

spl_pca = pca(spl_asv, ncom = 10, center = T, scale = T)
plot(spl_pca)
plotIndiv(spl_pca, group = spl_meta, ind.names = F, legend = T)

# plsda
set.seed(4)
spl_plsda = plsda(spl_asv, spl_meta, ncomp = 10)
plotIndiv(spl_plsda, comp = 1:2,
          group = spl_meta, ind.names = F,
          ellipse = T, legend = T)
spl_perf_plsda = perf(spl_plsda, validation = 'Mfold', folds = 5, nrepeat = 50, auc = T, progressBar = T)
plot(spl_perf_plsda)
spl_perf_plsda$choice.ncomp
auroc(spl_plsda, roc.comp = 10)

# splsda
spl_keepx = seq(5,120, 5)
spl_splsda = tune.splsda(spl_asv, spl_meta, ncomp = 10, validation = 'Mfold', 
                         auc = T,
                         folds = 5, 
                         progressBar = T, 
                         dist = 'max.dist',
                         measure = 'overall',
                         test.keepX = spl_keepx,
                         nrepeat = 50,
                         cpus = 8)
spl_error = spl_splsda$error.rate
spl_ncomp = spl_splsda$choice.ncomp$ncomp
spl_select_keepx = spl_splsda$choice.keepX[1:spl_ncomp]
plot(spl_splsda)


spl_splsda_final = splsda(spl_asv, spl_meta, ncomp = spl_ncomp, keepX = spl_select_keepx)

fig.splsda.auc = auroc(spl_splsda_final)

fig.splsda.ord = plotIndiv(spl_splsda_final, comp = c(1,2),
          group = spl_meta,
          ind.names = T,
          ellipse = T,
          legend = T,
          title = 'sPLS-DA')

fig.splsda.compiled = 
ggarrange(
  fig.splsda.ord$graph + theme(legend.position = 'bottom'),
  fig.splsda.auc$graph.Comp1 + theme_bw() + theme(legend.position = 'bottom') + guides(color = guide_legend(ncol = 1)),
labels = 'auto')

ggsave('fig_splsda.pdf', height = 7.5, width = 7.5)


spl_perf_splsda = perf(spl_splsda_final, validation = 'Mfold', 
                       folds = 5, dist = 'max.dist', nrepeat = 10, progressBar = T)
spl_perf_splsda$error.rate

# spl_col = spl_meta %>% as.factor %>% as.numeric() %>% color.mixo()
# cim(spl_splsda_final,
#     row.sideColors = spl_col, 
#     cluster = 'col', 
#     transpose = T,
#     row.names = F,
#     col.names = F,
#     col.cex = 0.01,
#     legend = list(legend = spl_meta %>% as.factor() %>% levels()),
#     save = 'pdf',
#     name.save = 'heatmap')

# validation
spl_valid_asv = validset %>% otu_table %>% data.frame %>% + 0.1 %>% logratio.transfo(logratio = 'CLR')
spl_valid_meta = validset %>% sample_data %>% data.frame %>% pull(eth)


spl_valid_test = predict(spl_splsda_final, newdata = spl_valid_asv)
spl_output = spl_valid_test$class$max.dist[,spl_ncomp] 

# spl_ci = table(spl_output, spl_valid_meta) confusion matrix for chinese and indian
# spl_cim # confusion matrix for chinese indian malay
#1 - 41,24,23
#2 - 39, 24, 23
#3 - 41,24,23
#4 - 41,27,23

save.image()


# plotLoadings
spl_splsda_ci_ml  # splsda model for chinese indian (mainland)
spl_splsda_ci_my  # splsda model for chinese indian (malaysian)

spl_loading_my = plotLoadings(spl_splsda_ci_my, comp = 1, method = 'mean', contrib = 'max', ndisplay = 15) %>% 
  rownames_to_column() %>% 
  select(c(rowname, GroupContrib, importance)) 

library(ggsci)

fig.loading = bind_rows(spl_loading_ml, spl_loading_my, .id = 'model') %>% 
  mutate(model2 = ifelse(.$model == 1, 'Mainland', 'Malaysia')) %>% 
  ggplot(aes(x = rowname, y = importance, fill = GroupContrib)) + 
  geom_col() + 
  coord_flip() + 
  facet_wrap(~model2) + 
  scale_fill_d3(name= 'Ethnicity') +
  theme_bw() +
  xlab('') +
  ylab('Importance')

ggsave(file = 'fig_loading.png', height = 7.5, width = 7.5)

# reviewer 1 comment 9 and 16 PCA with ellipse
library(ggsci)
library(ggpubr)
fig.pca = ggarrange(
    ps.rev.genus.trim.clr2 %>% jd_plot_ord2(var = 'country2') + stat_ellipse(level = 0.60) + scale_color_d3(name = 'Country') + theme(legend.position = 'bottom') + guides(color = guide_legend(title.position = 'top', ncol = 1)),
    ps.rev.genus.trim.clr2 %>% jd_plot_ord2(var = 'eth') + stat_ellipse(level = 0.60) + scale_color_d3(name = 'Ethnicity') + theme(legend.position = 'bottom') + guides(color = guide_legend(title.position = 'top', ncol = 1)),
    ps.rev.genus.trim.clr2 %>% jd_plot_ord2(var = 'ce') + stat_ellipse(level = 0.60) + scale_color_d3(name = 'Ethnicity-Country') + theme(legend.position = 'bottom') + guides(color = guide_legend(title.position = 'top', ncol = 1)),
    labels = 'auto', nrow = 1)
ggsave('210527_fig_pca_ethnic.pdf', height = 7.5, width = 10)

save.image()

# reviewer 1 comment 11

jd_plot_alpha_stat2 = function (ps, var, measure = "Shannon", stat = "kruskal.test") {

    sampledata = phyloseq::sample_data(ps)
    sampledata = subset(sampledata, !is.na(sampledata[[var]]))
    ps = phyloseq::phyloseq(otu_table(ps), sampledata, tax_table(ps))
    
    rich = phyloseq::estimate_richness(ps)
    rich$Pielou = rich$Shannon/log(rich$Observed)
    
    sampledata = cbind(sampledata, rich)
    
    fig1 = ggplot2::ggplot(sampledata, aes_string(x = var, 
                                                  y = measure, fill = var)) + ggplot2::geom_boxplot() + 
      # ggpubr::stat_compare_means(method = stat) + 
      ggplot2::theme_bw() + 
      ggplot2::theme(legend.position = "blank") + 
      ggplot2::ylab(measure) + ggplot2::xlab(stringr::str_to_sentence(var)) + 
      ggplot2::geom_jitter(width = 0.1)
    
    return(fig1)

}

# stat calculation
jd_alpha_stat = function (ps, var, measure) {
  list = lapply(var, function(var) {
    sampledata = phyloseq::sample_data(ps)
    sampledata = subset(sampledata, !is.na(sampledata[[var]]))
    ps = phyloseq::phyloseq(otu_table(ps), sampledata, tax_table(ps))
    rich = phyloseq::estimate_richness(ps, measures = measure)
    data = sample_data(ps) %>% data.frame
    data[[measure]] = rich[[measure]]
    formula = as.formula(paste0(measure, " ~ ", var))
    statistic = ggpubr::compare_means(formula, data = data, p.adjust.method = 'BH')
    return(statistic)
  })
  return(list)
}

ps.rev.genus.rf %>% jd_alpha_stat(var = 'ce', measure = 'Shannon')

fig_alpha = ggarrange(
#chao1, eth
ps.rev.genus.rf %>% jd_plot_alpha_stat2(var = 'eth', measure = 'Chao1') + xlab('') + 
  geom_segment(aes(x = 1, xend = 2, y = 260, yend = 260)) +
  geom_segment(aes(x = 2, xend = 3, y = 270, yend = 270)) +
  geom_segment(aes(x = 1, xend = 3, y = 280, yend = 280)) +
  geom_text(x = 1.5, y = 265, label = '*') +
  geom_text(x = 2.5, y = 275, label = '*') +
  geom_text(x = 2, y = 285, label = '*') + stat_compare_means(method = 'kruskal.test', label.x = 1, label.y = 290),

# shannon, ethnicity
ps.rev.genus.rf %>% jd_plot_alpha_stat2(var = 'eth', measure = 'Shannon') + xlab('') + 
  # geom_segment(aes(x = 1, xend = 2, y = 4.5, yend = 4.5)) +
  geom_segment(aes(x = 2, xend = 3, y = 4.7, yend = 4.7)) +
  geom_segment(aes(x = 1, xend = 3, y = 4.9, yend = 4.9)) +
  # geom_text(x = 1.5, y = 4.6, label = 'ns') +
  geom_text(x = 2.5, y = 4.8, label = '*') +
  geom_text(x = 2, y = 5, label = '*') + stat_compare_means(method = 'kruskal.test', label.x = 1, label.y = 5.1),

# chao1, country
ps.rev.genus.rf %>% jd_plot_alpha_stat2(var = 'country2', measure = 'Chao1') + xlab('') + 
  geom_segment(aes(x = 1, xend = 2, y = 260, yend = 260)) +
  geom_segment(aes(x = 1, xend = 3, y = 270, yend = 270)) +
  geom_segment(aes(x = 1, xend = 4, y = 280, yend = 280)) +
  geom_segment(aes(x = 2, xend = 3, y = 290, yend = 290)) +
  geom_segment(aes(x = 3, xend = 4, y = 300, yend = 300)) +
  geom_text(x = 1.5, y = 265, label = '*') +
  geom_text(x = 2, y = 275, label = '*') +
  geom_text(x = 2.5, y = 285, label = '*') +
  geom_text(x = 2.5, y = 295, label = '*') +
  geom_text(x = 3.5, y = 305, label = '*') + stat_compare_means(method = 'kruskal.test', label.x = 1, label.y = 310),

# shannon, country
ps.rev.genus.rf %>% jd_plot_alpha_stat2(var = 'country2', measure = 'Shannon') + xlab('') + 
  geom_segment(aes(x = 2, xend = 4, y = 4.5, yend = 4.5)) +
  geom_segment(aes(x = 3, xend = 4, y = 4.7, yend = 4.7)) +
  geom_segment(aes(x = 1, xend = 4, y = 4.9, yend = 4.9)) +
  geom_text(x = 3, y = 4.6, label = '*') +
  geom_text(x = 3.5, y = 4.8, label = '*') +
  geom_text(x = 2.5, y = 5, label = '*') + stat_compare_means(method = 'kruskal.test', label.x = 1, label.y = 5.1),

# chao1, ce
ps.rev.genus.rf %>% jd_plot_alpha_stat2(var = 'ce', measure = 'Chao1') + xlab('') + 
  geom_segment(aes(x = 1, xend = 2, y = 260, yend = 260)) +
  geom_segment(aes(x = 1, xend = 4, y = 270, yend = 270)) +
  geom_segment(aes(x = 1, xend = 5, y = 280, yend = 280)) +
  geom_segment(aes(x = 2, xend = 5, y = 290, yend = 290)) +
  geom_segment(aes(x = 2, xend = 6, y = 300, yend = 300)) +
  geom_segment(aes(x = 3, xend = 1, y = 310, yend = 310)) +
  geom_segment(aes(x = 3, xend = 5, y = 320, yend = 320)) +
  geom_segment(aes(x = 4, xend = 5, y = 330, yend = 330)) +
  geom_segment(aes(x = 5, xend = 6, y = 340, yend = 340)) +
  geom_text(x = 1.5, y = 265, label = '*') +
  geom_text(x = 2.5, y = 275, label = '*') +
  geom_text(x = 3, y = 285, label = '*') +
  geom_text(x = 3.5, y = 295, label = '*') +
  geom_text(x = 4, y = 305, label = '*') +
  geom_text(x = 2, y = 315, label = '*') +
  geom_text(x = 4, y = 325, label = '*') +
  geom_text(x = 4.5, y = 335, label = '*') +
  geom_text(x = 5.5, y = 345, label = '*') + stat_compare_means(method = 'kruskal.test', label.x = 1, label.y = 350) + theme(axis.text.x = element_text(angle = 45, hjust = 1)),

# shannon, ce
ps.rev.genus.rf %>% jd_plot_alpha_stat2(var = 'ce', measure = 'Shannon') + xlab('') + 
  geom_segment(aes(x = 3, xend = 4, y = 4.5, yend = 4.5)) +
  geom_segment(aes(x = 2, xend = 3, y = 4.7, yend = 4.7)) +
  geom_segment(aes(x = 3, xend = 6, y = 4.9, yend = 4.9)) +
  geom_segment(aes(x = 4, xend = 5, y = 5.1, yend = 5.1)) +
  geom_segment(aes(x = 2, xend = 5, y = 5.3, yend = 5.3)) +
  geom_segment(aes(x = 5, xend = 6, y = 5.5, yend = 5.5)) +
  geom_segment(aes(x = 1, xend = 4, y = 5.7, yend = 5.7)) +
  geom_segment(aes(x = 1, xend = 2, y = 5.9, yend = 5.9)) +
  geom_segment(aes(x = 1, xend = 6, y = 6.1, yend = 6.1)) +
  geom_text(x = 3.5, y = 4.6, label = '*') +
  geom_text(x = 2.5, y = 4.8, label = '*') +
  geom_text(x = 4.5, y = 5, label = '*') +
  geom_text(x = 4.5, y = 5.2, label = '*') +
  geom_text(x = 3.5, y = 5.4, label = '*') +
  geom_text(x = 5.5, y = 5.6, label = '*') +
  geom_text(x = 2.5, y = 5.8, label = '*') +
  geom_text(x = 1.5, y = 6, label = '*') +
  geom_text(x = 3.5, y = 6.2, label = '*') + stat_compare_means(method = 'kruskal.test', label.x = 1, label.y = 6.3) + theme(axis.text.x = element_text(angle = 45, hjust = 1)),

ncol = 2, nrow = 3, labels = 'auto')

ggsave(file = 'fig_alpha.pdf', height = 13.5, width = 8)


save.image()

# reviewer 1 comment 15 heatmap

cor.otu = ps.rev.genus.trim.clr2 %>% otu_table %>% t %>% data.frame
cor.meta = ps.rev.genus.trim.clr2 %>% sample_data %>% data.frame

cor.bound = bind_cols(cor.otu, cor.meta['eth'], cor.meta['country2'], cor.meta['ce'])

#transform ethnicity to dummy numerics 1 Chinese 2 Indian 3 Malay
cor.mat = cor.bound %>% mutate(eth_Chinese = ifelse(.$eth == 'Chinese', 1, 0)) %>%
          mutate(eth_Indian = ifelse(.$eth == 'Indian', 1, 0)) %>%
          mutate(eth_Malay = ifelse(.$eth == 'Malay', 1, 0)) %>%
          mutate(country_CN = ifelse(.$country2 == 'China', 1, 0)) %>%
          mutate(country_IN = ifelse(.$country2 == 'India', 1, 0)) %>%
          mutate(country_ID = ifelse(.$country2 == 'Indonesia', 1, 0)) %>%
          mutate(country_MY = ifelse(.$country2 == 'Malaysia', 1, 0)) %>%
          mutate(ce_Chinese.CN = ifelse(.$ce == 'Chinese.China', 1, 0)) %>%
          mutate(ce_Chinese.MY = ifelse(.$ce == 'Chinese.Malaysia', 1, 0)) %>%
          mutate(ce_Indian.IN = ifelse(.$ce == 'Indian.India', 1, 0)) %>%
          mutate(ce_Indian.MY = ifelse(.$ce == 'Indian.Malaysia', 1, 0)) %>%
          mutate(ce_Malay.ID = ifelse(.$ce == 'Malay.Indonesia', 1, 0)) %>%
          mutate(ce_Malay.MY = ifelse(.$ce == 'Malay.Malaysia', 1, 0)) %>%
          select(!c(eth, country2, ce)) %>% cor(method = 'spearman') 

#get significance of each correlation
cor.mat.sig = cor.mtest(cor.mat, conf.level = 0.95)

pdf('fig_corrplot.pdf', height = 10, width = 10)
par(mfrow = c(2, 2))
corrplot(cor.mat[1:31,125:137],
         p.mat = cor.mat.sig$p[1:31,125:137],
         method = 'color',
         insig = 'label_sig',
         sig.level = .05, 
         pch.cex = .9, 
         pch.col = "black",
         cl.pos = 'r',
         cl.ratio = 0.5)

corrplot(cor.mat[32:62,125:137],
         p.mat = cor.mat.sig$p[32:62,125:137],
         method = 'color',
         insig = 'label_sig',
         sig.level = .05, 
         pch.cex = .9, 
         pch.col = "black",
         cl.pos = 'r',
         cl.ratio = 0.5)

corrplot(cor.mat[63:93,125:137],
         p.mat = cor.mat.sig$p[63:93,125:137],
         method = 'color',
         insig = 'label_sig',
         sig.level = .05, 
         pch.cex = .9, 
         pch.col = "black",
         cl.pos = 'r',
         cl.ratio = 0.5)

corrplot(cor.mat[94:124,125:137],
         p.mat = cor.mat.sig$p[94:124,125:137],
         method = 'color',
         insig = 'label_sig',
         sig.level = .05, 
         pch.cex = .9, 
         pch.col = "black",
         cl.pos = 'r',
         cl.ratio = 0.5)

dev.off()
par(mfrow = c(1,1))

# calculate partial correlation stats adjusted for potential confounders

library(ppcor)

cor.bound2 = bind_cols(cor.otu, cor.meta[,c('eth', 'country2', 'ce', 'kit', 'vregion', 'preservative')])

#transform ethnicity to dummy numerics 1 Chinese 2 Indian 3 Malay
cor.mat2 = cor.bound2 %>% mutate(eth_Chinese = ifelse(.$eth == 'Chinese', 1, 0)) %>%
  mutate(eth_Indian = ifelse(.$eth == 'Indian', 1, 0)) %>%
  mutate(eth_Malay = ifelse(.$eth == 'Malay', 1, 0)) %>%
  select(!c(eth, ce))

cor.mat2$country2 = cor.mat2$country2 %>% recode('China' = 1, 'India' = 2, 'Indonesia' = 3, 'Malaysia' = 4)
cor.mat2$kit = cor.mat2$kit %>% recode('qiagen_qiaampdnastoolminikit' = 1, 
                        'phenol_chloroform' = 2,
                        'qiagen_unspecified' = 3,
                        'roche_magnapurelctotalnucleicacidisolationkit' = 4,
                        'mobio_powersoildnaisolationkit' = 5,
                        'qiagen_qiaamppowerfecaldnakit' = 6,
                        'mpbio_fastdnaspinkit' = 7,
                        'mobio_powerfecaldnaisolationkit' = 8)  
cor.mat2$preservative = cor.mat2$preservative %>% recode('none' = 1, 'omnigene' = 2, 'rnalater' = 3)
cor.mat2$vregion = cor.mat2$vregion %>% recode('v34' = 1, 'v4' = 2)
  
pcor.control = cor.mat2[,c('country2', 'kit', 'preservative', 'vregion')]
pcor.eth = c('eth_Chinese', 'eth_Indian', 'eth_Malay')
pcor.asv = cor.mat2[,1:124] %>% names

pcor.df = lapply(pcor.eth, function(x) {
  lapply(pcor.asv, function(y) {
    pcor.test(x = cor.mat2[[x]], 
              y = cor.mat2[[y]],
              z = pcor.control,
              method = 'spearman')
  })
}) %>% bind_rows()

pcor.df.name = lapply(pcor.eth, function(x) {
  lapply(pcor.asv, function(y) {
    paste0(x, '_', y)
  })
}) %>% unlist()

pcor.df$comparison = pcor.df.name

write.csv(
pcor.df %>% filter(p.value <= 0.05) %>% arrange(desc(abs(estimate))),
file = 'pcor_ethnicity.csv', row.names = F, quote = F)

save.image()

# LDM
library(LDM)

ldm_asv = ps.rev.genus.trim %>% otu_table %>% t %>% data.frame 
ldm_meta = ps.rev.genus.trim.clr2 %>% sample_data %>% data.frame

ldm_fit = ldm(formula = ldm_asv | (country2 + preservative + vregion + kit) ~ eth, data = ldm_meta, n.perm.max = 0)

ldm_fit$VE.global.freq.submodels # variance explained 
ldm_fit$VE.otu.freq.submodels  # variance explained by each asv
ldm_fit$F.global.freq # F statistics

# residual plot for raw data
ldm_scree_freq = c(ldm_fit$VE.global.freq.submodels / ldm_fit$VE.df.submodels, ldm_fit$VE.global.freq.residuals)
ldm_color = c('red', rep('black', length(ldm_scree_freq) - 1))
plot(ldm_scree_freq / sum(ldm_scree_freq), main = 'Frequency Scale', xlab = 'Component', ylab = 'Proportion of total sum of squares', col = ldm_color)

# residual plot for transformed data
ldm_scree_tran = c(ldm_fit$VE.global.tran.submodels / ldm_fit$VE.df.submodels, ldm_fit$VE.global.tran.residuals)
ldm_color = c('red', rep('black', length(ldm_scree_freq) - 1))
plot(ldm_scree_tran / sum(ldm_scree_tran), main = 'Frequency Scale', xlab = 'Component', ylab = 'Proportion of total sum of squares', col = ldm_color)

# significance ldm test with permutations
ldm_res = ldm(formula = ldm_asv | (country2 + preservative + vregion + kit) ~ eth, data = ldm_meta, fdr.nominal = 0.05)

# explore some outputs
ldm_res$global.tests.stopped
ldm_res$otu.tests.stopped
ldm_res$p.global.omni # significance
ldm_res$detected.otu.omni # significant taxa associated with ethnicity
ldm_res$p.global.omni # global pval
ldm_res$q.otu.omni %>% t # FDR corrected

# summary table
ldm_df = data.frame(tax = ps.rev.genus.trim %>% tax_table() %>% data.frame %>% pull(Genus),
           pval = ldm_res$p.global.omni,
           qval = ldm_res$q.otu.omni %>% t) %>% arrange(qval)

# summary table 2
ldm_w1 = match(ldm_res$detected.otu.omni, colnames(ldm_res$q.otu.omni))
ldm_o = ldm_w1[order(ldm_res$p.otu.omni[1, ldm_w1])]

ldm_sum = data.frame(mean.freq = signif(ldm_res$mean.freq[ldm_o],3),
                     pval = signif(ldm_res$p.otu.omni[1,ldm_o],3),
                     adj.pval = signif(ldm_res$q.otu.omni[1,ldm_o],3),
                     directionIndian = t(ifelse(ldm_res$beta['ethIndian',]>0, '+', '-'))[ldm_o,],
                     directionMalay =  t(ifelse(ldm_res$beta['ethMalay',]>0, '+', '-'))[ldm_o,],
                     ASV = colnames(ldm_res$q.otu.omni)[ldm_o],
                     row.names = NULL)
# write.csv(ldm_sum, file = 'ldm_summary.csv', quote = F, row.names = F)

# boxplot of direction
ldm_adj.data1 = adjust.data.by.covariates(formula = ~ country2 + preservative + vregion + kit, 
                                          data = ldm_meta, 
                                          otu.table = ldm_asv, 
                                          center.otu.table = F)

# png(filename = 'fig_ldm_bp.png', units = 'in', height = 7.5, width = 10, res = 300)
# par(mfrow = c(2,4))
# for(i in 1:ldm_res$n.detected.otu.omni) {
#   boxplot(ldm_adj.data1$y.freq[,ldm_o[i]] ~ ldm_meta$eth,
#           main = ldm_res$detected.otu.omni[i],
#           ylab = 'Relative abundance',
#           xlab = 'Ethnicity')
#   
# }
# dev.off()
# par(mfrow = c(1,1))

# sadeq comment ligilactobacillus
bind_cols(
  ps.rev.genus.trim.clr2 %>% otu_table() %>% data.frame %>% t %>% data.frame %>% pull('ASV13_Bifidobacterium'),
  ps.rev.genus.trim.clr2 %>% sample_data %>% data.frame %>% pull(BioProject),
  ps.rev.genus.trim.clr2 %>% sample_data %>% data.frame %>% pull(eth)
) %>% ggplot(aes(x = ...2, y = ...1, color = ...3)) + geom_boxplot() + coord_flip() + theme_classic() +
  xlab('Study ID') + ylab('ASV13_Bifidobacterium Abundance (CLR-normalised data)') + scale_color_discrete(name = 'Ethnicity')


bind_cols(
ps.rev.genus.trim.clr2 %>% otu_table() %>% data.frame %>% t %>% data.frame %>% dplyr::select(c('ASV76_Ligilactobacillus', 
                                                                                               'ASV278_Olsenella', 
                                                                                               'ASV1096_Paludicola', 
                                                                                               'ASV522_Libanicoccus', 
                                                                                               'ASV13_Bifidobacterium', 
                                                                                               'ASV736_Solobacterium', 
                                                                                               'ASV714_Coprobacter', 
                                                                                               'ASV730_Enterorhabdus')),
ps.rev.genus.trim.clr2 %>% sample_data %>% data.frame %>% pull(eth)
) %>% pivot_longer(!...9) %>% ggplot(aes(x = ...9, y = value, color = ...9)) + 
             geom_boxplot() + 
             facet_wrap(~name) + 
             theme_bw() + 
             theme(legend.position = 'none') +
            xlab('Ethnicity') +
  ylab('CLR-normalised abundance') + ggtitle('Supplementary Figure S4: LDM significant taxa')
             
ggsave('fig_ldm_taxa.pdf', height = 7.5, width = 7.5)
save.image()

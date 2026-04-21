#####################################################################################################################
#### Code for paper: 
### Title:  Investigating host–parasite associations in cuckoo bumblebees using joint species distribution models  
#### AUTHOR LIST: Francesco Cerini, Elisa Serafini, Daniele Delle Monache, Andrea Chiocchio, Giuseppe Martino, Antonino Siclari, Daniele Canestrelli
#####################################################################################################################

rm(list = ls())

#packages
require(tidyverse)
require(data.table)
require(openxlsx)
require(Hmsc)
require(bipartite)
require(pals)
require(circlize)
require(phytools)
require(ggpubr)
require(phangorn)
require(patchwork)
library(circlize)
library(gridExtra)
library(grid)
library(reshape2) 
require(corrplot)
library(ape)
library(pheatmap)
library(gridExtra)
library(grid)
library(ggtree)
library(gridGraphics) 
library(usdm)
## load Bombus data
b_data<-read.xlsx("Data/Dataset_final_bombus2021_github.xlsx",
                  detectDates = T)

## 0. Data wrangling, cleaning and initial plots ####
#quick table to see how rare are the rare species
table(b_data$Genetic_identification) %>% as.data.frame()

## create dataframe to visualize species per site and abundance and change names to specify the Psythirus subgenus
b_data_sel<-b_data %>% group_by(Site, Genetic_identification) %>% 
  summarise(Abundance = n()) %>% 
  mutate(
    Genetic_identification = str_replace(
      Genetic_identification,
      "B\\. (rupestris|vestalis|barbutellus|sylvestris)",
      "B. P. \\1"
    )
  )

## plot to see how the species are distributed in the sites

# First define species color palette ONCE
species_levels <- sort(unique(b_data_sel$Genetic_identification))  # or a custom order
species_colors <- setNames(kelly(length(species_levels)), species_levels)

sp_sites_p<-ggplot(b_data_sel %>%
                     group_by(Site) %>%
                     mutate(TotalAbundance = sum(Abundance)) %>%
                     ungroup() %>% 
                     mutate(Site_new = as.numeric(sub("^S", "",Site))), 
                   aes(x = reorder(Site_new, -TotalAbundance),
                       y = log(Abundance),
                       fill = Genetic_identification))+
  geom_bar(stat = "identity",
           position = "stack",
           col = "black")+
  scale_fill_manual(values=species_colors)+
  # scale_fill_viridis_d(option = "C")+
  theme_bw(base_size = 16)+
  labs(y = "Abundance (Log scale)",
       x = "Site",
       fill = "Species")+
  theme(legend.text = element_text(face = "italic"))

ggsave(sp_sites_p,
       file = "Results/Figure_2_Species_in_sites_ab.png",
       width = 25, height = 20, units = "cm", dpi = 600)


## 1. create matrix for hmsc (Y) ####
b_data_mat<-b_data_sel %>% ungroup() %>% 
  dplyr::filter(!Genetic_identification %in% c("B. ruderatus","B. sylvarum")) %>% 
  pivot_wider(id_cols = Site,
              names_from = Genetic_identification,
              values_from = Abundance,
              values_fill = 0) %>% 
  ungroup() %>% column_to_rownames("Site") %>% as.matrix()


## 2. Phylogeny #####
## Load Phylogenetic tree
tree<-ape::read.tree("Data/Bombus_def.nwk")

# change tip names
old_names <- tree$tip.label
plot(tree)
new_names<-c("B. lapidarius", "B. pratorum", "B. lucorum", "B. terrestris", "B. hortorum", "B. pascuorum", 
             "B. P. rupestris", "B. P. vestalis", "B. P. sylvestris", "B. P. barbutellus")

# apply
tree$tip.label <- new_names
plot(tree)

# check
sort(tree$tip.label)
sort(colnames(b_data_mat))


## 2.1. Species traits ####
sp_tr<-read.xlsx("Data/Bombus_sp_traits2021.xlsx") %>%
  mutate(Social = as.factor(Social)) %>% 
  column_to_rownames("Species")

str(sp_tr)
all(c(rownames(sp_tr)) == c(colnames(b_data_mat)))
sp_tr <- sp_tr[colnames(b_data_mat), ]
all(c(rownames(sp_tr)) == c(colnames(b_data_mat)))


#check collinearity between tongue and inter tegular distance 
ggscatter(sp_tr ,
         x = "ITD_mm", 
         y = "Tongue_mm",
         add = "reg.line", size = 3, shape = 1)+
  stat_cor()+
  stat_regline_equation(label.y = 8)

#clearly tongue and ITC (a proxy for body size) are collinear, we opt to choose just one

## 3. Environment #### 
## Load site environmental data
env <- read.xlsx("Data/Sites_variables_final_clc_github.xlsx") %>%
  mutate(
    across(where(is.character), as.factor),   # convert all character columns to factors
    across(where(is.logical), as.factor)      # also convert logicals to factors if any
  )

#relate sites to species abundances and landcover classes
b_data_sel_lc<- merge(b_data_sel, env, by = "Site") %>% 
  select(Site, Genetic_identification, Abundance, LABEL3)
lc_sites<-b_data_sel_lc %>% select(Site, LABEL3) %>% 
  unique()
table(lc_sites$LABEL3)

### 3.1 Collinearity and redundancy test pipline ####
# extract numeric variables we intend to use in the model
env_num_mod<-env %>% 
  dplyr::select(c("Mean_t_ws_feb_apr",
                  "sd_t_ws_feb_apr",
                  "Altitude",
                  "height_mean",
                  "height_cvrob"))

cor_mat <- cor(env_num_mod, use = "pairwise.complete.obs")
vif_result <- vifstep(env_num_mod,
                      th = 10)  # threshold can be 5 or 10
vif_result
# we remove the seasonality (standard deviation of mean temperature)

## 4. HMSC ####
## Define the spatial Random factor, the sites, using the coordinates. First we create the coordinates object
xycoords<-env %>% dplyr::select(Site, Latitude, Longitude) %>% 
  column_to_rownames("Site") %>% as.matrix()

## then we define the spatial random level
rL1<-HmscRandomLevel(sData = xycoords)

## Now specify the study design, that is: which are the factors (the Sites) over which the environmental variables vary
studyDesign<-data.frame(Site = as.factor(env$Site))

#### 4.1. Model formula ####
## Model formula
XFormula_ms = ~ LABEL3 + Altitude + height_mean + height_cvrob + 
  poly(Mean_t_ws_feb_apr, degree = 2, raw = TRUE) 

# formula modelling influence of species traits on species responses to environmental variables, we used social strategy and the proxy for body size (ITD)
TrFormula = ~ Social + ITD_mm

# model with traits
model_ms_tr<-Hmsc(Y = b_data_mat, 
                  XData = env, 
                  phyloTree = tree,
                  TrData = sp_tr,
                  TrFormula = TrFormula,
                  XFormula = XFormula_ms, 
                  studyDesign = studyDesign, 
                  ranLevels = list(Site = rL1),
                  distr = "lognormal poisson")

## let's run the simulations to get estimates
## first let's set the model parameters
nChains = 4
thin = 500
samples = 1000
transient = 500*thin
verbose = 500*thin

## RUN THE MODEL with long chains for inference
ms_m_tr <- sampleMcmc(model_ms_tr,
                   samples = samples,
                   thin = thin,
                   transient = transient,
                   nChains = nChains,
                   verbose = verbose,
                   nParallel = 4
) #run time circa 8 hours 

# save model
saveRDS(ms_m_tr, file = "Results/hmsc_macroscale_model_traits.rds")

#load model
ms_m_tr<-readRDS("Results/hmsc_macroscale_model_traits.rds")
ms_m_tr$XFormula
ms_m_tr$XData
ms_m_tr$studyDesign
colnames(ms_m_tr$X)
head(ms_m_tr$X)

#### 4.2. Model diagnostics ####
mpost_ms_tr <- convertToCodaObject(ms_m_tr,
                                spNamesNumbers = c(T,F),
                                covNamesNumbers = c(T,F))

#plotting chains to visually check convergence and lack of autocorrelation in the chains
pdf("Results/macroscale_m_beta_plot_tr.pdf", width = 8, height = 6)
pdf("Results/macroscale_m_omega_plot_tr.pdf", width = 8, height = 6) 
par(mar = c(5, 5, 2, 2))
plot(mpost_ms_tr$Beta)
plot(mpost_ms_tr$Omega[[1]])
dev.off()

# Beta (species niches. responses to environmental variables),
psrf.beta_tr <-gelman.diag(mpost_ms_tr$Beta, multivariate = T)$psrf
# Omega (residual or environmentally constrained co-occurrence) parameters
psrf.omega_tr<-gelman.diag(mpost_ms_tr$Omega[[1]], multivariate = F)$psrf

# dataframe for plotting
psrf.df_ms_tr <- data.frame(
  value = c(psrf.beta_tr, psrf.omega_tr),
  parameter = rep(c("Beta", "Omega"),
                  c(length(psrf.beta_tr), length(psrf.omega_tr)))
)

# Create a plot object
psrf.p_tr <- ggplot(psrf.df_ms_tr, aes(x = value)) +
  geom_histogram(alpha = 0.6, position = "identity", 
                 col = "black",
                 bins = 20) +
  # scale_fill_manual(values = c("Beta" = "steelblue", "Omega" = "tomato")) +
  facet_wrap(~parameter)+
  theme_bw(base_size = 16)+
  labs(
    x = "PSRF value",
    y = "Count"
  )

#save plot
ggsave(psrf.p_tr,
       file = "Results/Figure_S1_psrf_model_tr.png",
       width = 20, height = 10, units = "cm", dpi = 600)

## check model explanatory power
predY_tr=computePredictedValues(ms_m_tr, 
                             expected = FALSE)
MF_tr = evaluateModelFit(hM = ms_m_tr, predY = predY_tr)
summary(MF_tr$SR2)
summary(MF_tr$RMSE)
summary(MF_tr$O.AUC)
summary(MF_tr$O.TjurR2)
summary(MF_tr$C.SR2)
summary(MF_tr$C.RMSE)

#### 4.4. Variance Partitioning ####
# let us group all the landcover codes column of the matrix into the "Land Cover" group to observe how much it contributes to explain the variance of the model

groupnames_ms = c("Land-cover", "Altitude", "Canopy","Temperature")
group_ms = c(
  rep(1, 9),  # Land cover
  2,           # altitude
  3, 3,        # canopy structure
  4, 4         # climate (polynomial terls)
)


#calculate VP
VP_ms_tr = computeVariancePartitioning(ms_m_tr, group = group_ms, 
                                    groupnames = groupnames_ms)

# Melt the data for ggplot
vp_df_ms_tr <- melt(VP_ms_tr$vals) 
colnames(vp_df_ms_tr) <- c("Factor", "Species", "Proportion")

# Custom ggplot of Variance Partitioning
vp_plot_ms_tr<-ggplot(vp_df_ms_tr,
                   aes(x = Species, y = Proportion, fill = Factor)) +
  geom_bar(stat = "identity", col = "black") +
  scale_fill_brewer(palette = "Pastel2") +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    legend.position = "right",  # try "bottom" if still overlaps
    legend.title = element_blank(),
  ) +
  labs(
    # title = "Variance Partitioning (Abundance model)",
    x = "Species",
    y = "Proportion of Variance Explained"
  )

# save plot
ggsave(vp_plot_ms_tr,
       file = "Results/Figure_3_vp_plot_tr.png",
       width = 20, height = 15, units = "cm", dpi = 600)



#### 4.5 Beta parameters (Species niches) ####
postBeta_ms_tr<-getPostEstimate(ms_m_tr, parName = "Beta")

#Get species and covariate names from model
spNames<- ms_m_tr$spNames
covNames<-ms_m_tr$covNames

#pull the data from postBeta for plotting
mbeta_ms_tr = postBeta_ms_tr$mean
sbetaP_ms_tr = postBeta_ms_tr$support

#Choose the support level for Beta plot
supportLevel <- 0.95

#format the data as sign and narrow to the covariates that have the desired support level
toPlot_ms_tr = sign(mbeta_ms_tr)
toPlot_ms_tr = toPlot_ms_tr * ((sbetaP_ms_tr > supportLevel) + 
                           (sbetaP_ms_tr < (1 - supportLevel)) > 0)

#format the data as matrix and add column and row names
betaMat_ms_tr = matrix(toPlot_ms_tr, nrow = ms_m_tr$nc, ncol = ncol(ms_m_tr$Y))
colnames(betaMat_ms_tr)<- spNames
rownames(betaMat_ms_tr) <- covNames

# create dataframe
betaMat_ms_tr<- as.data.frame(betaMat_ms_tr)
#update covariate names
# env vars
rownames(betaMat_ms_tr) <-c("Intercept",
                         "Broad-leaved forest",
                         "Complex cultivation patterns",
                         "Coniferous forest",
                         "Mixed forest",
                         "Non-irrigated arable land",
                         "Pastures",
                         "Sclerophyllous vegetation",
                         "Transitional woodland-shrub",
                         "Altitude",
                         "Mean Canopy Height",
                         "Canopy CV",
                         "Mean Temperature 1st",
                         "Mean Temperature 2nd")

#reformat for ggplot
betaMatmelt_ms_tr<-as.data.frame(melt(as.matrix(betaMat_ms_tr))) %>% 
  rename(Factor = Var1,
         Species = Var2)

#plot phylotree with ggtree
tree_p<-ggtree(tree, 
               branch.length = "none")

n_species <- length(tree$tip.label)

# Tree with fixed y limits
tree_p <- ggtree(tree,
                 branch.length = "none"
                 ) +
  scale_y_continuous(limits = c(0.9, n_species + 0.2))

#extract order of the plotted tips from the tree
tip_order<-get_taxa_name(tree_p)

## put the beta dataset with species following the display order of the phylo tree
betaMatmelt_ms_tr$Species <- factor(betaMatmelt_ms_tr$Species, levels = rev(tip_order))

beta_plot_ms_tr<-ggplot(betaMatmelt_ms_tr,
                     aes(x = Factor, y = Species, fill = factor(value))) +
  geom_tile(color = "black") +
  scale_fill_manual(
    values = c(
      "-1" = "dodgerblue3",
      "0"  = "white",
      "1"  = "firebrick3"
    ),
    name = "Sign",
    labels = c("-1" = "Negative", "0" = "Neutral", "1" = "Positive")
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white"),
    legend.key.width = unit(0.7, 'cm'),
    legend.key.height = unit(0.7, 'cm'),  # smaller legend box
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(face = "italic")
  ) +
  labs(
    x = NULL,
    y = "Species"
  )

# Combine
combined_ms<- tree_p + plot_spacer() + beta_plot_ms_tr + plot_layout(widths = c(5, -0.65 ,8), guides = "collect")

# save plot
ggsave(combined_ms,
       file = "Results/Figure_4_Phylotree_heatmap_tr.png",
       width = 23, height = 15, units = "cm", dpi = 300)


# # prediction plot to interpret betas
# Gradient = constructGradient(ms_m_tr, focalVariable = "Mean_t_ws_feb_apr")
# predY = predict(ms_m_tr, Gradient = Gradient, expected = TRUE)
# 
# plotGradient(ms_m_tr, Gradient, pred = predY, 
#              measure = "Y", index = 1,showData = TRUE)
# plotGradient(ms_m_tr, Gradient, pred = predY, 
#              measure = "Y", index = 2,showData = TRUE)
# plotGradient(ms_m_tr, Gradient, pred = predY, 
#              measure = "Y", index = 3,showData = TRUE)
# plotGradient(ms_m_tr, Gradient, pred = predY, 
#              measure = "Y", index = 4,showData = TRUE)
# plotGradient(ms_m_tr, Gradient, pred = predY, 
#              measure = "Y", index = 5,showData = TRUE)
# plotGradient(ms_m_tr, Gradient, pred = predY, 
#              measure = "Y", index = 6,showData = TRUE)
# plotGradient(ms_m_tr, Gradient, pred = predY, 
#              measure = "Y", index = 7,showData = TRUE)
# plotGradient(ms_m_tr, Gradient, pred = predY, 
#              measure = "Y", index = 8,showData = TRUE)
# plotGradient(ms_m_tr, Gradient, pred = predY, 
#              measure = "Y", index = 9,showData = TRUE)
# plotGradient(ms_m_tr, Gradient, pred = predY, 
#              measure = "Y", index = 10,showData = TRUE)

# 4.5.1. Phylogenetic signal ####
summary(mpost_ms_tr$Rho)[2]

#4.5.2. Species traits signal ####
postGamma = getPostEstimate(ms_m_tr, parName ="Gamma")

VP_ms_tr$R2T$Beta

VP_ms_tr$R2T$Y

#supp. plot of traits
#Get traits and covariate names from model
trNames<- ms_m_tr$trNames

#pull the data from postGamma for plotting
mgamma_ms_tr = postGamma$mean
sgamma_ms_tr = postGamma$support

#Choose the support level for Beta plot
supportLevel <- 0.95

#format the data as sign and narrow to the covariates that have the desired support level
toPlot_ms_tr_gamma = sign(mgamma_ms_tr)
toPlot_ms_tr_gamma = toPlot_ms_tr_gamma * ((sgamma_ms_tr > supportLevel) + 
                                 (sgamma_ms_tr < (1 - supportLevel)) > 0)

#format the data as matrix and add column and row names
gammaMat_ms_tr = matrix(toPlot_ms_tr_gamma, 
                        nrow = ms_m_tr$nc,
                        ncol = 4)

colnames(gammaMat_ms_tr)<- trNames
rownames(gammaMat_ms_tr) <- covNames

# create dataframe
gammaMat_ms_tr<- as.data.frame(gammaMat_ms_tr)
#update covariate names
# env vars
rownames(gammaMat_ms_tr) <-c("Intercept",
                            "Broad-leaved forest",
                            "Complex cultivation patterns",
                            "Coniferous forest",
                            "Mixed forest",
                            "Non-irrigated arable land",
                            "Pastures",
                            "Sclerophyllous vegetation",
                            "Transitional woodland-shrub",
                            "Altitude",
                            "Mean Canopy Height",
                            "Canopy CV",
                            "Mean Temperature 1st",
                            "Mean Temperature 2nd")

#reformat for ggplot
gammaMatmelt_ms_tr<-as.data.frame(melt(as.matrix(gammaMat_ms_tr))) %>% 
  rename(Factor = Var1,
         Species = Var2)

gamma_plot_ms_tr<-ggplot(gammaMatmelt_ms_tr,
                        aes(x = Factor, y = Species, fill = factor(value))) +
  geom_tile(color = "black") +
  scale_fill_manual(
    values = c(
      "-1" = "dodgerblue3",
      "0"  = "white",
      "1"  = "firebrick3"
    ),
    name = "Sign",
    labels = c("-1" = "Negative", "0" = "Neutral", "1" = "Positive")
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white"),
    legend.key.width = unit(0.7, 'cm'),
    legend.key.height = unit(0.7, 'cm'),  # smaller legend box
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.y = element_blank(),
  ) +
  labs(
    x = NULL,
    y = "Species"
  )


# save plot
ggsave(gamma_plot_ms_tr,
       file = "Results/Figure_S2_traits.png",
       width = 23, height = 15, units = "cm", dpi = 300)


#### 4.6. Omega parameters (species x species residual association) ####
OmegaCor_ms_tr = computeAssociations(ms_m_tr)

#redefine support level (95% posterior distribution)
supportLevel =0.95

# 1 = insignificant (<0.95), 0 = significant (>=0.95)
p_mat_ms_tr <- ifelse(OmegaCor_ms_tr[[1]]$support >= 0.95, 1, 0) 
# Set up plotting layout

png("Results/Figure_5_species_correlations_tr.png",
    width = 2000, height = 1500, res = 300)

corrplot(OmegaCor_ms_tr[[1]]$mean,
         tl.col = "black",
         method = "circle",
         type = "upper",
         order = "original",
         diag = FALSE,
         p.mat = p_mat_ms_tr,
         insig = "pch",
         pch = 8,
         pch.col = "red",
         pch.cex = .8,
         font = 3)

dev.off()

#### END ####

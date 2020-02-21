## White TE. Structural colours reflect individual quality: a meta-analysis

## Reset
  rm(list=ls())

## Libraries
  library(metafor)
  library(PRISMAstatement)
  library(ape)
  library(rotl)
  library(metaAidR)
  library(tidyverse)

## Functions
  # MLMR function with common parameters, for simplicity
  mlmr <- function(dat, variable){
          rma.mv(yi, 
                 vi,
                 mods = ~ variable - 1,
                 random = list(~ 1 | obs, ~1 | study_id, ~1 | species),
                 R = list(species = phylo),
                 method = "REML", 
                 data = dat)
  }
  
## Data
  coldat <- read.csv("../data/data.csv", stringsAsFactors = FALSE)
  tree <- read.tree("../data/tree.tre")

## Prep phylogenetic data
  # Replate some synonymous species names
  coldat[which(coldat$species == 'Tetrao_tetrix'), ]$species <- "Lyrurus_tetrix"
  coldat[which(coldat$species == 'Cyanistes_caeruleus'), ]$species <- "Cyanistes_caeruleus_caeruleus"
  
  # Compute correlation matrix
  phylo <- vcv(tree, corr = TRUE) 
  
  # Plot tree
  png('../figs/fig_phylogeny.png', width = 21, height = 21, units = 'cm', res = 300)
    plot(tree, cex = .8, label.offset = .1, no.margin = TRUE)
  dev.off()

## Estimate and subset effect sizes
  coldat <- escalc("ZCOR", ri = r, ni = n, data = coldat, append = TRUE)
  
## Search strategy and publication bias
  
## PRISMA statement
  prisma(found = 3684,
         found_other = 6,
         no_dupes = 3482, 
         screened = 3482, 
         screen_exclusions = 3430, 
         full_text = 52,
         full_text_exclusions = 11,
         qualitative = 41,
         quantitative = 41,
         width = 800, 
         height = 800)
  
# Funnel plot and Eggers test
  # Null model
  all_null <- rma.mv(yi, 
                     vi,
                     random = list(~ 1 | obs, ~1 | study_id, ~1 | species),
                     R = list(species = phylo),
                     method = "REML", 
                     data = coldat)
  summary(all_null)
  I2(all_null, coldat$vi)

  # Eggers
  residuals <- residuals(all_null)
  precision <- sqrt(1/all_null$vi)
  model <- rma(yi = residuals, sei = 1/precision) 
  regtest(model, model = "lm")
  
  # Funnel
  png('../figs/fig_funnel.png', width = 21, height = 21, units = 'cm', res = 300)
    par(mar = c(6, 5, 4, 2))
      funnel(model, 
             xlab = 'Effect size (Fisher\'s Z)',
             ylab = 'Standard error',
             xlim = c(-3, 3))
  dev.off()
  
## Summary info
  length(coldat$r)  # No. effects
  length(unique(coldat$study_id))  # No. studies
  length(unique(coldat$species))  # No. species
  # Proportion of effects per class
  coldat %>% 
    group_by(class) %>% 
    summarise(n = n()) %>% 
    mutate(freq = n / sum(n))
  unique(coldat$species)  # Species list
  
## Moderators 
  # Quality
  m_quality <- mlmr(coldat, coldat$quality_measure)
  summary(m_quality)
  I2(m_quality, coldat$vi)
  
  # Experimental vs observational
  m_exp <- mlmr(coldat, coldat$exp_obs)
  summary(m_exp)
  I2(m_exp, coldat$vi)
  
  # Colour variable
  m_var <- mlmr(coldat, coldat$col_var)
  summary(m_var)
  I2(m_var, coldat$vi)
  
  # Sex
  m_sex <- mlmr(coldat, coldat$sex)
  summary(m_sex)
  I2(m_sex, coldat$vi)

  # Iridescence
  m_irid <- mlmr(coldat, coldat$iridescent)  
  summary(m_irid)
  I2(m_irid, coldat$vi)
  
  # Control present
  m_control <- mlmr(coldat, coldat$control)
  summary(m_control)
  I2(m_control, coldat$vi)
  
## Plotting
  # Tidy & compile directional effects
  effects <- data.frame(effect = row.names(rbind(all_null$b,
                                                   m_quality$b,
                                                   m_var$b,
                                                   m_sex$b,
                                                   m_exp$b,
                                                   m_irid$b,
                                                   m_control$b)),
                        b = c(all_null$b,
                              m_quality$b,
                              m_var$b,
                              m_sex$b,
                              m_exp$b,
                              m_irid$b,
                              m_control$b),
                        lci = c(all_null$ci.lb,
                                m_quality$ci.lb,
                                m_var$ci.lb,
                                m_sex$ci.lb,
                                m_exp$ci.lb,
                                m_irid$ci.lb,
                                m_control$ci.lb),
                        uci = c(all_null$ci.ub,
                                m_quality$ci.ub,
                                m_var$ci.ub,
                                m_sex$ci.ub,
                                m_exp$ci.ub,
                                m_irid$ci.ub,
                                m_control$ci.ub))
  effects$name <- c('intercept model',
                    'age',
                    'body condition',
                    'immune function',
                    'parasite resistance',
                    'brightness',
                    'chroma',
                    'hue',
                    'composite',
                    'female',
                    'male',
                    'not distinghusied',
                    'experimental',
                    'observational',
                    'iridescent (vs not)',
                    'included (vs not)')
  effects$n <- c(nrow(coldat),
                 length(which(coldat$quality_measure == 'age')),
                 length(which(coldat$quality_measure == 'condition')),
                 length(which(coldat$quality_measure == 'immune')),
                 length(which(coldat$quality_measure == 'parasite')),
                 length(which(coldat$col_var == 'brightness')),
                 length(which(coldat$col_var == 'chroma')),
                 length(which(coldat$col_var == 'hue')),
                 length(which(coldat$col_var == 'pca')),
                 length(which(coldat$sex == 'female')),
                 length(which(coldat$sex == 'male')),
                 length(which(coldat$sex == 'zcombined')),
                 length(which(coldat$exp_obs == 'exp')),
                 length(which(coldat$exp_obs == 'obs')),
                 length(which(coldat$iridescent == 1)),
                 length(which(coldat$control == 1)))
  
  # Back-transform to r
  effects$r <- transf.ztor(effects$b)
  effects$r.lci <- transf.ztor(effects$lci)
  effects$r.uci <- transf.ztor(effects$uci)
  
  # Output
  write.csv(effects, '../output/summary_effects.csv', row.names = FALSE)
  
  # Forest plot
  png('../figs/fig_forest.png', width = 21, height = 31, units = 'cm', res = 300)
  par(lwd = 2, mar = c(4, 12, 2, 6), oma = c(0, 0, 0, 0))
  
  # Directional
    plot(1:nrow(effects) ~ rev(effects$r),
         yaxt = 'n',
         xlab = 'Pearsons\'s r',
         ylab = '',
         ylim = c(0.5, nrow(effects) + 0.5),
         xlim = c(-0.6, 0.6),
         pch = 18)
    arrows(x0 = rev(effects$r), x1 = rev(effects$r.uci), 
           y0 = 1:nrow(effects), y1 = 1:nrow(effects), angle = 90, length = 0.1)
    arrows(x0 = rev(effects$r), x1 = rev(effects$r.lci), 
           y0 = 1:nrow(effects), y1 = 1:nrow(effects), angle = 90, length = 0.1)
    abline(v = 0)
    abline(h = 15.6, lty = 2)
    axis(2, at = c(1:nrow(effects)), 
         labels = rev(effects$name), 
         las = 1, 
         tick = FALSE, 
         cex.axis = 0.9, 
         font = 1)
    # Labels
    text(-1.1, 15.4, 'Quality', cex = 0.9, font = 2, xpd = TRUE)
    text(-1.1, 11.4, 'Component', cex = 0.9, font = 2, xpd = TRUE)
    text(-1.1, 7.4, 'Sex', cex = 0.9, font = 2, xpd = TRUE)
    text(-1.1, 4.4, 'Study', cex = 0.9, font = 2, xpd = TRUE)
    text(-1.1, 2.4, 'Optics', cex = 0.9, font = 2, xpd = TRUE)
    text(-1.1, 1.4, 'Control', cex = 0.9, font = 2, xpd = TRUE)
    
    # Sample sizes
    text(0.7, 16.5, 'n', cex = 0.9, font = 2, xpd = TRUE)
    text(0.7, 16, '186', cex = 0.9, font = 2, xpd = TRUE)
    text(0.7, 15, '37', cex = 0.9, font = 2, xpd = TRUE)
    text(0.7, 14, '102', cex = 0.9, font = 2, xpd = TRUE)
    text(0.7, 13, '11', cex = 0.9, font = 2, xpd = TRUE)
    text(0.7, 12, '36', cex = 0.9, font = 2, xpd = TRUE)
    text(0.7, 11, '60', cex = 0.9, font = 2, xpd = TRUE)
    text(0.7, 10, '57', cex = 0.9, font = 2, xpd = TRUE)
    text(0.7, 9, '50', cex = 0.9, font = 2, xpd = TRUE)
    text(0.7, 8, '19', cex = 0.9, font = 2, xpd = TRUE)
    text(0.7, 7, '29', cex = 0.9, font = 2, xpd = TRUE)
    text(0.7, 6, '146', cex = 0.9, font = 2, xpd = TRUE)
    text(0.7, 5, '11', cex = 0.9, font = 2, xpd = TRUE)
    text(0.7, 4, '52', cex = 0.9, font = 2, xpd = TRUE)
    text(0.7, 3, '134', cex = 0.9, font = 2, xpd = TRUE)
    text(0.7, 2, '67', cex = 0.9, font = 2, xpd = TRUE)
    text(0.7, 1, '28', cex = 0.9, font = 2, xpd = TRUE)
    
  dev.off()

  
################################################################################
##                  IMPORT LIBRARIES AND SPILLOVER DATA                       ##
################################################################################

# preliminary

library(tidyverse)
cIn = "C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/spillovers/"

# chose computation: subcortical / cortical

cSubcortical = "sottocort/" # Subcortical
cCortical = "cort_aggr/" # Cortical
cRegion = cCortical

# set inputs

vcFilenames = str_subset(list.files(paste0(cIn, cRegion)), "pairspill")
cLags = gsub(".*l([0-9]+)_.*", "\\1", vcFilenames[1])
vID = sub(".*id([0-9]+).*", "\\1", vcFilenames)
if(cRegion == cSubcortical){
  vOutl = c("178243", "198653", "115825", "550439", "562345", "706040",
            "169747", "158035", "581450", "412528")
} else if (cRegion == cCortical){
  vOutl = c("115825", "198653", "178243", "706040", "412528", "581450", "550439",
            "463040", "172130", "145834")
}
vID = vID[!vID %in% vOutl]


# read spillovers

lPairSpill = list()
if(cRegion == cSubcortical){
  mNet = 
    read_csv(paste0(cIn, cRegion, "/netspill_subcort_l", cLags, ".csv"))
  for(i in 1 : length(vID)){
    cID = vID[i]
    lPairSpill[[cID]] = 
      read_csv(paste0(cIn, cRegion, "pairspill_subcort_l", 
                      cLags, "_id", cID, ".csv"), show_col_types = FALSE)
    rownames(lPairSpill[[cID]]) = colnames(lPairSpill[[cID]])
  }
} else if (cRegion == cCortical){
  mNet = read_csv(paste0(cIn, cRegion, "/netspill_cort_l", cLags, ".csv"))
  for(i in 1 : length(vID)){
    cID = vID[i]
    lPairSpill[[cID]] = 
      read_csv(paste0(cIn, cRegion, "pairspill_cort_l", 
                      cLags, "_id", cID, ".csv"), show_col_types = FALSE)
    rownames(lPairSpill[[cID]]) = colnames(lPairSpill[[cID]])
  }
}

################################################################################
##                             NET SPILLOVERS                                 ##
################################################################################

# plot net spillovers

cTitle = 
  ifelse(cRegion == cSubcortical, "Subcortical", "Cortical")

mNet |> 
  as.data.frame() |>
  pivot_longer(cols = - ID, names_to = "var", values_to = "val") |> 
  separate(var, into = c("reg", "hemi"), sep = "_") |>
  ggplot(aes(x = hemi, y = val, col = reg, fill = reg)) +
  geom_hline(yintercept = 0, color = "grey90", size = 0.1) +
  geom_boxplot(alpha = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  labs(x = "Hemisphere / Region", y = "Net Spillover", title = cTitle) +
  ylim(- 2.5, 4) 

# anova net spillovers

library(rstatix)

mNetTidy = 
  mNet |>  
  as.data.frame() |>  # Keep row index
  pivot_longer(cols = -ID, names_to = "var", values_to = "val") |> 
  separate(var, into = c("reg", "hemi"), sep = "_")

hyp = 2 # 1: outliers, 2: normality, 3: eq variances, 4: sphericity

if (hyp == 1) { # rejected
  mNetTidy |> 
    group_by(hemi, reg) |> 
    identify_outliers(val)
} else if (hyp == 2) { # rejected
  mNetTidy |> 
    group_by(hemi, reg) |> 
    shapiro_test(val)
  # ggpubr::ggqqplot("val", ggtheme = theme_bw()) +
  # facet_grid(hemi ~ reg)
} else if (hyp == 3) { # rejected
  mNetTidy |> 
    group_by(hemi) |> 
    levene_test(val ~ reg)
} else if (hyp == 4) { # rejected
  box_m(data = mNetTidy[, "val", drop = FALSE], group = mNetTidy$reg)
} else {
  message("No valid hypothesis selected")
}

mNetTidy |> 
  anova_test(
    dv = val, wid = ID, 
    between = reg, within = hemi) |>
  get_anova_table() |>
  # write_csv(paste0(
  #   "C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/draft/net_spillovers/manova_",
  #   ifelse(cRegion == cSubcortical, "sottocort", "cort"),
  #   "_", cLags, ".csv"))
  print()

dPInt = 
  mNetTidy |> 
  anova_test(
    dv = val, wid = ID, 
    between = reg, within = hemi) |>
  get_anova_table() |>
  as_tibble() |>
  slice(3) |>
  pull(p) 


# post-hoc anova net spillovers

group_by_if = function(.data, ...){
  if(dPInt < 0.05){
    group_by(.data, ...)
  } else{
    group_by(.data)
  }
}

mNetTidy |> # main effect
  group_by(hemi) |>
  anova_test(dv = val, wid = ID, between = reg) |>
  get_anova_table() |>
  adjust_pvalue(method = "bonferroni") |>
  # write_csv(paste0(
  #   "C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/draft/net_spillovers/posthoc_maineff_",
  #   ifelse(cRegion == cSubcortical, "sottocort", "cort"),
  #   "_", cLags, ".csv"))
  print(n = 50)

mNetTidy |> # pairwise diff
  group_by_if(hemi) |>
  pairwise_t_test(val ~ reg, p.adjust.method = "bonferroni") |>
  # write_csv(paste0(
  #   "C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/draft/net_spillovers/posthoc_pairdiff_",
  #   ifelse(cRegion == cSubcortical, "sottocort", "cort"),
  #   "_", cLags, ".csv"))
  print(n = 50)

mNetTidy |> # pairwise less
  group_by_if(hemi) |> # turn-off when insignificant interaction
  t_test(val ~ reg, p.adjust.method = "bonferroni", alternative = "less") |>
  # write_csv(paste0(
  #   "C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/draft/net_spillovers/posthoc_pairless_",
  #   ifelse(cRegion == cSubcortical, "sottocort", "cort"),
  #   "_", cLags, ".csv"))
  print(n = 50)

mNetTidy |> # pairwise greater
  group_by_if(hemi) |> # turn-off when insignificant interaction
  t_test(val ~ reg, p.adjust.method = "bonferroni", alternative = "greater") |>
  # write_csv(paste0(
  #   "C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/draft/net_spillovers/posthoc_pairgreater_",
  #   ifelse(cRegion == cSubcortical, "sottocort", "cort"),
  #   "_", cLags, ".csv"))
  print(n = 50)

################################################################################
##                    NET SPILLOVERS SENSITIVITY ANALYSIS                     ##
################################################################################

vcFilenames_all = list.files(paste0(cIn, cRegion))
vLags = 
  sub(".*_l([0-9]+)_id.*", "\\1",
      vcFilenames_all[grepl("_l[0-9]+_id", vcFilenames_all)]) |> 
  unique()

lNetSpill_cort <- lNetSpill_subcort <- list()
for(p in 1 : length(vLags)){
  cLags = vLags[p]
  lNetSpill_subcort[[p]] = read_csv(paste0(cOut, "sottocort/netspill_subcort_l", cLags, ".csv"))
  lNetSpill_cort[[p]] = read_csv(paste0(cOut, "cort_aggr/netspill_cort_l", cLags, ".csv"))
}

do.call("rbind",
          lapply(seq_along(lNetSpill_cort),
                 FUN = function(i) {
                   lNetSpill_cort[[i]] |>
                     mutate(lag = as.numeric(vLags[i])) |>
                     select(- ID) |>
                     apply(2, mean)})) |>
  as.data.frame() |>
  mutate(lag = as.factor(lag)) |>
  pivot_longer(!lag, names_to = 'network', values_to = 'net') |>
  separate(network, into = c("hemi", "network"), sep = "_") |>
  
  ggplot2::ggplot(ggplot2::aes(x = network, y = net, fill = net, group = lag)) +
  ggplot2::geom_col(position = 'dodge') +
  facet_grid(. ~ hemi, switch = "y", scales = "free") +

  ggplot2::scale_fill_gradient2(
    low = "red4", mid = "white", high = "blue4",
    midpoint = 0, # limits = c(-1, 1),
    name = "Value"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    strip.text.y.left = ggplot2::element_text(angle = 90, size = 14, vjust = 0.5),
    strip.placement = "outside",
    strip.background = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1, size = 12 * 1.1),
    axis.title.x = ggplot2::element_blank(),
    legend.position = "none",
    plot.title = ggplot2::element_text(size = 12)) +
  ggtitle("Sensitivity Analysis: Net spillovers on Cortical nets")

################################################################################
##                          PAIRWISE SPILLOVERS                               ##
################################################################################

# rstudioapi::navigateToFile("C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/Desktop/p2_spill.R")
# rstudioapi::navigateToFile("C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/programs/p_read_spill.R")

cTitle = 
  ifelse(cRegion == cSubcortical, "Subcortical", "Cortical")

# plot matrices

mAvgPair = lPairSpill[[1]] / length(lPairSpill)
for(i in 2 : length(lPairSpill)){
  mAvgPair = mAvgPair + lPairSpill[[i]] / length(lPairSpill) 
}

vOrd <- c("DMN R", "DMN L", "FPN R", "FPN L", "LIM R", "LIM L",
          "VAN R", "VAN L", "DAN R", "DAN L", "SMN R", "SMN L", 
          "VIS R", "VIS L")
  
mAvgPair |> 
  as.matrix() |>
  reshape2::melt() |>
  mutate(Var1 = gsub("H$", "", Var1),
         Var2 = gsub("H$", "", Var2)) |>
  mutate(Var1 = gsub("_", " ", Var1),
         Var2 = gsub("_", " ", Var2)) |>
  mutate(Var1 = factor(Var1, levels = rev(vOrd)),
         Var2 = factor(Var2, levels = rev(vOrd))) |>
  ggplot(aes(Var2, Var1, fill = value * 100)) + # x cols (v2), y rows (v1)
  geom_tile() +
  scale_fill_gradient2(low = "#00B0F6", mid = "white", high = "#F8766D",
                                midpoint = 0, name = "Net") + 
  theme_minimal() +
  labs(x = "To", y = "From", title = cTitle) +
  theme(
    axis.text.x = element_text(size = 12 * 1.1, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12 * 1.1, vjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12))

# single observations

f_get_nth <- 
  function(v, n, type = c("max", "min")) {
    type <- match.arg(type) 
    if (type == "max") {
      sorted_indices <- order(v, decreasing = TRUE)
    } else {
      sorted_indices <- order(v, decreasing = FALSE)
    }
    return(sorted_indices[n])
  }

vSim = 
  lapply(X = lPairSpill,
         FUN = function(a, b) cor(c(as.matrix(a)), c(as.matrix(b))),
         b = mAvgPair) |>
  unlist()

for(i in seq(1, 151, 30)){
  id = names(lPairSpill)[f_get_nth(v = vSim, n = i)]
  p =
    lPairSpill[[id]] |>
    as.matrix() |>
    reshape2::melt() |>
    mutate(Var1 = gsub("H$", "", Var1),
           Var2 = gsub("H$", "", Var2)) |>
    mutate(Var1 = gsub("_", " ", Var1),
           Var2 = gsub("_", " ", Var2)) |>
    mutate(Var1 = factor(Var1, levels = rev(vOrd)),
           Var2 = factor(Var2, levels = rev(vOrd))) |>
    ggplot(aes(Var2, Var1, fill = value * 100)) + # x cols (v2), y rows (v1)
    geom_tile() +
    scale_fill_gradient2(low = "#00B0F6", mid = "white", high = "#F8766D",
                         midpoint = 0, name = "Net") + 
    theme_minimal() +
    labs(x = "To", y = "From", 
         title = paste0(cTitle, " (", i, "^: ", id, ")")) +
    theme(
      axis.text.x = element_text(size = 12 * 1.1, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 12 * 1.1, vjust = 0.5),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12))
  print(p)
}


# anova

f_extract_pair = 
  function(x){
    vOrd = c("DMN R", "DMN L", "FPN R", "FPN L", "LIM R", "LIM L", 
             "VAN R", "VAN L", "DAN R", "DAN L", "SMN R", "SMN L", "VIS R", "VIS L")
    x = as.matrix(x) 
    colnames(x) = gsub("_", " ", gsub("H", "", colnames(x)))
    rownames(x) = gsub("_", " ", gsub("H", "", rownames(x)))
    if(cRegion == cCortical){
      x = x[vOrd, vOrd]
    }
    vOut = x[upper.tri(x, diag = FALSE)]
    mInd = which(upper.tri(x, diag = FALSE), arr.ind = TRUE)
    names(vOut) = paste(rownames(x)[mInd[, 1]], 
                        colnames(x)[mInd[, 2]],
                        sep = "/")
    return(vOut)
    
  }

mPairVec = t(sapply(X = lPairSpill, FUN = f_extract_pair))

mPairVec |>
  as_tibble() |>
  rownames_to_column(var = "ID") |>
  pivot_longer(!ID, names_to = "int", values_to = "sp") %>%
  aov(sp ~ int, data = .) |>
  broom::tidy() |>
  write_csv(paste0(
    "C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/draft/pairwise_spillovers/anova_",
    ifelse(cRegion == cSubcortical, "sottocort", "cort"),
    "_", cLags, ".csv"))
  # residuals() |> ks.test("pnorm")  # rejected
  # rstatix::levene_test(formula = sp ~ int) # rejected


mPairVec |> # post-hoc diff
  as_tibble() |>
  rownames_to_column(var = "ID") |>
  pivot_longer(!ID, names_to = "int", values_to = "sp") |>
  rstatix::pairwise_t_test(sp ~ int, p.adjust.method = "bonferroni") |>
  write_csv(paste0(
    "C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/draft/pairwise_spillovers/posthoc_pairdiff_",
    ifelse(cRegion == cSubcortical, "sottocort", "cort"),
    "_", cLags, ".csv"))
  # View()

mPairVec |> # post-hoc less
  as_tibble() |>
  rownames_to_column(var = "ID") |>
  pivot_longer(!ID, names_to = "int", values_to = "sp") |> 
  rstatix::t_test(sp ~ int, alternative = "less", p.adjust.method = "bonferroni") |>
  write_csv(paste0(
    "C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/draft/pairwise_spillovers/posthoc_pairless_",
    ifelse(cRegion == cSubcortical, "sottocort", "cort"),
    "_", cLags, ".csv"))
  # print(n = 50)

mPairVec |> # post-hoc greater
  as_tibble() |>
  rownames_to_column(var = "ID") |>
  pivot_longer(!ID, names_to = "int", values_to = "sp") |> 
  rstatix::t_test(sp ~ int, alternative = "greater", p.adjust.method = "bonferroni") |>
  write_csv(paste0(
    "C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/draft/pairwise_spillovers/posthoc_pairgreater_",
    ifelse(cRegion == cSubcortical, "sottocort", "cort"),
    "_", cLags, ".csv"))
  # print(n = 50)

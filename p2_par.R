library(tidyverse)
library(ggpubr)
library(latex2exp)

##############################################################################
##                                IMPORT                                    ##
##############################################################################

cIn = "C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/estimates/"
cDirecS = "C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/brain_data/CSV/"
cDirecC = "C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/brain_data/HCP_Net_BOLD/"
vcFilenamesS = list.files(cDirecS)
vcFilenamesC = list.files(cDirecC)
vID = gsub("\\D", "", vcFilenamesS)

cIn = "C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/estimates/"
for(j in 1 : length(vID)){
  cID = vID[j]
  if(j == 1){
    tEstS = 
      read.csv(paste0(cIn, "sottocort/est_sottocort_l51015_", cID, ".csv")) 
    tEstC = 
      read.csv(paste0(cIn, "cort_aggr/est_cortaggr_l51015_", cID, ".csv")) 
  } else {
    tEstS = 
      tEstS |>
      left_join(
        read.csv(paste0(cIn, "sottocort/est_sottocort_l51015_", cID, ".csv")),
        by = "par"
      )
    tEstC = 
      tEstC |>
      left_join(
        read.csv(paste0(cIn, "cort_aggr/est_cortaggr_l51015_", cID, ".csv")),
        by = "par"
      )
  }
}
colnames(tEstC) <- colnames(tEstS) <- c("par", vID)

vOutS = c("178243", "198653", "115825", "550439", "562345", "706040",
          "169747", "158035", "581450", "412528")

vOutC = c("115825", "198653", "178243", "706040", "412528", "581450", "550439",
          "463040", "172130", "145834")

##############################################################################
##                              COMPARISON H                                ##
##############################################################################

tCort = 
  tEstC |>
  pivot_longer(!par, names_to = "id", values_to = "val") |>
  mutate(par = gsub(pattern = "_mean", replacement = "", x = par)) |>
  filter(!id %in% vOutC)

tSubcort =
  tEstS |>
  pivot_longer(!par, names_to = "id", values_to = "val") |>
  filter(!id %in% vOutS)

rm(tEstC, tEstS)
vNetCort = c("DMN", "FPN", "LIM", "VAN", "DAN", "SMN", "VIS")
vNetSubcort = c("HIP", "AMY", "pTHA", "NAc", "GP", "aTHA", "PUT", "CAU")
tCofCort = read_csv("C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/draft/ts_analysis/cof_cortical.csv")
tCofSubcort = read_csv("C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/draft/ts_analysis/cof_subcortical.csv") |>
  mutate(id = as.integer(id))

H_subcort_both = 
  tSubcort |>
  filter(str_detect(par, "H_")) |>
  mutate(par = str_remove(par, "H_"), id = as.integer(id)) |>
  rename(net = par, mde = val) |>
  left_join(tCofSubcort, by = c("net", "id")) |>
  rename(cof = H) |>
  select(net, id, mde, cof) |>
  mutate(id = as.integer(id))

H_cort_both = 
  tCort |>
  filter(str_detect(par, "^H_")) |>
  separate(par, into = c("a", "side", "net"), sep = "_") |>
  unite("net", a, net, side, sep = "_") |>
  rename(mde = val) |>
  mutate(id = as.double(id),
         net = str_remove(net, "^H_")) |>
  left_join(tCofCort, by = c("net", "id")) |>
  rename(cof = H) |>
  select(net, id, mde, cof) |>
  mutate(id = as.integer(id)) 

# correlation
H_subcort_both |>
  group_by(net) |>
  summarise(cor = cor(cof, mde))

H_cort_both |>
  group_by(net) |>
  summarise(cor = cor(cof, mde))

# sign concordance
H_subcort_both |>
  mutate(mde = mde - 1/2,
         cof = cof - 1/2) |>
  mutate(conc = sign(mde) * sign(cof)) |>
  group_by(net) |>
  summarise(sum(conc) / n())
H_cort_both |>
  mutate(mde = mde - 1/2,
         cof = cof - 1/2) |>
  mutate(conc = sign(mde) * sign(cof)) |>
  group_by(net) |>
  summarise(sum(conc) / n())

# cosine similarity
H_subcort_both |>
  summarise((mde %*% cof) / ((mde %*% mde) * (cof %*% cof)))
H_cort_both |>
  summarise((mde %*% cof) / ((mde %*% mde) * (cof %*% cof)))

# distance
H_subcort_both |>
  summarise(mse = mean((mde - cof) ^ 2),
            eucl = sqrt(sum((mde - cof) ^ 2)))
H_cort_both |>
  summarise(mse = mean((mde - cof) ^ 2),
            eucl = sqrt(sum((mde - cof) ^ 2)))

# variability
summary(H_subcort_both[, c("mde", "cof")])
summary(H_cort_both[, c("mde", "cof")])

##############################################################################
##                                 UNIV                                     ##
##############################################################################

# Cortical
library(dplyr)
library(ggplot2)
library(tidyr)

tCort |>
  filter(!id %in% vOutC) |>
  filter(grepl("^nu_|^a_|^H_", par)) |>
  separate(par, into = c("par", "net", "hemi"), sep = "_") |>
  mutate(net_hemi = interaction(net, hemi, sep = "_"),
         par = recode(par,
                      "nu" = "nu",       
                      "a" = "alpha",
                      "H" = "H")) |>
  group_by(par, net, hemi) |>
  filter(val > quantile(val, 0.25) - 1.5 * IQR(val),
         val < quantile(val, 0.75) + 1.5 * IQR(val)) |>
  mutate(net = factor(net, levels = vNetSubcort)) |>
  select(par, val, net_hemi) |>
  group_by(par, net_hemi) |>
  summarise(avg = round(mean(val), 2)) |> 
  pivot_wider(id_cols = net_hemi, values_from = avg, names_from = par)

# tests of hypothesis for markovianity
sig_lvl = 0.05
for(k in 1 : 7){
  tCort |>
    filter(str_detect(par, vNetCort[k]) & str_detect(par, "^H_")) |>
    separate(par, into = c("par", "side", "LIM"), sep = "_") |>
    select(id, side, val) |>
    mutate(side = as.factor(side)) |>
    pivot_wider(id_cols = id, values_from = val, names_from = side) |>
    pivot_longer(cols = c(LH, RH), names_to = "side", values_to = "val") |>
    ggplot(aes(x = val, fill = side, color = side)) +
    geom_density(alpha = 0.4) +
    theme_minimal()
  vEst = 
    tCort |>
    filter(str_detect(par, vNetCort[k]) & str_detect(par, "^H_")) |>
    separate(par, into = c("par", "side", vNetCort[k]), sep = "_") |>
    select(id, side, val) |>
    mutate(side = as.factor(side)) |>
    pivot_wider(id_cols = id, values_from = val, names_from = side) |>
    summarise(sd_LH = sd(LH), sd_RH = sd(RH),
              mean_LH = mean(LH), mean_RH = mean(RH))
  print(abs((vEst[3] - 0.5) / vEst[1]) < qnorm(p = 1 - sig_lvl))
  print(abs((vEst[4] - 0.5) / vEst[2]) < qnorm(p = 1 - sig_lvl))
}

tCort |>
  filter(!id %in% vOutC) |>
  filter(grepl("^nu_|^a_|^H_", par)) |>
  separate(par, into = c("par", "hemi", "net"), sep = "_") |>
  mutate(net_hemi = interaction(net, hemi, sep = "_"),
         par = recode(par,
                      "nu" = "nu",       
                      "a" = "alpha",
                      "H" = "H")) |>
  group_by(par, net, hemi) |>
  filter(val > quantile(val, 0.25) - 1.5 * IQR(val),
         val < quantile(val, 0.75) + 1.5 * IQR(val)) |>
  mutate(net = factor(net, levels = vNetCort)) |>
  ggplot(aes(x = net_hemi, y = val, fill = net)) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_discrete(labels = function(x) gsub(".*_", "", x)) +
  facet_grid(par ~ net, scales = "free", space = "fixed", 
             labeller = label_parsed) +
  geom_hline(data = . %>% filter(par == "H") %>% distinct(par, net), 
             aes(yintercept = 0.5), linetype = "dashed", color = "black",
             linewidth = 0.7) +
  theme_linedraw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none")



# Subcortical

tSubcort |>
  filter(!id %in% vOutS) |>
  filter(grepl("^nu_|^a_|^H_", par)) |>
  separate(par, into = c("par", "net", "hemi"), sep = "_") |>
  mutate(net_hemi = interaction(net, hemi, sep = "_"),
         par = recode(par,
                      "nu" = "nu",       
                      "a" = "alpha",
                      "H" = "H")) |>
  group_by(par, net, hemi) |>
  filter(val > quantile(val, 0.25) - 1.5 * IQR(val),
         val < quantile(val, 0.75) + 1.5 * IQR(val)) |>
  mutate(net = factor(net, levels = vNetSubcort)) |>
  select(par, val, net_hemi) |>
  group_by(par, net_hemi) |>
  summarise(avg = round(mean(val), 2)) |> 
  pivot_wider(id_cols = net_hemi, values_from = avg, names_from = par)

tSubcort |>
  filter(!id %in% vOutS) |>
  filter(grepl("^nu_|^a_|^H_", par)) |>
  separate(par, into = c("par", "net", "hemi"), sep = "_") |>
  mutate(net_hemi = interaction(net, hemi, sep = "_"),
         par = recode(par,
                      "nu" = "nu",       
                      "a" = "alpha",
                      "H" = "H")) |>
  group_by(par, net, hemi) |>
  filter(val > quantile(val, 0.25) - 1.5 * IQR(val),
         val < quantile(val, 0.75) + 1.5 * IQR(val)) |>
  mutate(net = factor(net, levels = vNetSubcort)) |>
  ggplot(aes(x = net_hemi, y = val, fill = net)) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_discrete(labels = function(x) gsub(".*_", "", x)) +
  facet_grid(par ~ net, scales = "free", space = "fixed", 
             labeller = label_parsed) +
  geom_hline(data = . %>% filter(par == "H") %>% distinct(par, net), 
             aes(yintercept = 0.5), linetype = "dashed", color = "black",
             linewidth = 0.7) +
  theme_linedraw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none")

tSubcort |>
  filter(str_detect(par, "^H")) |>
  mutate(net = str_remove(par, "^H_"),
         id = as.double(id)) |>
  rename(H = val) |>
  dplyr::select(id, net, H) |>
  filter(!id %in% vOutS) |>
  separate(net, into = c("net", "hemi"), sep = "_") |>
  group_by(net) |>
  ggplot(aes(x = H, fill = hemi)) +
  # geom_density(alpha = 0.7, adjust = 1.5, color = "black") +
  geom_histogram(binwidth = 0.03, alpha = 0.7, position = "identity", color = "black") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 0.5) +
  facet_wrap(~ net) +  # Creates a separate plot for each 'net'
  labs(x = latex2exp::TeX("H"),
       y = "Density",
       fill = "Hemisphere") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"))

tCort |>
  filter(str_detect(par, "^H")) |>
  mutate(net = str_remove(par, "^H_"),
         id = as.double(id)) |>
  rename(H = val) |>
  dplyr::select(id, net, H) |>
  filter(!id %in% vOutC) |>
  separate(net, into = c("hemi", "net"), sep = "_") |>
  group_by(net) |>
  ggplot(aes(x = H, fill = hemi)) +
  # geom_density(alpha = 0.7, adjust = 1.5, color = "black") +
  geom_histogram(binwidth = 0.03, alpha = 0.7, position = "identity", color = "black") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 0.5) +
  facet_wrap(~ net) +  # Creates a separate plot for each 'net'
  labs(x = latex2exp::TeX("H"),
       y = "Density",
       fill = "Hemisphere") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"))

##############################################################################
##                                MULTIV                                    ##
##############################################################################

# Cortical

tCortSplit = 
  tCort |>
  group_by(par) |>
  summarise(val = mean(val)) |>
  filter(grepl("^rho|^eta", par)) |>
  separate(par, into = c("par", "V"), sep = "_", extra = "merge") |>
  separate(V, into = c("V1", "V2"), sep = "/") |>
  separate(V1, into = c("hemi_i", "net_i"), sep = "_") |>
  separate(V2, into = c("hemi_j", "net_j"), sep = "_") |>
  pivot_wider(names_from = "par", values_from = "val")
  
vNamesC = 
  expand.grid(hemi = c("LH", "RH"),
              net = vNetCort) |>
  with(paste(net, hemi, sep = "_"))

mRhoC =
  matrix(NA, nrow = length(vNamesC), ncol = length(vNamesC))
rownames(mRhoC) <- colnames(mRhoC) <- vNamesC
for(i in 1 : nrow(mRhoC)){
  for(j in 1 : ncol(mRhoC)){
    if (i == j) {
      mRhoC[i, j] =  1
    } else if (i != j) {
      name_i = colnames(mRhoC)[i]
      name_j = colnames(mRhoC)[j]
      vci = strsplit(name_i, "_")[[1]]
      vcj = strsplit(name_j, "_")[[1]]
      r1 = 
        tCortSplit |>
        filter(hemi_i == vci[2], hemi_j == vcj[2],
               net_i == vci[1], net_j == vcj[1]) |>
        pull(rho)
      r2 = 
        tCortSplit |>
        filter(hemi_i == vcj[2], hemi_j == vci[2],
               net_i == vcj[1], net_j == vci[1]) |>
        pull(rho)
      mRhoC[i, j] = ifelse(length(r2) == 0, r1, r2)
    }
  }
}

mRhoC |>
  as.data.frame() |>
  rownames_to_column("row") |>
  pivot_longer(-row, names_to = "col", values_to = "value") |>
  mutate(
    row_net = sub("_.*", "", row),
    row_hemi = sub(".*_", "", row),
    col_net = sub("_.*", "", col),
    col_hemi = sub(".*_", "", col),
    row_label = paste(row_net, row_hemi, sep = "\n"),
    col_label = paste(col_net, col_hemi, sep = "\n"),
    row_label = factor(row_label, levels = 
                         as.vector(t(outer(vNetCort, c("LH", "RH"),
                                           paste, sep = "\n")))),
    col_label = factor(col_label, levels = 
                         as.vector(t(outer(vNetCort, c("LH", "RH"), paste, sep = "\n")))),
    row_net = factor(row_net, levels = vNetCort),
    col_net = factor(col_net, levels = vNetCort)
  ) |>
  ggplot(aes(x = col_label, y = row_label, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(- 1, 1)) +
  facet_grid(row_net ~ col_net, space = "free", scales = "free", drop = FALSE) +
  labs(x = NULL, y = NULL, fill = TeX("$\\rho$")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 10),
    axis.text.y = element_text(angle = 0, hjust = 1, size = 10),
    strip.text.x = element_text(size = 14),
    strip.text.y = element_text(size = 14),
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines")
  ) +
  scale_x_discrete(labels = unique(df$col_hemi)) +
  scale_y_discrete(labels = unique(df$row_hemi))

# Subcortical

tSubcortSplit = 
  tSubcort |>
  group_by(par) |>
  summarise(val = mean(val)) |>
  filter(grepl("^rho|^eta", par)) |>
  separate(par, into = c("par", "V"), sep = "_", extra = "merge") |>
  separate(V, into = c("V1", "V2"), sep = "/") |>
  separate(V1, into = c("net_i", "hemi_i"), sep = "_") |>
  separate(V2, into = c("net_j", "hemi_j"), sep = "_") |>
  pivot_wider(names_from = "par", values_from = "val")

vNamesS = 
  expand.grid(hemi = c("LH", "RH"),
              net = vNetSubcort) |>
  with(paste(net, hemi, sep = "_"))

mRhoS =
  matrix(NA, nrow = length(vNamesS), ncol = length(vNamesS))
rownames(mRhoS) <- colnames(mRhoS) <- vNamesS
for(i in 1 : nrow(mRhoS)){
  for(j in 1 : ncol(mRhoS)){
    if (i == j) {
      mRhoS[i, j] =  1
    } else if (i != j) {
      name_i = colnames(mRhoS)[i]
      name_j = colnames(mRhoS)[j]
      vci = strsplit(name_i, "_")[[1]]
      vcj = strsplit(name_j, "_")[[1]]
      r1 = 
        tSubcortSplit |>
        filter(hemi_i == vci[2], hemi_j == vcj[2],
               net_i == vci[1], net_j == vcj[1]) |>
        pull(rho)
      r2 = 
        tSubcortSplit |>
        filter(hemi_i == vcj[2], hemi_j == vci[2],
               net_i == vcj[1], net_j == vci[1]) |>
        pull(rho)
      mRhoS[i, j] = ifelse(length(r2) == 0, r1, r2)
    }
  }
}

mRhoS |>
  as.data.frame() |>
  rownames_to_column("row") |>
  pivot_longer(-row, names_to = "col", values_to = "value") |>
  mutate(
    row_net = sub("_.*", "", row),
    row_hemi = sub(".*_", "", row),
    col_net = sub("_.*", "", col),
    col_hemi = sub(".*_", "", col),
    row_label = paste(row_net, row_hemi, sep = "\n"),
    col_label = paste(col_net, col_hemi, sep = "\n"),
    row_label = factor(row_label, levels = 
                         as.vector(t(outer(vNetSubcort, c("LH", "RH"),
                                           paste, sep = "\n")))),
    col_label = factor(col_label, levels = 
                         as.vector(t(outer(vNetSubcort, c("LH", "RH"), paste, sep = "\n")))),
    row_net = factor(row_net, levels = vNetSubcort),
    col_net = factor(col_net, levels = vNetSubcort)
  ) |>
  ggplot(aes(x = col_label, y = row_label, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(- 1, 1)) +
  facet_grid(row_net ~ col_net, space = "free", scales = "free", drop = FALSE) +
  labs(x = NULL, y = NULL, fill = TeX("$\\rho$")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 10),
    axis.text.y = element_text(angle = 0, hjust = 1, size = 10),
    strip.text.x = element_text(size = 14),
    strip.text.y = element_text(size = 14),
    panel.grid = element_blank(),
    panel.spacing = unit(0, "lines")
  ) +
  scale_x_discrete(labels = unique(df$col_hemi)) +
  scale_y_discrete(labels = unique(df$row_hemi))

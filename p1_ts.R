library(tidyverse)


################################################################################
#                                  IMPORT                                      #
################################################################################

library(foreach)
library(doParallel)

cl <- makeCluster(8)
registerDoParallel(cl)

cDirecS = "C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/brain_data/CSV/"
cDirecC = "C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/brain_data/HCP_Net_BOLD/"

vcFilenamesS = list.files(cDirecS)
vcFilenamesC = list.files(cDirecC)

vcNetworks = c("Vis", "SomMot", "DorsAttn", "SalVentAttn", "Limbic", "Cont", "Default")
vID = gsub("\\D", "", vcFilenamesS)
lAll =
  foreach(i = 1 : length(vcFilenamesS),
          .packages = c("tidyverse")) %dopar% {

            cID = vID[i]

            # import sottocort
            tInS =
              readr::read_csv(paste0(cDirecS, cID, ".nii.gz.csv"), show_col_types = FALSE) |>
              dplyr::select(- "...1")
            # if(i == 1){
            vcSottocort = colnames(tInS)[!str_detect(string = colnames(tInS),
                                                     pattern = str_c(vcNetworks, collapse = "|"))]
            # }

            # import cort
            tInC = readr::read_csv(paste0(cDirecC, cID, ".nii.gz.csvnetwork.csv"),
                                   show_col_types = FALSE) |>
              dplyr::select(- "...1")

            return(list("S" = tInS[, vcSottocort], "C" = tInC))
          }

stopCluster(cl)

tS = data.frame(matrix(NA, nrow = 0, ncol = length(colnames(lAll[[1]]$S))))
colnames(tS) <- colnames(lAll[[1]]$S)
tC = data.frame(matrix(NA, nrow = 0, ncol = length(colnames(lAll[[1]]$C))))
colnames(tC) <- colnames(lAll[[1]]$C)

for(i in 1 : length(lAll)){
  tS = rbind(tS,
             lAll[[i]]$S |>
               mutate(id = vID[i],
                      t = row_number()))
  tC = rbind(tC,
             lAll[[i]]$C |>
               mutate(id = vID[i],
                      t = row_number()))
}

colnames(tS)[1 : 16] = gsub("\\-(rh|lh)$", "_\\U\\1", colnames(lAll[[i]]$S), perl = TRUE)
colnames(tC)[1 : 14] =  sub("^(LH|RH)_(.*)", "\\2_\\1", gsub("_mean", "", colnames(lAll[[i]]$C)))

rm(lAll)


################################################################################
#                                  STUDY                                       #
################################################################################

tC |>
  select(- c("id", "t")) |> apply(MARGIN = 2, FUN = mean)

library(vars)
library(Spillover)

# plot series

tmax = 3600; nd = 77
tS |>
  filter(t < tmax,
         id == vID[nd]) |>
  pivot_longer(!c("t", "id"),
               values_to = "val",
               names_to = "reg") |>
  separate(reg, into = c("reg", "hemi")) |>
  rename(Hemisphere = hemi) |>
  mutate(id = as.factor(id),
         Hemisphere = as.factor(Hemisphere)) |>
  group_by(id, Hemisphere) |>
  ggplot(aes(x = t, y = val, col = Hemisphere, group = interaction(id, Hemisphere))) +
  geom_line() +
  facet_wrap(~ reg) +
  ylab("Value") +
  xlab("Time (seconds)") +
  ggtitle("Subcortical structures") +
  theme_minimal()

nd = 151
tC |>
  filter(t < tmax,
         id == vID[nd]) |>
  pivot_longer(!c("t", "id"),
               values_to = "val",
               names_to = "reg") |>
  separate(reg, into = c("reg", "hemi")) |>
  rename(Hemisphere = hemi) |>
  mutate(id = as.factor(id),
         Hemisphere = as.factor(Hemisphere)) |>
  group_by(id, Hemisphere) |>
  ggplot(aes(x = t, y = val, col = Hemisphere, group = interaction(id, Hemisphere))) +
  geom_line() +
  facet_wrap(~ reg) +
  ylab("Value") +
  xlab("Time (seconds)") +
  ggtitle("Cortical structures (aggregates)") +
  theme_minimal()

# stationarity test

library(tseries)

tS |> # h0: unit root, all rejected
  group_by(id) %>%
  summarise(across(starts_with("HIP") | starts_with("AMY") | starts_with("pTHA") | starts_with("aTHA") |
                     starts_with("NAc") | starts_with("GP") | starts_with("PUT") | starts_with("CAU"),
                   ~ adf.test(.x)$p.value,
                   .names = "adf_{.col}")) %>%
  ungroup() |>
  View()
  # write.csv("C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/draft/ts_analysis/adf_subcortical.csv")

tC |> # h0: unit root, all rejected
  group_by(id) %>%
  summarise(across(starts_with("VIS") | starts_with("SMN") | starts_with("DAN") | starts_with("VAN") |
                     starts_with("LIM") | starts_with("FPN") | starts_with("DMN"),
                   ~ adf.test(.x)$p.value,
                   .names = "adf_{.col}")) %>%
  ungroup() |>
  View()
  # write.csv("C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/draft/ts_analysis/adf_cortical.csv")

# fractality test

vOutS = c("178243", "198653", "115825", "550439", "562345", "706040",
          "169747", "158035", "581450", "412528")

source("C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2fou/functions/fWXY23.R")
tHS =
  tS |>
  pivot_longer(!c("id", "t"), names_to = "net", values_to = "val") |>
  arrange(id, net, t) |>
  group_by(id, net) |>
  summarise(H = H.hat(val),
            se = H.se(val)) |>
  mutate(
    t = (H - 0.5) / se,
    less = t < qnorm(0.05)
  )

tHS = read_csv("C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/draft/ts_analysis/cof_subcortical.csv")

tHS |>
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

tHS |>
  dplyr::select(id, net, less) |>
  filter(!id %in% vOutS) |>
  separate(net, into = c("net", "hemi"), sep = "_") |>
  group_by(net) |>
  summarise(
    rough = sum(less == TRUE),
    nonrough = sum(less == FALSE),
    proportion_rough = mean(less == TRUE),
    .groups = 'drop'
  ) # |>
  # write_csv("C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/draft/ts_analysis/cof_subcortical_summary.csv")


tHC =
  tC |>
  pivot_longer(!c("id", "t"), names_to = "net", values_to = "val") |>
  arrange(id, net, t) |>
  group_by(id, net) |>
  summarise(H = H.hat(val),
            se = H.se(val)) |>
  mutate(
    t = (H - 0.5) / se,
    less = t < qnorm(0.05),
    equal = abs(t) > qnorm(0.025)
  )

tHC = read_csv("C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/draft/ts_analysis/cof_cortical.csv")

vOutC = c("115825", "198653", "178243", "706040", "412528", "581450", "550439",
          "463040", "172130", "145834")
tHC |>
  dplyr::select(id, net, H) |>
  filter(!id %in% vOutC) |>
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

tHC |>
  dplyr::select(id, net, less) |>
  filter(!id %in% vOutC) |>
  separate(net, into = c("net", "hemi"), sep = "_") |>
  group_by(net) |>
  summarise(
    rough = sum(less == TRUE),
    nonrough = sum(less == FALSE),
    proportion_rough = mean(less == TRUE),
    .groups = 'drop'
  ) # |>
#   write_csv("C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/draft/ts_analysis/cof_cortical_summary.csv")

# comparison with MDE estimates of H
cIn = "C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/estimates/"
for(j in 1 : length(vID)){
  cID = vID[j]
  if(j == 1){
    tHS_post =
      read.csv(paste0(cIn, "sottocort/est_sottocort_l51015_", cID, ".csv")) |>
      filter(grepl("^H_", par))
    tHC_post =
      read.csv(paste0(cIn, "cort_aggr/est_cortaggr_l51015_", cID, ".csv")) |>
      filter(grepl("^H_", par))
  } else {
    tHS_post =
      tHS_post |>
      left_join(
        read.csv(paste0(cIn, "sottocort/est_sottocort_l51015_", cID, ".csv")) |>
          filter(grepl("^H_", par)),
        by = "par"
      )
    tHC_post =
      tHC_post |>
      left_join(
        read.csv(paste0(cIn, "cort_aggr/est_cortaggr_l51015_", cID, ".csv")) |>
          filter(grepl("^H_", par)),
        by = "par"
      )
  }

}
colnames(tHC_post) <- colnames(tHS_post) <- c("par", vID)

# tHS_post |>
tHC_post |>
  pivot_longer(!par, names_to = "id", values_to = "H") |>
  rename(net = par) |>
  mutate(net = gsub("^H_", "", net)) |>
  # filter(!id %in% vOutS) |>
  filter(!id %in% vOutC) |>
  # separate(net, into = c("net", "hemi"), sep = "_") |>
  mutate(net = gsub("_mean", "", net)) |>
  separate(net, into = c("hemi", "net"), sep = "_") |>
  group_by(net) |>
  ggplot(aes(x = H, fill = hemi)) +
  # geom_density(alpha = 0.7, adjust = 1.5, color = "black") +
  geom_histogram(binwidth = 0.03, alpha = 0.7, position = "identity", color = "black") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 0.5) +
  facet_wrap(~ net) +  # Creates a separate plot for each 'net'
  labs(x = latex2exp::TeX("\\hat{H}"),
       y = "Density",
       fill = "Hemisphere") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"))

tHS_post |>
  pivot_longer(!par, names_to = "id", values_to = "H") |>
  rename(net = par) |>
  mutate(net = gsub("^H_", "", net),
         id = as.integer(id)) |>
  left_join(tHS, by = c("net", "id")) |>
  select(net, id, H.x, H.y) |>
  rename(post = H.x, pre = H.y) |>
  mutate(net = str_sub(net, 1, - 4)) |>
  ggplot(aes(x = pre, y = post, color = net)) +
  geom_point() +
  ylim(0, 1) +
  xlim(0, 1) +
  ggtitle("Subcortical")

tHC_post |>
  pivot_longer(!par, names_to = "id", values_to = "H") |>
  rename(net = par) |>
  mutate(net = gsub("^H_|_mean", "", net),
         id = as.integer(id)) |>
  separate(net, into = c("part1", "part2"), sep = "_") |>
  unite(net, part2, part1, sep = "_") |>
  left_join(tHC, by = c("net", "id")) |>
  select(net, id, H.x, H.y) |>
  rename(post = H.x, pre = H.y) |>
  mutate(net = str_sub(net, 1, - 4)) |>
  ggplot(aes(x = pre, y = post, color = net)) +
  geom_point() +
  ylim(0, 1) +
  xlim(0, 1) +
  ggtitle("Cortical")

# corr lh-rh

t_corr_S =
  tS |>
  pivot_longer(!c("t", "id"),
               values_to = "val",
               names_to = "reg") |>
  separate(reg, into = c("reg", "hemi")) |>
  pivot_wider(id_cols = c("id", "t", "reg"),
              names_from = "hemi", values_from = "val") |>
  group_by(id, reg) |>
  summarise(rho = cor(rh, lh))

t_corr_C =
  tC |>
  pivot_longer(!c("t", "id"),
               values_to = "val",
               names_to = "reg") |>
  separate(reg, into = c("hemi", "reg")) |>
  pivot_wider(id_cols = c("id", "t", "reg"),
              names_from = "hemi", values_from = "val") |>
  group_by(id, reg) |>
  summarise(rho = cor(RH, LH))

# plot corr

t_corr_S |>
  ggplot(aes(x = reg, y = rho, col = id)) +
  geom_point() +
  ggtitle("Sottocorticali: correlazioni LH - RH") +
  theme(legend.position = "none")

t_corr_C |>
  ggplot(aes(x = reg, y = rho, col = id)) +
  geom_point() +
  ggtitle("Corticali: correlazioni LH - RH") +
  theme(legend.position = "none")


################################################################################
#                                OUTLIERS                                      #
################################################################################
# Campi, Marta, et al. "Signature Isolation Forest." arXiv preprint arXiv:2403.04405 (2024).

library(isotree)

# import

tSignS = read_csv("C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/brain_data/sign_d3_sottocort.csv")
tSignC = read_csv("C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/brain_data/sign_d3_cort.csv")

# outlier detection (consider standard deviations)

oIfS = isolation.forest(data = tSignS, output_score = TRUE)
v10iS = order(oIfS$scores, decreasing = TRUE)[1:10]
oIfC = isolation.forest(data = tSignC, output_score = TRUE)
v10iC = order(oIfC$scores, decreasing = TRUE)[1:10]

# plot scores

scores_S <- data.frame(Index = 1:length(oIfS$scores), Score = oIfS$scores)
top10_S <- data.frame(Index = v10iS, Score = oIfS$scores[v10iS], Label = vID[v10iS])

scores_C <- data.frame(Index = 1:length(oIfC$scores), Score = oIfC$scores)
top10_C <- data.frame(Index = v10iC, Score = oIfC$scores[v10iC], Label = vID[v10iC])

p1 <- ggplot(scores_S, aes(x = Index, y = Score)) +
  geom_point(color = "blue") +
  geom_point(data = top10_S, aes(x = Index, y = Score), color = "red") +
  geom_text(data = top10_S, aes(x = Index, y = Score, label = Label),
            position = position_nudge(y = 0.02), size = 3) +
  labs(title = "Subcortical", y = "Isolation Forest Score") +
  theme_minimal() # +
  # theme(axis.text = element_blank(), axis.title.x = element_blank(),
  #       panel.grid.major = element_line(color = "gray", linetype = "dashed"))

p2 <- ggplot(scores_C, aes(x = Index, y = Score)) +
  geom_point(color = "blue") +
  geom_point(data = top10_C, aes(x = Index, y = Score), color = "red") +
  geom_text(data = top10_C, aes(x = Index, y = Score, label = Label),
            position = position_nudge(y = 0.02), size = 3) +
  labs(title = "Cortical", y = "Isolation Forest Score") +
  theme_minimal() # +
  # theme(axis.text = element_blank(), axis.title.x = element_blank(),
  #       panel.grid.major = element_line(color = "gray", linetype = "dashed"))

library(gridExtra)
grid.arrange(p1, p2, ncol = 2)

# plot outlier time series

vOutS = vID[v10iS]

for(i in 1 : length(vOutS)){
  p =
    tS |>
    filter(id ==
             vID[77]) |>   # least outlier
             # vOutS[i]) |>
    pivot_longer(!c("t", "id"),
                 values_to = "val",
                 names_to = "reg") |>
    separate(reg, into = c("reg", "hemi")) |>
    rename(Hemisphere = hemi) |>
    mutate(id = as.factor(id),
           Hemisphere = as.factor(Hemisphere)) |>
    group_by(id, Hemisphere) |>
    ggplot(aes(x = t, y = val, col = Hemisphere, group = interaction(id, Hemisphere))) +
    geom_line() +
    facet_wrap(~ reg) +
    ylab("Value") +
    xlab("Time (seconds)") +
    ggtitle(paste0("Subcortical structures (ID: ", vOutS[i], ")")) +
    theme_minimal()
  print(p)
}

vOutC = vID[v10iC]

for(i in 1 : length(vOutC)){
  p =
    tC |>
    filter(id ==
             # vID[151]) |>   # least outlier
             vOutC[i]) |>
    pivot_longer(!c("t", "id"),
                 values_to = "val",
                 names_to = "reg") |>
    separate(reg, into = c("reg", "hemi")) |>
    rename(Hemisphere = hemi) |>
    mutate(id = as.factor(id),
           Hemisphere = as.factor(Hemisphere)) |>
    group_by(id, Hemisphere) |>
    ggplot(aes(x = t, y = val, col = Hemisphere, group = interaction(id, Hemisphere))) +
    geom_line() +
    facet_wrap(~ reg) +
    ylab("Value") +
    xlab("Time (seconds)") +
    ggtitle(paste0("Cortical structures (ID: ", vOutC[i], ")")) +
    theme_minimal()

  print(p)
}

# outliers characteristics

vOutliers = unique(c(vID[v10iC], vID[v10iS]))

readxl::read_excel("C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/brain_data/cognitive_variables.xlsx") |>
  rename(crist = FA1, fluid = FA2) |>
  left_join(
    y =
      read_delim("C:/Users/dugor/OneDrive - Universita' degli Studi di Roma Tor Vergata/phd/2. neuro/brain_data/volatility_brain_project.csv",
                 delim = ";") |>
      rename(ID = Subject),
    by = "ID"
  ) |>
  mutate(outlier = ID %in% unique(c(vID[v10iC], vID[v10iS]))) |>
  select(- ID) |>
  group_by(outlier) |>
  summarise(
    across(where(is.numeric),
           list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE)),
    .names = "{.col}_{.fn}"),
    .groups = "drop"
    ) |>
  pivot_longer(!outlier, names_to = "var", values_to = "val") |>
  mutate(
    type = ifelse(grepl("mean", var), "mean", "sd"),
    var = gsub("_mean|_sd", "", var)
  ) |>
  # pivot_wider(id_cols = c("outlier", "var"), names_from = type, values_from = val) |>
  pivot_wider(id_cols = "var", names_from = c(type, outlier), values_from = val) |>
  mutate(t_stat =
           (mean_FALSE - mean_TRUE) /
           sqrt(sd_FALSE ^ 2 / 13 + sd_TRUE ^ 2 / (173 - 13)),
         nu = (sd_FALSE ^ 2 / 13 + sd_TRUE ^ 2 / (173 - 13)) ^ 2 /
           (sd_FALSE ^ 4 / (13 ^ 2 * 12) + sd_TRUE ^ 4 / (160 ^ 2 * 159)),
         pval = 2 * (1 - pt(q = abs(t_stat), df = nu)),
         h0 = pval > 0.05)



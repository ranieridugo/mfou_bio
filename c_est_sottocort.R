source("/home/dugo/brain/functions/f_estimator_mde.R")
vcLibraries = c("readr", "dplyr", "stringr", "tibble", "tidyr", "parallel")
lapply(vcLibraries, require, character.only = TRUE)
vcFilenames = list.files("/home/dugo/brain/data/CSV")
ncores = 84

f_namedvec_to_df <- function(vec){
  data.frame(
    par = 
      if(is.null(names(vec))){"na"} else {names(vec)},
    val = vec
  )
}

for(lag in c(0)){ # c(- 1, 0, 1, 2, 3)){
  iLMax = 5 + lag; vLAdd = c(10, 15) + lag
  
  f_iteration = 
    function(i){
      cFile = vcFilenames[i]
      cID = gsub("\\D", "", cFile)
      
      tIn =
        readr::read_csv(paste0("/home/dugo/brain/data/CSV/", cID, ".nii.gz.csv"), show_col_types = FALSE) |>
        select(- "...1") |> 
        as.data.frame()
      colnames(tIn) = gsub("\\-(rh|lh)$", "_\\U\\1", colnames(tIn), perl = TRUE)
      
      vcNetworks = c("Vis", "SomMot", "DorsAttn", "SalVentAttn", "Limbic", "Cont", "Default")
      vNames = colnames(tIn)
      vcSottocort = vNames[!str_detect(string = vNames, pattern = str_c(vcNetworks, collapse = "|"))]
      
      vSottocort <-
        tryCatch(
          f_mde(mX = tIn[, vcSottocort], delta = 1, iLMax = iLMax, vLAdd = vLAdd,
                type = "cau"),
          error = function(e) NA  
        )
      
      write_csv(f_namedvec_to_df(vSottocort), 
                paste0("/home/dugo/brain/estimates/est_sottocort_l", 
                       iLMax, paste0(vLAdd, collapse = ""), "_", cID, ".csv"))
      
    }
  
  mclapply(c(1 : length(vcFilenames)), f_iteration, mc.cores = ncores)
}

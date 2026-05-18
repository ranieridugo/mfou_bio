source("/home/dugo/brain/functions/f_estimator_mde.R")
vcLibraries = c("readr", "dplyr", "stringr", "tibble", "tidyr", "parallel")
lapply(vcLibraries, require, character.only = TRUE)
vcFilenames = list.files("/home/dugo/brain/data/cort_aggr/")
ncores = 84

f_namedvec_to_df <- function(vec){
  data.frame(
    par = 
      if(is.null(names(vec))){"na"} else {names(vec)},
    val = vec
  )
}

for(lag in c(- 1, 1, 2, 3, 4, 5)){
  iLMax = 5 + lag; vLAdd = c(10, 15) + lag
  
  f_iteration = 
    function(i){
      cFile = vcFilenames[i]
      cID = gsub("\\D", "", cFile)
      
      tIn =
        readr::read_csv(paste0("/home/dugo/brain/data/cort_aggr/", cFile),
                        show_col_types = FALSE) |>
        select(- "...1") |> 
        as.data.frame()
      
      vOut <-
        tryCatch(
          f_mde(mX = tIn, delta = 1, iLMax = iLMax, vLAdd = vLAdd, type = "cau"),
          error = function(e) NA  
        )
      
      write_csv(f_namedvec_to_df(vOut),
                paste0("/home/dugo/brain/estimates/cort_aggr/est_cortaggr_l",
                       paste0(iLMax, paste0(vLAdd, collapse = ""), collapse = ""),
                       "_", cID, ".csv"))
    }
  
  mclapply(c(1 : length(vcFilenames)), f_iteration, mc.cores = ncores)
}
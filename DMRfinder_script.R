# Arguments -------------------------------------------------------------------
option_list = list(
  optparse::make_option(opt_str = c("-d", "--datadir"),
                        default = "data",
                        help = "data directory [default= %default]"),
  optparse::make_option(opt_str = c("-o", "--outdir"),
                        default = "results",
                        help = "output directory [default= %default]"),
  optparse::make_option(opt_str = c("-n", "--cores"),
                        default = 12,
                        help = "Number of threads [default= %default]"),
  optparse::make_option(opt_str = c("-l", "--log"),
                        default = "log",
                        help = "log directory [default= %default]"),
  optparse::make_option(opt_str = c("-m", "--metadata"),
                        default = "metadata.xlsx",
                        help = "path to metadata excel sheet [default= %default]"),
  optparse::make_option(opt_str = c("-g", "--grouping"),
                        default = "Condition1",
                        help = "Condition to group the samples by in the metadata-sheet [default= %default]")
)

opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

# Packages installation/loading -------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
})


# Setup  ------------------------------------------------------------------
future::plan(future::multisession, workers = opt$cores)
dir.create(path = opt$log, showWarnings = F)
dir.create(path = opt$outdir, showWarnings = F)

# Run -------------------------------------------------------------------------
files_w_full_path <- list.files(path = opt$datadir,
                                pattern=".bam", 
                                full.names = T)


files_simplified <- files_w_full_path %>% 
  basename() %>% 
  str_remove(pattern="_S[0-9]{1,3}_R1_001.UMI_trimmed.fq_trimmed_bismark_bt2.deduplicated.bam")
  #str_remove(pattern="_S[0-9]_R1_001.UMI_trimmed.fq_trimmed_bismark_bt2.deduplicated.bam")

## Generate coverage files --------------------------------------------------
print("Generating coverage files...")
furrr::future_map2(.x = files_w_full_path, 
                   .y=files_simplified, 
                   .f =~system2(command = "samtools",
                                args = c("view", "-h" ,.x, "|", 
                                         "python3 extract_CpG_data.py",
                                         "-i -", 
                                         "-o", file.path(opt$outdir,paste0(.y, ".cov")),
                                         "-v"),
                                stdout = file.path(opt$log, paste0(.y, "_cov_out.txt")),
                                stderr = file.path(opt$log, paste0(.y, "_cov_err.txt"))
                                
                   )
)


## Combine CpG sites -------------------------------------------------------
print("Combining CpG sites...")
system2(command="python3",
        args= c("combine_CpG_sites.py -v -o", file.path(opt$outdir,"combined.csv"), 
                file.path(opt$outdir, paste0(basename(files_simplified), ".cov"))),
        stdout = file.path(opt$log, "combine_out.txt"),
        stderr = file.path(opt$log, "combine_err.txt"))


## Find DMRS ---------------------------------------------------------------

meta_data_sheet <- readxl::read_xlsx(opt$metadata)
treatments <- meta_data_sheet %>%
  pull(opt$grouping) %>% unique()

# pull IDS for each treatment
grp_1_IDs <- meta_data_sheet %>% 
  filter(.data[[opt$grouping]]==treatments[1]) %>% 
  pull(`Sample ID (Library ID)`)

grp_2_IDs <- meta_data_sheet %>% 
  filter(.data[[opt$grouping]]==treatments[2]) %>% 
  pull(`Sample ID (Library ID)`)

print("Finding DMRs...")
system2(command="Rscript",
        args=c("findDMRs.r",
               "-i", file.path(opt$outdir, "combined.csv"),
               "-o", file.path(opt$outdir, "results.csv"),
               "-v", 
               "-n",
               "Control,Exptl",
               paste(grp_1_IDs, collapse = ","),
               paste(grp_2_IDs, collapse = ",")),
        stdout = file.path(opt$log, "findDMRS_out.txt"),
        stderr = file.path(opt$log, "findDMRS_err.txt"))

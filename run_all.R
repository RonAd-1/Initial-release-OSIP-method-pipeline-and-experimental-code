# run_all.R - Master Entry Point Pipeline

cat("=================================================================\n")
cat("                OSIP Execution Pipeline                          \n")
cat("=================================================================\n\n")

# Flag file used to track whether environment setup has been run
SETUP_FLAG <- ".setup_complete"

# 1. Determine if this is the initial run
if (!file.exists(SETUP_FLAG)) {
  cat("[1/2] Initial run detected. Running environment setup...\n\n")
  
  # Run dependency check and installation
  source("setup.R")
  
  # Create the hidden marker file to signify initial setup is finished
  writeLines(as.character(Sys.time()), SETUP_FLAG)
  cat("\n[✓] Initial setup complete. Flag file '.setup_complete' created.\n\n")
} else {
  cat("[1/2] Setup already completed on this system. Skipping 'setup.R'.\n")
  cat("      (Delete '.setup_complete' if you need to force re-running setup)\n\n")
}

# 2. Execute main comparison analysis
cat("[2/2] Running main analysis pipeline...\n\n")
source("main_comparison_units_borders.R")

cat("\n=================================================================\n")
cat("                Pipeline Execution Finished                      \n")
cat("=================================================================\n")
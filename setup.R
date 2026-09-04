# setup.R - One-time dependency installation and version validation script

cat("=================================================================\n")
cat("          OSIP Environment & Dependency Setup                     \n")
cat("=================================================================\n\n")

# 1. Verify requirements.txt exists
req_file <- "requirements.txt"
if (!file.exists(req_file)) {
  stop("Error: 'requirements.txt' file not found in working directory.")
}

# 2. Read and parse requirements.txt
lines <- readLines(req_file, warn = FALSE)
lines <- trimws(lines)
lines <- lines[lines != "" & !startsWith(lines, "#")]

parse_requirement <- function(line) {
  # Match pattern like "MatchIt (>= 4.7.2)" or simple "MatchIt"
  pattern <- "^([a-zA-Z0-9.]+)\\s*(?:\\(>=\\s*([0-9.-]+)\\))?$"
  matches <- regmatches(line, regexec(pattern, line))[[1]]
  
  if (length(matches) == 0) {
    warning("Could not parse requirement line: ", line)
    return(NULL)
  }
  
  list(
    package = matches[2],
    min_version = if (matches[3] != "") matches[3] else NULL
  )
}

# To this (Base R - works with zero dependencies loaded):
requirements <- Filter(Negate(is.null), lapply(lines, parse_requirement))

# 3. Check and install dependencies
to_install <- character(0)
to_update  <- character(0)

installed_pkgs <- installed.packages()[, "Package"]

for (req in requirements) {
  pkg <- req$package
  min_ver <- req$min_version
  
  if (!(pkg %in% installed_pkgs)) {
    to_install <- c(to_install, pkg)
  } else if (!is.null(min_ver)) {
    curr_ver <- packageVersion(pkg)
    if (curr_ver < numeric_version(min_ver)) {
      cat(sprintf("[-] Package '%s' is outdated (Installed: %s, Required: >= %s)\n", 
                  pkg, curr_ver, min_ver))
      to_update <- c(to_update, pkg)
    }
  }
}

pkgs_needed <- unique(c(to_install, to_update))

# 4. Execute Installation
if (length(pkgs_needed) > 0) {
  cat(sprintf("\nFound %d package(s) requiring installation or update.\n", length(pkgs_needed)))
  cat("Packages:", paste(pkgs_needed, collapse = ", "), "\n\n")
  
  install.packages(pkgs_needed, repos = "https://cloud.r-project.org")
  
  cat("\n[✓] Dependency check and installation complete!\n")
} else {
  cat("\n[✓] All required packages are installed and up to date!\n")
}

cat("=================================================================\n")
# setup.R

# setup.R - One-time dependency installation script

cat("Checking project dependencies against requirements.txt...\n")

if (!file.exists("requirements.txt")) {
  stop("requirements.txt file not found in the root directory!")
}

# Read package names from file
required_packages <- readLines("requirements.txt")
required_packages <- trimws(required_packages)
required_packages <- required_packages[required_packages != "" & !startsWith(required_packages, "#")]

# Determine which packages are missing
missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]

# Install missing packages
if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
  cat("All missing dependencies successfully installed!\n")
} else {
  cat("All required packages are already installed.\n")
}

## Quick Start

### Prerequisites
Ensure you have **R (≥ 4.2.0)** and **RStudio** installed on your system.

### Installation & Execution

**1. Clone the Repository**  
Clone this repository to your local machine and navigate into the project directory.

**2. Run the Full comparison**  
Open RStudio or an R terminal in the project root directory and execute:

```R
source("run_all.R")
```

> **Note:** On the initial run, `run_all.R` automatically reads `requirements.txt`, verifies and installs missing CRAN dependencies via `setup.R`, and creates a local `.setup_complete` marker file. All subsequent runs will automatically bypass the setup step and execute the main pipeline (`main_comparison_units_borders.R`) directly.

If for you are interested in the OSIP only comparison, you can run instead:

```R
source("run_osip_only.R")
```

In both cases, you can either run both strict and robust modes, 
or only one of them.

---

## Modifying the Program for Your Needs

If you want to experiment with or extend the comparison framework, there are two primary files you should edit: `dataset_configs.R` and `main_comparison_units_borders.R`.

### 1. `dataset_configs.R`
This file contains all configuration parameters and macros used throughout the execution pipeline, including OSIP framework parameters, dataset-specific settings, and naming conventions.

Common modifications include:

* **Adding a New Dataset:** To run the comparison on a new dataset, update this file by registering the new dataset similarly to existing ones. Add its name to the `DATASET_CHOICES` list and define a corresponding `CONFIG` list.
* **Adding a New Method:** To integrate a custom matching or partitioning method, update all relevant macro objects (e.g., `Methods`, `METHOD_LABELS`).

### 2. `main_comparison_units_borders.R`
This is the core script driving the analysis pipeline. It handles importing dependencies, executing all comparison methods (including OSIP), and aggregating the results into tables and diagnostic plots.

The pipeline generates three primary outputs to monitor and evaluate results:

1. **General Output Log (`output_path`):**  
   Contains complete execution details, intermediate stages, and diagnostic sanity checks. By default, output files are named using the pattern:  
   `main_<dataset_name>_output_delta_<delta_value>_power_<dist_power>.txt`  
   The target directory and file name can be modified via the `output_path` variable at the top of the file.

2. **Console Output & Diagnostic Plots:**  
   Key execution metrics and plots are printed directly to the active console and graphics device for real-time tracking.

3. **Structured Summary Output Files:**  
   Generates dedicated output files for post-hoc analysis, such as overall performance comparisons (`consolidated_results`) and balance metrics (`love_plot_table`).

> **Note:** When adding a new dataset or method, minor updates may also be required in helper functions called by `main_comparison_units_borders.R` (e.g., `load_and_prep_data`, `compare_methods`, `complete_comparison`).

---

### 3. Adding a New Package Dependency
If your custom analysis or method relies on a new R package:

1. **Update `requirements.txt`:** Add the package name and optional minimum version (e.g., `pkgName (>= 1.0.0)`).
2. **Re-run Setup:** Execute `source("setup.R")` (or delete `.setup_complete` and run `source("run_all.R")`) to install the new package across environment instances.
3. **Import Package:** Add `library(pkgName)` to the **1. Load Required Libraries** section at the top of `main_comparison_units_borders.R`.
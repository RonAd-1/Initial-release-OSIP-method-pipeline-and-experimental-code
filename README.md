# Initial-release-OSIP-method-pipeline-and-experimental-code
This code is attached to the paper "Optimal Study Interior Partitioning (OSIP): An Optimization Approach to Observational Study Design under Fine Balance Constraints" by Ron Adar and Asaf Levin

## Quick Start

### Prerequisites
Ensure you have **R (≥ 4.2.0)** and **RStudio** installed on your system.

---

### Installation & Execution

#### 1. Clone the Repository
Clone this repository to your local machine:
```bash
git clone [https://github.com/RonAd-1/Initial-release-OSIP-method-pipeline-and-experimental-code.git)
cd your-osip-repo

#### 2. Install Dependencies (One-Time Setup)
Open RStudio or an R terminal in the project root directory and run the setup script. 
This reads requirements.txt and automatically installs any missing CRAN dependencies:

source("setup.R")

#### 3.Running the Analysis
Once the dependencies are installed, execute the main analysis pipeline:

R
source("main_comparison_units_borders.R")
# --- dataset_configs.R ---
                       
# In the following PARAMS list there are two fields you may want to edit: 

# -----------------------------------------------------------------------------
# 1. DELTA_VALUES: Specifies the delta threshold(s) for the OSIP framework.
# -----------------------------------------------------------------------------
# Option A: Single value
# DELTA_VALUES <- c(0.20)

# Option B: Regular sequence with constant step size
# DELTA_VALUES <- seq(0.15, 0.20, by = 0.05)

# Option C: Custom sequence with non-constant step sizes
# DELTA_VALUES <- c(0.15, 0.17, 0.22, 0.24)

# 2. SEARCHES_TO_USE: If you want to use only a subset of the heuristics you 
# can specify that here. The default is using all three. 

# All other fields are either technical parameters (e.g., DEBUG_STATUS) or 
# settings related to downsampling larger datasets. 
# None of these were used for the datasets in the paper. 
# MAX_VAL_FOR_PLOT is a technical parameter that restricts 
# the upper limit of a plot's axis you do not need to modify it unless you
# encounter a specific visualization issue.


PARAMS <- list(
  DELTA_VALUES = c(0.20),
  SEARCHES_TO_USE = c("local_search", "enhanced_ls", "simulated_annealing"),
  SEED_IN <- 25,
  SAMPLE = c(TRUE,FALSE),
  SAMPLE_SIZE = 200,
  N_CONTROL_SAMPLE = 750,
  DEBUG_STATUS = c(TRUE, FALSE),
  MAX_VAL_FOR_PLOT = NULL
  )

# Define the Enum-like object
Methods <- list(
  unmatched       = "unmatched",
  cem             = "cem",
  full_matching   = "full_matching",
  optimal_1_1     = "optimal_1_1",
  cardinality_1_1 = "cardinality_1_1",
  genetic_1_1     = "genetic_1_1",
  
  osip_step1_strict = "osip_step1_strict",  
  osip_step2_strict = "osip_step2_strict",
  osip_step2_strict_balanced = "osip_step2_strict_balanced",
  
  osip_step1_robust = "osip_step1_robust",  
  osip_step2_robust = "osip_step2_robust",
  osip_step2_robust_balanced = "osip_step2_robust_balanced",
  
  quintile        = "quintile",
  refined_quintile = "refined_quintile"
)

# Define labels as a simple vector, then assign the Enum values as names
METHOD_LABELS <- setNames(
  c("Unmatched", "CEM", "Optimal (1:1)", "Cardinality (1:1)", "OSIP_Step1_Strict", "OSIP_Step2_Strict",
    "OSIP_Step2_Strict_Balanced", 
    "OSIP_Step1_Robust", "OSIP_Step2_Robust","OSIP_Step2_Robust_Balanced",
    "Genetic (1:1)", "Full_matching", "Quintile (1:1)", "Refined Quintile (1:1)"),
  c(Methods$unmatched, Methods$cem, Methods$optimal_1_1, Methods$cardinality_1_1, 
    Methods$osip_step1_strict, Methods$osip_step2_strict, Methods$osip_step2_strict_balanced, 
    Methods$osip_step1_robust, Methods$osip_step2_robust, Methods$osip_step2_robust_balanced,
    Methods$genetic_1_1, Methods$full_matching, Methods$quintile, Methods$refined_quintile)
)

METHODS_1_TO_1_MATCHING_OSIP <- c(
  METHOD_LABELS[Methods$osip_step1_strict],
  METHOD_LABELS[Methods$osip_step2_strict],
  METHOD_LABELS[Methods$osip_step1_robust],
  METHOD_LABELS[Methods$osip_step2_robust]
)

# 1. Define Groups using the Enum Keys to pull the current Labels
METHODS_1_TO_1_MATCHING <- c(
  METHOD_LABELS[Methods$optimal_1_1],
  METHOD_LABELS[Methods$cardinality_1_1],
  METHOD_LABELS[Methods$genetic_1_1],
  METHOD_LABELS[Methods$osip_step2_strict],
  METHOD_LABELS[Methods$osip_step2_robust],
  METHOD_LABELS[Methods$quintile],
  METHOD_LABELS[Methods$refined_quintile]
)

METHODS_NON_1_TO_1_MATCHING <- c(
  METHOD_LABELS[Methods$unmatched],
  METHOD_LABELS[Methods$cem],
  METHOD_LABELS[Methods$full_matching]
)

# CONSTANTS for Dataset Selection
DATASET_CHOICES <- list(
  LALONDE = "lalonde",
  NSW_MIXTAPE = "nsw_mixtape",
  RHC = "rhc",
  LINDNER = "lindner",
  JOBS = "jobs",
  IDHP = "idhp",
  NHEFS = "nhefs"
)

# JOBS_CONFIG

JOBS_CONFIG <- list(
  TREATMENT_VAR      = "treat",
  OUTCOME_VAR        = "depress2",
  ID_VAR = "id",
  NUMERIC_COVARIATES = c("sex", "marital", "nonwhite"), # Adding 'occp' here if you want it,
  FACTOR_COVARIATES  = c("age", "educ", "income", "econ_hard", "depress1"),
  ALL_COVARIATES     =  c("age", "educ", "income", "econ_hard", "depress1", "sex", "marital", "nonwhite")
)

# IDHP_CONFIG 

# 3. Configuration for the IDHP dataset (X1-X25 format)
IDHP_CONFIG <- list(
  # Core variables
  TREATMENT_VAR = "treatment",  # adjust if different (e.g., "treat")
  OUTCOME_VAR = "y_factual",    # the observed outcome
  ID_VAR = "id", 
  
  # Covariates for matching (X1-X25)
  ALL_COVARIATES = paste0("X", 1:25),
  
  # CEM_SUBSET_COVARIATES = c(paste0("X", 1:6), "X7", "X8"),
  
  CEM_SUBSET_COVARIATES = paste0("X", 1:6),
  
  # Subsets of covariates based on type
  NUMERIC_COVARIATES = paste0("X", 1:6),
  FACTOR_COVARIATES = paste0("X", 7:25),
  
  # Additional outcome/evaluation variables (not used in matching, but for evaluation)
  COST_VAR = "c_factual",   # observed cost
  TRUE_OUTCOME_CONTROL = "mu0",  # true potential outcome under control
  TRUE_OUTCOME_TREATED = "mu1"   # true potential outcome under treatment
)

LALONDE_CONFIG <- list(
  # Core variables
  TREATMENT_VAR = "treat",
  OUTCOME_VAR   = "re78",
  ID_VAR        = "id", 
  
  # Covariates - Matches names in image_14a73a.png
  # Note: 'race' is a single factor here, not split into black/hisp
  ALL_COVARIATES = c("age", "educ", "race", "married", "nodegree", "re74", "re75"),
  
  # Subsets for standardizing and modeling
  NUMERIC_COVARIATES = c("age", "educ", "re74", "re75"),
  
  # 'race' is explicitly a factor with 3 levels in this version
  FACTOR_COVARIATES = c("race", "married", "nodegree")
)

# 2. Configuration for the NSW Mixtape dataset
# (Assumes 'race' is replaced by 'black' and 'hispan' 0-1 columns)
NSW_MIXTAPE_CONFIG <- list(
  # Core variables - often the same
  TREATMENT_VAR = "treat",
  OUTCOME_VAR = "re78",
  ID_VAR = "id", 
  
  # Covariates: 'race' is replaced by 'black' and 'hisp'
  ALL_COVARIATES = c("age", "educ", "black", "hisp", "marr", "nodegree", "re74", "re75"),
  
  # Subsets of covariates based on type
  # Note: 'black' and 'hispan' are often treated as numeric (0/1) or factors depending on the analysis. 
  # Here we'll treat them as numeric 0/1 for simplicity in the NUMERIC list.
  NUMERIC_COVARIATES = c("age", "educ", "black", "hisp", "re74", "re75"),
  FACTOR_COVARIATES = c("marr", "nodegree") # Removed 'race', added 'black' and 'hisp' to numeric
)

# ============================================================
# Configuration for the Right Heart Catheterization (RHC) dataset
# ============================================================

RHC_CONFIG <- list(
  
  # ... (Core variables remain the same) ...
  TREATMENT_VAR = "swang1",     
  OUTCOME_VAR   = "dth30",      
  ID_VAR        = "id",         
  
  # ----------------------------------------------------------
  # Updated ALL_COVARIATES (Uses the DUMMY variables)
  # ----------------------------------------------------------
  ALL_COVARIATES = c(
    # Demographics
    "age", "sex_male", "race_black", # <--- **FIXED NAMES HERE**
    
    # Vitals / Continuous
    "meanbp1", 
    
    # Comorbidities (History)
    "chfhx", "malighx"
  ),
  
  # ----------------------------------------------------------
  # Updated NUMERIC_COVARIATES
  # ----------------------------------------------------------
  NUMERIC_COVARIATES = c(
    "age", 
    "meanbp1" 
  ),
  
  # ----------------------------------------------------------
  # Updated FACTOR_COVARIATES (The binary dummies)
  # ----------------------------------------------------------
  FACTOR_COVARIATES = c(
    "sex_male",      
    "race_black",    
    "chfhx",         
    "malighx"        
  )
)

LINDNER_CONFIG <- list(
  NAME          = "lindner",
  TREATMENT_VAR = "abcix",        # 1 = Abciximab, 0 = Usual Care
  OUTCOME_VAR   = "log_cardbill",  # Log of cardiac-related costs
  ID_VAR        = "id",           # You may need to create this: lindner$id <- 1:nrow(lindner)
  
  # ----------------------------------------------------------
  # ALL_COVARIATES (Dummies and Continuous)
  # ----------------------------------------------------------
  ALL_COVARIATES = c(
    "stent",    # Binary: Stent deployed
    "height",   # Numeric
    "female",   # Binary
    "diabetic", # Binary
    "acutemi",  # Binary: Acute MI within 7 days
    "ejecfrac", # Numeric: Ejection fraction
    "ves1proc"  # Numeric: Vessels involved
  ),
  
  # ----------------------------------------------------------
  # NUMERIC_COVARIATES
  # ----------------------------------------------------------
  NUMERIC_COVARIATES = c(
    "height", 
    "ejecfrac", 
    "ves1proc"
  ),
  
  # ----------------------------------------------------------
  # FACTOR_COVARIATES (Binary Dummies)
  # ----------------------------------------------------------
  FACTOR_COVARIATES = c(
    "stent", 
    "female", 
    "diabetic", 
    "acutemi"
  )
)


NHEFS_CONFIG <- list(
  # Core variables
  TREATMENT_VAR = "qsmk",       # 1 = quit smoking, 0 = did not quit
  OUTCOME_VAR   = "wt82_71",    # Weight change (kg) between 1971 and 1982
  ID_VAR        = "id",       # NHEFS unique participant identifier
  
  # Covariates — standard set from Hernán & Robins (2020)
  # Captures demographics, smoking history, baseline health behaviours, and weight
  ALL_COVARIATES = c(
    "sex", "race", "age", "school",
    "smokeintensity", "smokeyrs",
    "exercise", "active",
    "wt71"
  ),
  
  # Continuous covariates
  # Note: smokeintensity, smokeyrs, and wt71 are right-skewed —
  # strong candidates for robust Mahalanobis based on cross-dataset findings
  NUMERIC_COVARIATES = c("age", "school", "smokeintensity", "smokeyrs", "wt71"),
  
  # Ordered/binary factors
  # exercise: 0=much, 1=moderate, 2=little  — treat as factor (3 levels)
  # active:   0=very, 1=moderate, 2=inactive — treat as factor (3 levels)
  # sex and race are binary (0/1)
  FACTOR_COVARIATES = c("sex", "race", "exercise", "active")
)
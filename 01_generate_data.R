# Generate the simulation data used in the paper (Rscript 01_generate_data.R).
# Install OpenMx and expm first. Run from this repository's root directory.
# CONTROLS: edit the defaults below, or pass --effects=Weak,Mid,Large,
# --replications=1:500, --workers=1, --output-dir=generated, --overwrite=false,
# --retain-only=false. Default saves the original full trajectories.
# A small run is: Rscript 01_generate_data.R --effects=Weak --replications=1
# Default output: generated/Simulations_<effect>/{Dataset,Trues} <rep>.RDS.
# All 100000 observations are generated for each of 60 subjects. Optionally set
# --retain-only=true to save only the 200 observations used at each interval
# (.5, 1, 5): 460 rows per subject, preserving every observation used for fitting.
# Optional retention uses about 1.7 MB per dataset (2.5 GB total); default full
# trajectories use 360 MB per dataset (540 GB total). Allow 1 GB RAM per worker.
# Retention saves storage, not simulation time. Use full mode if changing the
# fitting intervals or series length; existing files need --overwrite=true.
# Archived truth matrices are included because their original RNG seed was not
# recorded. Subject-level RNG, burn-in, time indexing, and noise rounding match
# the original generator. Do not change them when reproducing the study.
# To create a NEW study, use --new-truth=true --truth-seed=12345 (changes results).
# Functions can also be sourced without running the simulation.

defaults <- list(effects = c("Weak", "Mid", "Large"), replications = 1:500,
                 workers = 1L, output_dir = "generated", overwrite = FALSE,
                 retain_only = FALSE, new_truth = FALSE, truth_seed = 12345L)

read_options <- function(args = commandArgs(TRUE)) {
  opt <- defaults
  for (arg in args) {
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- gsub("-", "_", pair[1], fixed = TRUE)
    if (length(pair) != 2L || !key %in% names(opt)) stop("Unknown option: ", arg)
    value <- pair[2]
    opt[[key]] <- switch(key,
      effects = strsplit(value, ",", fixed = TRUE)[[1]],
      replications = {
        if (!grepl("^[0-9]+(:[0-9]+)?$", value)) stop("Use --replications=1 or 1:500")
        bounds <- as.integer(strsplit(value, ":", fixed = TRUE)[[1]])
        if (length(bounds) == 1L) bounds else seq.int(bounds[1], bounds[2])
      },
      workers = as.integer(value), truth_seed = as.integer(value),
      overwrite = match.arg(value, c("true", "false")) == "true",
      retain_only = match.arg(value, c("true", "false")) == "true",
      new_truth = match.arg(value, c("true", "false")) == "true", value)
  }
  stopifnot(all(opt$effects %in% c("Weak", "Mid", "Large")),
            all(opt$replications %in% 1:500), length(opt$workers) == 1L,
            is.finite(opt$workers), opt$workers >= 1L)
  opt
}

# Original stationary group/subgroup drift generator; preserved for new studies.
create.amat <- function(nvar = 6, prop = .20, allow.neg = TRUE, EF = .30,
                        EF.sd = .03, SEF = 1, subgroups = 2, sub.prop = .15) {
  repeat {
    A.mat <- diag(rnorm(nvar, -.50, .05), nvar)
    indices <- which(A.mat == 0)
    selected <- sample(indices, prop * length(indices))
    mult <- if (allow.neg) ifelse(rbinom(length(selected), 1, .50) == 1, 1, -1) else 1
    A.mat[selected] <- rnorm(length(selected), EF * mult, EF.sd)
    if (max(Re(eigen(A.mat)$values)) < -.20) break
  }
  if (subgroups <= 1) return(A.mat)
  A.mats <- list()
  GA.mat <- A.mat
  G.indices <- which(GA.mat == 0)
  n.paths <- floor(sub.prop * length(G.indices))
  for (j in seq_len(subgroups)) {
    if (n.paths > length(G.indices)) stop("Not enough unique subgroup paths.")
    repeat {
      A.mats[[j]] <- GA.mat
      selected <- sample(G.indices, n.paths)
      mult <- if (allow.neg) ifelse(rbinom(length(selected), 1, .50) == 1, 1, -1) else 1
      A.mats[[j]][selected] <- rnorm(length(selected), SEF * mult, EF.sd)
      if (max(Re(eigen(A.mats[[j]])$values)) < -.20) {
        G.indices <- G.indices[!G.indices %in% selected]
        break
      }
    }
  }
  A.mats[[subgroups + 1L]] <- GA.mat
  A.mats
}

# The five-decimal rounding is part of the original data-generating process.
compute_dtQ_identity <- function(A, dt) {
  I <- diag(1, nrow(A))
  ahash <- kronecker(A, I) + kronecker(I, A)
  K <- expm::expm(ahash * dt) - diag(1, nrow(ahash))
  round(matrix(solve(ahash, K %*% as.vector(I)), nrow(A), ncol(A)), 5)
}

simulate_subject <- function(sim_id, xp, A_ct, deltat = .01, times = 100000L,
                             retain_only = FALSE) {
  # Original PSOCK workers used clusterSetRNGStream before this per-subject seed.
  RNGkind("L'Ecuyer-CMRG")
  set.seed(1e6 * sim_id + xp)
  ne <- nrow(A_ct)
  M <- OpenMx::mxMatrix
  osc <- OpenMx::mxModel("OUMod",
    M("Full", ne, ne, free = TRUE, values = expm::expm(A_ct * deltat), name = "A"),
    M("Zero", ne, ne, name = "B"),
    M("Diag", ne, ne, free = FALSE, values = 1, name = "C",
      dimnames = list(paste0("X", 1:ne), paste0("F", 1:ne))),
    M("Zero", ne, ne, name = "D"),
    M("Symm", ne, ne, free = FALSE, values = compute_dtQ_identity(A_ct, deltat),
      name = "Q", lbound = 0),
    M("Diag", ne, ne, free = FALSE, values = 1e-5, name = "R"),
    M("Full", ne, 1, free = FALSE, values = 0, name = "x0"),
    M("Diag", ne, ne, free = FALSE, values = 1, name = "P0"),
    M("Zero", ne, 1, name = "u"),
    OpenMx::mxExpectationStateSpace("A", "B", "C", "D", "Q", "R", "x0", "P0", "u"))
  dat <- OpenMx::mxGenerateData(osc, nrows = times + 1000L)
  dat <- dat[1001:(times + 1000L), ]
  dat$Time <- (0:(times - 1L)) * deltat
  dat$id <- as.integer(xp)
  if (retain_only) {
    # Apply the original floating-point modulo rule BEFORE discarding any rows.
    keep <- sort(unique(unlist(lapply(c(.5, 1, 5), function(dt)
      head(which(abs(dat$Time %% dt) < 1e-8), 200L)))))
    dat <- dat[keep, ]
  }
  dat
}

generate_one <- function(task, truth, output_dir, overwrite = FALSE, retain_only = TRUE) {
  effect <- task$effect
  sim <- task$sim
  folder <- file.path(output_dir, paste0("Simulations_", effect))
  path <- file.path(folder, paste0("Dataset ", sim, ".RDS"))
  true_path <- file.path(folder, paste0("Trues ", sim, ".RDS"))
  if (!overwrite && file.exists(path) && file.exists(true_path)) return(path)
  xx <- truth[[effect]][[sim]]
  parts <- lapply(1:60, function(xp)
    simulate_subject(sim, xp, xx[[ceiling(xp / 20)]], retain_only = retain_only))
  dataset <- do.call(rbind, parts)
  rownames(dataset) <- NULL
  saveRDS(xx, true_path)
  saveRDS(dataset, path, compress = FALSE)
  message(effect, " replication ", sim, " complete")
  path
}

generate_data <- function(opt = read_options()) {
  for (pkg in c("OpenMx", "expm")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Install package '", pkg, "'.")
  }
  OpenMx::mxOption(NULL, "Number of Threads", 1)
  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
  truth <- setNames(vector("list", length(opt$effects)), opt$effects)
  for (effect in opt$effects) {
    truth[[effect]] <- if (opt$new_truth) {
      # Independent effect seeds make subset runs agree with a full new study.
      RNGkind("Mersenne-Twister", "Inversion", "Rejection")
      set.seed(opt$truth_seed + match(effect, c("Weak", "Mid", "Large")) - 1L)
      lapply(seq_len(max(opt$replications)), function(i) create.amat(6, .10,
        allow.neg = TRUE, EF = .50, EF.sd = .01,
        SEF = c(Weak = .30, Mid = .60, Large = .90)[[effect]],
        subgroups = 3, sub.prop = .10))
    } else readRDS(file.path("data", "truth", paste0(effect, ".rds")))
    dir.create(file.path(opt$output_dir, paste0("Simulations_", effect)), showWarnings = FALSE)
  }
  grid <- expand.grid(effect = opt$effects, sim = opt$replications, stringsAsFactors = FALSE)
  tasks <- split(grid, seq_len(nrow(grid)))
  workers <- min(opt$workers, length(tasks))
  if (workers == 1L) {
    paths <- lapply(tasks, generate_one, truth = truth, output_dir = opt$output_dir,
                    overwrite = opt$overwrite, retain_only = opt$retain_only)
  } else {
    cl <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, OpenMx::mxOption(NULL, "Number of Threads", 1))
    parallel::clusterExport(cl, c("simulate_subject", "compute_dtQ_identity", "generate_one"),
                            envir = environment(generate_one))
    paths <- parallel::parLapplyLB(cl, tasks, generate_one, truth = truth,
                                  output_dir = opt$output_dir, overwrite = opt$overwrite,
                                  retain_only = opt$retain_only)
  }
  message("Datasets written to: ", normalizePath(opt$output_dir, winslash = "/"))
  invisible(paths)
}

if (sys.nframe() == 0L) generate_data()

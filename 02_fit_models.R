# Fit the paper's S-GIMME and continuous-time PAM models.
# Run from this folder: Rscript 02_fit_models.R
# First install ctgimme from CRAN and the other dependencies (see README).
# CONTROLS below select effect conditions, intervals, replications, data/output
# directories and CPU workers. Example, one replication at all three intervals:
# Rscript 02_fit_models.R --effects=Weak --replications=1 --cores=4
# --data-dir=generated --output-dir=results can point to other locations.
# Inputs: 01's Simulations_<effect>/Dataset <r>.RDS and Trues <r>.RDS.
# Outputs: original fitted models, per-job metrics.rds, memberships.csv, logs,
# and combined memberships.csv for 03_figures.R --memberships=results/memberships.csv.
# The combined file includes every completed job in output-dir across invocations.
# The settings below are the archived study settings; changing them changes the
# experiment. Each model family runs in a fresh R process, as in the study.
# Completed jobs resume only when inputs, settings and software versions match.
# Incomplete jobs stop with their log path; use a NEW output directory to retry.
# Original fitting seeds were not recorded. seed=NULL preserves that behavior;
# set an integer to seed NEW fitting processes (not historical seed recovery).
# Parallel OpenMx restart draws can still vary because workers have separate RNGs.

effects <- c("Weak", "Mid", "Large")
replications <- 1:500
deltas <- c(0.5, 1, 5)
data_dir <- "generated"
output_dir <- "results"
cores <- 4L
seed <- NULL
settings <- list(sig.thrsh = 0.55, sub.sig.thrsh = 0.55, Galpha = 0.05,
                 S.Galpha = 0.05, Ialpha = 0.01, ben.hoch = TRUE,
                 subgroup.method = "pam", max.subgroups = 6,
                 ME.var = diag(1e-5, 6), PE.var = diag(1, 6))

Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           BLIS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  hit <- args[startsWith(args, paste0("--", name, "="))]
  if (length(hit)) sub(paste0("^--", name, "="), "", tail(hit, 1)) else default
}
parse_ids <- function(x) {
  if (!grepl("^[0-9]+(:[0-9]+)?(,[0-9]+(:[0-9]+)?)*$", x))
    stop("Use --replications=1, 1:500, or a comma-separated combination.")
  unlist(lapply(strsplit(x, ",", fixed = TRUE)[[1]], function(part) {
    bounds <- as.integer(strsplit(part, ":", fixed = TRUE)[[1]])
    if (length(bounds) == 2L) seq.int(bounds[1], bounds[2]) else bounds
  }), use.names = FALSE)
}
effects <- strsplit(option("effects", paste(effects, collapse = ",")), ",", fixed = TRUE)[[1]]
replications <- unique(parse_ids(option("replications", paste(replications, collapse = ","))))
deltas <- as.numeric(strsplit(option("deltas", paste(deltas, collapse = ",")), ",", fixed = TRUE)[[1]])
data_dir <- option("data-dir", data_dir)
output_dir <- option("output-dir", output_dir)
cores <- as.integer(option("cores", cores))
seed_arg <- option("seed", "")
if (nzchar(seed_arg)) seed <- as.integer(seed_arg)
stopifnot(all(effects %in% c("Weak", "Mid", "Large")), length(effects) > 0,
          all(replications %in% 1:500), length(replications) > 0,
          all(deltas %in% c(0.5, 1, 5)), length(deltas) > 0,
          length(cores) == 1L, !is.na(cores), cores >= 1L)
required <- c("OpenMx", "expm", "gimme", "ctgimme", "mclust")
for (p in required) if (!requireNamespace(p, quietly = TRUE)) stop("Install ", p, "; see README.")
versions <- vapply(required, function(p) as.character(utils::packageVersion(p)), "")
if (versions[["gimme"]] != "0.7.18")
  stop("Use gimme 0.7.18; see README for installation.")
source("R/metrics.R")
suppressPackageStartupMessages(library(OpenMx))
OpenMx::mxOption(NULL, "Number of Threads", 1)

subsample <- function(data, dt) {
  data <- data[order(data$id, data$Time), ]
  data <- data[abs(data$Time %% dt) < 1e-8, ]
  data <- do.call(rbind, lapply(split(data, data$id), function(x) head(x, 200)))
  rownames(data) <- NULL
  data
}
align_membership <- function(labels, ids) {
  ids <- as.integer(ids)
  if (length(ids) != 60L || length(labels) != 60L || anyNA(ids) ||
      !identical(sort(ids), 1:60) || anyNA(labels))
    stop("Membership must contain exactly one nonmissing label for every subject 1:60.")
  labels[order(ids)]
}

# Internal child mode: separate processes release all fitting resources between jobs.
job_file <- option("job", "")
if (nzchar(job_file)) {
  job <- readRDS(job_file)
  phase <- option("phase", "")
  if (!is.null(job$seed)) set.seed(job$seed)
  data <- subsample(readRDS(job$dataset), job$delta)
  stopifnot(identical(sort(unique(data$id)), as.numeric(1:60)) ||
              identical(sort(unique(data$id)), 1:60), all(table(data$id) == 200L),
            all(paste0("X", 1:6) %in% names(data)))
  if (phase == "gimme") {
    suppressPackageStartupMessages(library(gimme))
    install_gimme_empty_subgroup_hotfix()
    series <- lapply(1:60, function(i) data[data$id == i, paste0("X", 1:6), drop = FALSE])
    gimme::gimme(data = series, out = file.path(job$directory, "gimme"), subgroup = TRUE)
  } else if (phase == "ctgimme") {
    suppressPackageStartupMessages({library(expm); library(gimme); library(ctgimme)})
    fit <- do.call(ctgimme::ctgimme, c(list(varnames = paste0("X", 1:6), dataframe = data,
                    id = "id", time = "Time", cores = job$cores, directory = job$directory), job$settings))
    membership <- attr(fit, "ctgimme.membership")
    if (is.null(membership)) membership <- fit$membership
    membership <- align_membership(membership, names(membership))
    saveRDS(membership, file.path(job$directory, "ct_membership.rds"))
  } else stop("Unknown internal phase.")
  quit(save = "no")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
run_config <- list(settings = settings, cores = cores, versions = versions,
                   R = as.character(getRversion()), seed = seed)
config_path <- file.path(output_dir, "run-config.rds")
# Only traverse condition and replication folders, not the many model files.
condition_dirs <- list.dirs(output_dir, recursive = FALSE, full.names = TRUE)
job_dirs <- unlist(lapply(condition_dirs, list.dirs, recursive = FALSE, full.names = TRUE),
                   use.names = FALSE)
markers <- file.path(job_dirs, "completed.rds")
markers <- markers[file.exists(markers)]
if (file.exists(config_path)) {
  if (!identical(readRDS(config_path), run_config))
    stop("Changed fitting settings/software: choose a new output-dir.")
} else {
  # Also validate completed jobs from an earlier invocation lacking the root file.
  for (path in markers) if (!identical(readRDS(path)[names(run_config)], run_config))
    stop("Existing jobs use different fitting settings/software; choose a new output-dir.")
  saveRDS(run_config, config_path)
}
script <- normalizePath("02_fit_models.R", winslash = "/", mustWork = TRUE)
script_exe <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
all_memberships <- setNames(lapply(markers, function(path) {
  if (!identical(readRDS(path)[names(run_config)], run_config))
    stop("Completed job has different settings/software: ", dirname(path))
  read.csv(file.path(dirname(path), "memberships.csv"))
}), normalizePath(dirname(markers), winslash = "/", mustWork = TRUE))
write_memberships <- function() {
  rows <- if (length(all_memberships)) do.call(rbind, all_memberships) else
    data.frame(signal = character(), effect = numeric(), delta = numeric(),
               replication = integer(), subject = integer(), truth = integer(),
               ct_group = integer(), dt_group = integer())
  write.csv(rows, file.path(output_dir, "memberships.csv"), row.names = FALSE)
}
write_memberships()
effect_values <- c(Weak = 0.3, Mid = 0.6, Large = 0.9)
delta_codes <- c(`0.5` = "0_5", `1` = "1_0", `5` = "5_0")
for (effect in effects) for (replication in replications) {
  dataset <- normalizePath(file.path(data_dir, paste0("Simulations_", effect),
                                     paste0("Dataset ", replication, ".RDS")), winslash = "/", mustWork = TRUE)
  truth_path <- normalizePath(file.path(dirname(dataset), paste0("Trues ", replication, ".RDS")),
                             winslash = "/", mustWork = TRUE)
  hashes <- unname(tools::md5sum(c(dataset, truth_path)))
  for (delta in deltas) {
    condition <- paste(effect, delta_codes[[as.character(delta)]], sep = "_")
    directory <- file.path(output_dir, condition, paste0("ctsgtest", replication))
    signature <- list(inputs = hashes, delta = delta, settings = settings, cores = cores,
                      versions = versions, R = as.character(getRversion()), seed = seed)
    complete <- file.path(directory, "completed.rds")
    if (file.exists(complete)) {
      if (!identical(readRDS(complete), signature)) stop("Changed inputs/settings: choose a new output-dir.")
      rows <- read.csv(file.path(directory, "memberships.csv"))
    } else {
      if (dir.exists(directory) && length(list.files(directory, all.files = TRUE, no.. = TRUE)))
        stop("Incomplete job at ", directory, ". Review its logs; retry in a new output-dir.")
      dir.create(directory, recursive = TRUE, showWarnings = FALSE)
      directory <- normalizePath(directory, winslash = "/", mustWork = TRUE)
      message(condition, ", replication ", replication)
      job_file <- file.path(directory, "job.rds")
      job_seed <- if (is.null(seed)) NULL else seed + match(effect, names(effect_values)) * 10000L +
        replication * 10L + match(delta, c(0.5, 1, 5))
      saveRDS(list(dataset = dataset, delta = delta, directory = directory, settings = settings,
                   cores = cores, seed = job_seed), job_file)
      for (phase in c("gimme", "ctgimme")) {
        log <- file.path(directory, paste0(phase, ".log"))
        status <- system2(script_exe, c(shQuote(script), shQuote(paste0("--job=", job_file)),
                                        paste0("--phase=", phase)), stdout = log, stderr = log)
        if (status != 0) stop(phase, " failed; see ", log)
      }
      ct <- readRDS(file.path(directory, "ct_membership.rds"))
      dt <- read.csv(file.path(directory, "gimme", "summaryFit.csv"))
      # Match the archived metric helper: use the last integer in each filename.
      id_parts <- regmatches(dt$file, gregexpr("[0-9]+", dt$file, perl = TRUE))
      dt_ids <- vapply(id_parts, function(x) if (length(x)) as.integer(tail(x, 1)) else NA_integer_, 1L)
      dt_groups <- align_membership(dt$sub_membership, dt_ids)
      rows <- data.frame(signal = effect, effect = unname(effect_values[effect]), delta = delta,
                         replication = replication, subject = 1:60, truth = rep(1:3, each = 20),
                         ct_group = ct, dt_group = dt_groups)
      metrics <- compute_simulation_metrics(
        models_dir = file.path(directory, "Models", "Individuals"),
        gstruc_path = file.path(directory, "GStruc.RDS"), sgstruc_dir = directory,
        subjects = 60, comp_mats = readRDS(truth_path), membership = ct,
        subjects_per_group = 20, fit_csv_path = file.path(directory, "gimme", "summaryFit.csv"),
        sig_level = 0.05, digits = 2, bias_aggregation = "pooled")
      saveRDS(metrics, file.path(directory, "metrics.rds"))
      write.csv(rows, file.path(directory, "memberships.csv"), row.names = FALSE)
      saveRDS(signature, complete)
    }
    all_memberships[[directory]] <- rows
    write_memberships()
  }
}
write_memberships()
message("Finished selected jobs. Figure input: ", file.path(output_dir, "memberships.csv"))

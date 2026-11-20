# =============================================================================
# Heatmaps or Text? Sample-size sensitivity analysis (Section 4.2)
# Simulates the wrong-advice cases where the first answer was correct and
# estimates power for the primary test: text vs heatmap (alpha = .05).
# Package: lme4. Runtime: a few minutes.
# =============================================================================
suppressMessages(library(lme4))
library(parallel)
CORES <- if (.Platform$OS.type == "windows") 1 else detectCores()

FORMATS  <- c("none", "heatmap", "text", "both")
WILLIAMS <- rbind(c(1, 2, 4, 3),   # rows: the four orders; columns: blocks 1-4
                  c(2, 3, 1, 4),
                  c(3, 4, 2, 1),
                  c(4, 1, 3, 2))

sim_once <- function(N, OR, p_first_correct = 0.75, p_switch_heatmap = 0.30,
                     sd_part = 0.8, sd_case = 1.0) {
  order <- as.vector(replicate(N / 4, sample(4)))
  d <- expand.grid(case = 1:16, part = seq_len(N))
  d$block  <- (d$case - 1) %/% 4 + 1                 # 4 wrong-advice cases per block
  d$format <- factor(FORMATS[WILLIAMS[cbind(order[d$part], d$block)]], levels = FORMATS)
  d <- d[rbinom(nrow(d), 1, p_first_correct) == 1, ]

  eff <- c(none = 0, heatmap = 0, text = log(OR), both = log(OR) / 2)
  u <- rnorm(N, 0, sd_part)
  v <- rnorm(16, 0, sd_case)
  eta <- qlogis(p_switch_heatmap) + eff[as.character(d$format)] + u[d$part] + v[d$case]
  d$switched <- rbinom(nrow(d), 1, plogis(eta))
  d$block_c  <- d$block - 2.5

  fit <- tryCatch(suppressMessages(suppressWarnings(
    glmer(switched ~ format + block_c + (1 | part) + (1 | case), data = d,
          family = binomial, control = glmerControl(optimizer = "bobyqa")))),
    error = function(e) NULL)
  if (is.null(fit)) return(NA)

  b <- fixef(fit)
  V <- as.matrix(vcov(fit))
  est <- b["formattext"] - b["formatheatmap"]
  se  <- sqrt(V["formattext", "formattext"] + V["formatheatmap", "formatheatmap"] -
                2 * V["formattext", "formatheatmap"])
  as.numeric(2 * pnorm(-abs(est / se)) < 0.05)
}

n_runs <- 500
grid <- expand.grid(N = c(24, 32), OR = c(2.5, 3.0, 3.5))
for (i in seq_len(nrow(grid))) {
  hits <- unlist(mclapply(seq_len(n_runs), function(k) {
    set.seed(10000 * i + k)                          # same results on any machine
    sim_once(grid$N[i], grid$OR[i])
  }, mc.cores = CORES))
  ok <- hits[!is.na(hits)]
  pw <- mean(ok)
  se <- sqrt(pw * (1 - pw) / length(ok))
  grid$power[i]       <- round(pw, 3)
  grid$ci_low[i]      <- round(pw - 1.96 * se, 3)
  grid$ci_high[i]     <- round(pw + 1.96 * se, 3)
  grid$failed_fits[i] <- n_runs - length(ok)
  grid$text_rate[i]   <- round(plogis(qlogis(0.30) + log(grid$OR[i])), 2)
}
print(grid)
write.csv(grid, "power_results.csv", row.names = FALSE)

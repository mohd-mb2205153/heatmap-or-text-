# =============================================================================
# Heatmaps or Text? Dummy data and planned analysis (Section 5)
# Part A simulates dummy data with the study design; the values are invented.
# Part B runs the planned analysis on the three CSV files, so it can run
# unchanged on the real data exported in the same format.
# Packages: lme4, emmeans, ordinal, performance, DHARMa
# =============================================================================
suppressPackageStartupMessages({
  library(lme4)
  library(emmeans)
  library(ordinal)
})
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)
dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

FORMATS  <- c("none", "heatmap", "text", "both")          # both = heatmap plus text
WILLIAMS <- rbind(c(1, 2, 4, 3), c(2, 3, 1, 4),            # rows: orders 1-4
                  c(3, 4, 2, 1), c(4, 1, 3, 2))            # columns: blocks 1-4

# ==== PART A. Simulate dummy data ============================================
simulate_study <- function(N = 24, n_careless = 0, seed = 2026) {
  stopifnot(N %% 4 == 0)
  set.seed(seed)
  likert <- function(m) as.integer(pmin(5, pmax(1, round(m))))

  # Cases: set k is always shown in block k. In each set, cases 1-8 are
  # metastases; the AI is wrong on 1-2 and 9-10; case 16 is the attention check.
  pos <- rep(1:16, 4)
  cases <- data.frame(case = sprintf("C%02d", 1:64), block = rep(1:4, each = 16))
  cases$truth      <- as.integer(pos <= 8)
  cases$ai_correct <- as.integer(!pos %in% c(1, 2, 9, 10))
  cases$ai_says    <- ifelse(cases$ai_correct == 1, cases$truth, 1L - cases$truth)
  cases$attention  <- as.integer(pos == 16)
  cases$difficulty <- rnorm(64) + 0.5 * (cases$ai_correct == 0)
  cases$difficulty[cases$attention == 1] <- -4
  case_re <- rnorm(64, 0, 0.6)

  # Participants. Each group of four gets the four orders in random sequence.
  p <- data.frame(part = sprintf("P%02d", 1:N))
  p$order <- as.vector(replicate(N / 4, sample(4)))
  n_hmc <- round(0.75 * N)
  p$site <- sample(rep(c("HMC", "Sidra"), c(n_hmc, N - n_hmc)))
  p$role <- ifelse(p$site == "HMC",
                   sample(c("consultant", "resident"), N, TRUE, prob = c(0.55, 0.45)),
                   "consultant")
  p$experience <- ifelse(p$role == "resident", sample(1:4, N, TRUE),
                         sample(5:25, N, TRUE))
  p$subspecialty <- ifelse(p$site == "Sidra",
    sample(c("breast_general", "other"), N, TRUE, prob = c(0.2, 0.8)),
    sample(c("breast_general", "other"), N, TRUE, prob = c(0.7, 0.3)))
  p$ln_frequency <- ifelse(p$subspecialty == "breast_general",
    sample(c("weekly", "less_often"), N, TRUE, prob = c(0.8, 0.2)),
    sample(c("weekly", "less_often"), N, TRUE, prob = c(0.3, 0.7)))
  p$digital <- sample(c("yes", "no"), N, TRUE)
  p$training <- ifelse(p$role == "resident", "Arab region",
    sample(c("Arab region", "South Asia", "Europe/N. America"), N, TRUE,
           prob = c(0.3, 0.4, 0.3)))
  p$first_language <- ifelse(p$training == "Europe/N. America",
    sample(c("English", "other"), N, TRUE, prob = c(0.6, 0.4)),
    ifelse(p$training == "Arab region", "Arabic", "other"))

  # Questionnaire answers, scored in Part B as described in Appendix A.2
  prof <- rnorm(N) + 0.8 * (p$first_language == "English")
  p$lex_words    <- rbinom(N, 40, plogis(1.5 + 0.8 * prof))   # of 40 words
  p$lex_nonwords <- rbinom(N, 20, plogis(1.2 + 0.6 * prof))   # of 20 nonwords
  tia <- rnorm(N, 0, 0.7)
  p$tia_5  <- likert(3.8 - tia + rnorm(N, 0, 0.6))            # reverse-keyed item
  p$tia_12 <- likert(3.1 + tia + rnorm(N, 0, 0.6))
  p$tia_18 <- likert(3.4 + tia + rnorm(N, 0, 0.6))
  lex_c    <- (100 * (p$lex_words / 40 + p$lex_nonwords / 20) / 2 - 78) / 10
  skill    <- rnorm(N, 0, 0.6) + 0.03 * (p$experience - 10)
  reliance <- rnorm(N, 0, 0.8)

  # Trials
  t <- do.call(rbind, lapply(seq_len(N), function(i) {
    d <- cases
    d$part <- p$part[i]
    d$format <- FORMATS[WILLIAMS[p$order[i], d$block]]
    d$order_in_block <- ave(seq_len(64), d$block,
                            FUN = function(x) sample(seq_along(x)))
    d
  }))
  t <- t[order(t$part, t$block, t$order_in_block), ]
  ip <- match(t$part, p$part)
  ic <- match(t$case, cases$case)

  init_ok <- rbinom(nrow(t), 1, plogis(1.4 + skill[ip] - t$difficulty))
  t$initial <- ifelse(init_ok == 1, t$truth, 1L - t$truth)
  t$conf_initial <- as.integer(cut(0.9 * init_ok - 0.4 * t$difficulty + rnorm(nrow(t)),
                                   c(-Inf, -1.2, -0.3, 0.5, 1.4, Inf)))

  # Effects on the log-odds scale (dummy values)
  b_acc <- c(none = 0, heatmap = 0.4, text = 0.5, both = 0.7)          # accept correct
  b_fol <- c(none = 0, heatmap = 0, text = log(3), both = log(3) / 2)  # follow wrong
  eta <- ifelse(t$ai_correct == 1,
                qlogis(0.50) + b_acc[t$format],
                qlogis(0.25) + b_fol[t$format] -
                  0.3 * lex_c[ip] * (t$format %in% c("text", "both")))
  eta <- eta + reliance[ip] + case_re[ic] - 0.4 * (t$conf_initial - 3)
  switch_to_ai <- rbinom(nrow(t), 1, plogis(eta))
  drift <- rbinom(nrow(t), 1, 0.02)
  disagree <- t$initial != t$ai_says
  t$final <- ifelse(disagree, ifelse(switch_to_ai == 1, t$ai_says, t$initial),
                    ifelse(drift == 1, 1L - t$initial, t$initial))

  if (n_careless > 0) {                    # only used to test the exclusion rules
    k <- t$part %in% tail(p$part, n_careless)
    t$initial[k] <- rbinom(sum(k), 1, 0.5)
    t$final[k] <- t$initial[k]
    ka <- k & t$attention == 1
    t$final[ka] <- 1L - t$truth[ka]
  }

  agree_final <- t$final == t$ai_says
  t$conf_final <- pmin(5L, pmax(1L, t$conf_initial +
    round(0.6 * agree_final - 0.6 * (t$final != t$initial) + rnorm(nrow(t), 0, 0.6))))
  t$time_s <- round(exp(log(30) + 0.15 * (t$format %in% c("text", "both")) +
    0.1 * t$difficulty + rnorm(N, 0, 0.3)[ip] + rnorm(nrow(t), 0, 0.4)), 1)
  t$time_s[sample(nrow(t), 6)] <- 2.0
  t$time_s[sample(nrow(t), 3)] <- 420

  # Trust scale after each block (8 items, 1-5; item 6 is negatively worded)
  tr <- unique(t[, c("part", "block", "format")])
  lat <- 0.7 * tia[match(tr$part, p$part)] +
    0.3 * (tr$format %in% c("text", "both")) + rnorm(nrow(tr), 0, 0.5)
  for (j in 1:8) {
    x <- likert(3.3 + lat + rnorm(nrow(tr), 0, 0.6))
    tr[[paste0("trust", j)]] <- if (j == 6) 6L - x else x
  }

  trials <- t[, c("part", "case", "block", "order_in_block", "format", "truth",
                  "ai_says", "ai_correct", "attention", "initial", "conf_initial",
                  "final", "conf_final", "time_s")]
  write.csv(trials, "data/trials.csv", row.names = FALSE)
  write.csv(p, "data/participants.csv", row.names = FALSE)
  write.csv(tr, "data/trust.csv", row.names = FALSE)
}

N_SIM <- 24; N_CARELESS <- 0              # N must be a multiple of four
simulate_study(N = N_SIM, n_careless = N_CARELESS, seed = 2026)

# ==== PART B. Planned analysis (dummy or real data) ===========================
trials       <- read.csv("data/trials.csv", stringsAsFactors = FALSE)
participants <- read.csv("data/participants.csv", stringsAsFactors = FALSE)
trust        <- read.csv("data/trust.csv", stringsAsFactors = FALSE)

# ---- B1. Questionnaire scores (Appendix A.2) --------------------------------
participants$lextale <- 100 * (participants$lex_words / 40 +
                                 participants$lex_nonwords / 20) / 2
participants$propensity <- rowMeans(cbind(6 - participants$tia_5, participants$tia_12,
                                          participants$tia_18), na.rm = TRUE)
items <- trust[, paste0("trust", 1:8)]
items$trust6 <- 6 - items$trust6
trust$score <- rowMeans(items)

# ---- B2. Design checks ------------------------------------------------------
check <- function(ok, msg) if (!isTRUE(ok)) warning(msg, call. = FALSE)
n_p <- length(unique(trials$part))
blk <- unique(trials[, c("part", "block", "format")])
key <- paste(blk$part, blk$block)
blk$prev_format <- ifelse(blk$block == 1, "first",
                          blk$format[match(paste(blk$part, blk$block - 1), key)])
check(!any(duplicated(trials[, c("part", "case")])), "A participant saw a case twice")
check(all(tapply(trials$ai_says, trials$case, function(x) length(unique(x)) == 1)),
      "AI advice differs across formats for some case")
check(all(tapply(trials$block, trials$case, function(x) length(unique(x)) == 1)),
      "A case set is not fixed to one block")
check(all(table(trials$case, trials$format) == n_p / 4),
      "Cases are not shown equally often under each format")
check(all(table(blk$format, blk$block) == n_p / 4),
      "Formats are not balanced across block positions")
nxt <- table(factor(blk$prev_format, FORMATS), factor(blk$format, FORMATS))
check(all(nxt[row(nxt) != col(nxt)] == n_p / 4),
      "Formats do not follow each other equally often")

# ---- B3. Pre-registered exclusions ------------------------------------------
blocks_done <- tapply(trials$block, trials$part, function(x) length(unique(x)))
att <- trials[trials$attention == 1, ]
att_fail <- tapply(att$final != att$truth, att$part, sum)
drop_parts <- union(names(blocks_done)[blocks_done < 2], names(att_fail)[att_fail >= 2])
d <- trials[!(trials$part %in% drop_parts) & trials$attention == 0, ]
bad_time <- d$time_s < 3 | d$time_s > 300
print(c("participants: < 2 blocks" = sum(blocks_done < 2),
        "participants: >= 2 failed attention checks" = sum(att_fail >= 2),
        "trials: < 3 s or > 5 min" = sum(bad_time)))
d <- d[!bad_time, ]

# ---- B4. Outcomes (reliance is only visible when first answer != advice) ----
d$disagree <- d$initial != d$ai_says
d$accept_correct <- ifelse(d$ai_correct == 1 & d$disagree,
                           as.integer(d$final == d$ai_says), NA)
d$reject_incorrect <- ifelse(d$ai_correct == 0 & d$disagree,
                             as.integer(d$final == d$initial), NA)

d <- merge(d, participants, by = "part", sort = FALSE)
d$format <- factor(d$format, levels = FORMATS)
d$prev_format <- factor(blk$prev_format[match(paste(d$part, d$block), key)],
                        levels = c("first", FORMATS))
d$experience_c <- d$experience - mean(participants$experience)
d$lextale_c <- (d$lextale - mean(participants$lextale)) / 10      # per 10 points
d$block_c <- d$block - 2.5
for (v in c("site", "subspecialty", "ln_frequency", "digital", "training"))
  d[[v]] <- factor(d[[v]])

# ---- B5. Descriptives (Table 2) ---------------------------------------------
describe <- function(y) {
  k <- !is.na(d[[y]])
  n <- tapply(d[[y]][k], d$format[k], length)
  s <- tapply(d[[y]][k], d$format[k], sum)
  data.frame(format = FORMATS, n = as.vector(n), events = as.vector(s),
             pct = round(100 * as.vector(s / n), 1))
}
desc <- merge(describe("accept_correct"), describe("reject_incorrect"),
              by = "format", suffixes = c("_accept", "_reject"), sort = FALSE)
print(desc)
write.csv(desc, "results/descriptives.csv", row.names = FALSE)

# ---- B6. Primary models (RQ1) -----------------------------------------------
PRIMARY  <- "format + block_c"
CONTROLS <- "experience_c + subspecialty + ln_frequency + digital + site"
CTRL <- quote(glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

fit_glmer <- function(f, data_name = "d")
  suppressWarnings(suppressMessages(eval(bquote(
    glmer(.(f), data = .(as.name(data_name)), family = binomial, control = .(CTRL))))))

converged <- function(m)
  length(c(m@optinfo$conv$lme4$messages, m@optinfo$warnings)) == 0

check_fit <- function(m, label) {
  sep <- names(fixef(m))[abs(fixef(m)) > 10]
  cat(sprintf("[%s] converged: %s | singular: %s | possible separation: %s\n",
              label, converged(m), isSingular(m),
              if (length(sep)) paste(sep, collapse = ", ") else "none"))
  if (!converged(m)) {                     # refit with all optimisers
    fx <- summary(suppressWarnings(suppressMessages(allFit(m, verbose = FALSE))))$fixef
    spread <- apply(fx[, grep("^format", colnames(fx)), drop = FALSE], 2,
                    function(x) diff(range(x)))
    cat("  allFit: largest spread of format estimates =", signif(max(spread), 2), "\n")
  }
  invisible(converged(m) && !isSingular(m))
}

# Random slopes for format first (Barr et al. 2013), else random intercepts
fit_primary <- function(y, rhs = PRIMARY, data_name = "d") {
  m <- fit_glmer(as.formula(paste(y, "~", rhs, "+ (1 + format | part) + (1 | case)")),
                 data_name)
  re <- "random slopes"
  if (!(converged(m) && !isSingular(m))) {
    m <- fit_glmer(as.formula(paste(y, "~", rhs, "+ (1 | part) + (1 | case)")),
                   data_name)
    re <- "random intercepts"
  }
  attr(m, "re_structure") <- re
  m
}

# Primary test: text vs heatmap for rejecting incorrect advice (alpha = .05).
# All other contrasts are Holm-corrected within each outcome.
CONTRASTS <- list("heatmap vs none" = c(-1, 1, 0, 0),
                  "text vs none"    = c(-1, 0, 1, 0),
                  "both vs none"    = c(-1, 0, 0, 1),
                  "text vs heatmap" = c(0, -1, 1, 0),
                  "both vs text"    = c(0, 0, -1, 1))

planned_contrasts <- function(m, outcome) {
  ct <- contrast(emmeans(m, ~ format), method = CONTRASTS, adjust = "none")
  ct <- as.data.frame(summary(ct, infer = TRUE, type = "response"))
  primary <- outcome == "reject_incorrect" & ct$contrast == "text vs heatmap"
  p_adj <- ct$p.value
  p_adj[!primary] <- p.adjust(ct$p.value[!primary], method = "holm")
  data.frame(contrast = ct$contrast, role = ifelse(primary, "primary", "Holm"),
             OR = round(ct$odds.ratio, 2), CI_low = round(ct$asymp.LCL, 2),
             CI_high = round(ct$asymp.UCL, 2), p_raw = signif(ct$p.value, 3),
             p_adj = signif(p_adj, 3))
}

predicted_pct <- function(m) {
  pr <- as.data.frame(summary(emmeans(m, ~ format), type = "response"))
  data.frame(format = pr$format, pct = round(100 * pr$prob, 1),
             low = round(100 * pr$asymp.LCL, 1), high = round(100 * pr$asymp.UCL, 1))
}

results <- list()
for (y in c("accept_correct", "reject_incorrect")) {
  m <- fit_primary(y)
  cat("\n=====", y, "| n trials =", nobs(m), "|", attr(m, "re_structure"), "\n")
  check_fit(m, paste("primary", y))
  results[[y]] <- list(model = m, contrasts = planned_contrasts(m, y),
                       predicted = predicted_pct(m))
  print(results[[y]]$contrasts)
  print(results[[y]]$predicted)
  write.csv(results[[y]]$contrasts, paste0("results/contrasts_", y, ".csv"),
            row.names = FALSE)
}

# ---- B7. Model checks ---------------------------------------------------------
# Is the spread of rates between participants and between cases typical of
# data simulated from the fitted model?
sim_check <- function(m, nsim = 500) {
  fr <- model.frame(m)
  y <- model.response(fr)
  sims <- simulate(m, nsim = nsim, seed = 1)
  do.call(rbind, lapply(c("part", "case"), function(g) {
    obs <- var(tapply(y, fr[[g]], mean))
    simv <- sapply(sims, function(s) var(tapply(s, fr[[g]], mean)))
    data.frame(between = g, observed = round(obs, 3), simulated = round(mean(simv), 3),
               p = round(min(1, 2 * min(mean(simv >= obs), mean(simv <= obs))), 3))
  }))
}
for (y in names(results)) {
  m <- results[[y]]$model
  cat("\n--- Checks:", y, "---\n")
  if (has_pkg("DHARMa")) tryCatch({
    res <- DHARMa::simulateResiduals(m, n = 1000, seed = 1)
    cat("Uniformity p =",
        signif(DHARMa::testUniformity(res, plot = FALSE)$p.value, 3), "\n")
    res_p <- DHARMa::recalculateResiduals(res, group = m@frame$part)
    cat("Dispersion by participant p =",
        signif(DHARMa::testDispersion(res_p, plot = FALSE)$p.value, 3), "\n")
  }, error = function(e) print(e))
  print(sim_check(m))
}

# ---- B8. RQ2 (exploratory): training background and English proficiency ----
report_rq2 <- function(m1) {
  if ("lextale_c" %in% all.vars(formula(m1)))
    print(summary(emtrends(m1, ~ format, var = "lextale_c"), infer = TRUE))
  if ("training" %in% all.vars(formula(m1)))
    print(summary(contrast(emmeans(m1, ~ format | training),
                           method = CONTRASTS["text vs heatmap"]),
                  infer = TRUE, type = "response"))
}
for (y in names(results)) {
  m0 <- results[[y]]$model
  cat("\n=====", y, ": RQ2 moderation\n")
  m1 <- fit_glmer(update(formula(m0), . ~ . + format * training + format * lextale_c))
  if (check_fit(m1, "joint moderation")) {
    print(anova(m0, m1))
    report_rq2(m1)
  } else for (mod in c("training", "lextale_c")) {     # one moderator at a time
    m1s <- fit_glmer(update(formula(m0), as.formula(paste(". ~ . + format *", mod))))
    check_fit(m1s, paste("moderator", mod))
    print(anova(m0, m1s))
    report_rq2(m1s)
  }
}

# ---- B9. Secondary outcomes -------------------------------------------------
d$conf_final_o <- factor(d$conf_final, levels = 1:5, ordered = TRUE)
d$conf_initial_c <- d$conf_initial - 3
d$ai_correct_f <- factor(d$ai_correct, levels = c(1, 0),
                         labels = c("correct", "incorrect"))
m_conf <- clmm(conf_final_o ~ format * ai_correct_f + conf_initial_c +
                 (1 | part) + (1 | case), data = d)
cat("\n===== Final confidence (ordinal mixed model)\n")
print(summary(m_conf)$coefficients)

m_time <- lmer(log(time_s) ~ format * ai_correct_f + block_c + (1 | part) + (1 | case),
               data = d)
cat("\n===== Decision time ratios vs no explanation\n")
print(round(exp(fixef(m_time)[grep("^format[a-z]+$", names(fixef(m_time)))]), 2))

cronbach <- function(x) {
  k <- ncol(x)
  k / (k - 1) * (1 - sum(apply(x, 2, var)) / var(rowSums(x)))
}
cat("\nCronbach's alpha, trust scale:", round(cronbach(items), 2), "\n")
trust <- merge(trust[!(trust$part %in% drop_parts), ],
               participants[, c("part", "propensity")], by = "part")
trust$format <- factor(trust$format, levels = FORMATS)
trust$propensity_c <- trust$propensity - mean(participants$propensity)
m_trust <- lmer(score ~ format + I(block - 2.5) + propensity_c + (1 | part),
                data = trust)
print(summary(contrast(emmeans(m_trust, ~ format), method = CONTRASTS, adjust = "holm"),
              infer = TRUE))

# ---- B10. Robustness and sensitivity analyses --------------------------------
d_cons <- droplevels(d[d$role == "consultant", ])
d_hmc  <- droplevels(d[d$site == "HMC", ])
sens <- list("adjusted for controls" = c(paste(PRIMARY, "+", CONTROLS), "d"),
             "carry-over"            = c(paste(PRIMARY, "+ prev_format"), "d"),
             "consultants only"      = c(PRIMARY, "d_cons"),
             "HMC only"              = c(PRIMARY, "d_hmc"))
for (y in names(results)) for (s in names(sens)) {
  f <- as.formula(paste(y, "~", sens[[s]][1], "+ (1 | part) + (1 | case)"))
  ms <- tryCatch(fit_glmer(f, sens[[s]][2]), error = function(e) NULL)
  if (is.null(ms)) { cat("\n", y, "/", s, ": could not be fitted\n"); next }
  cat("\n", y, "/", s, "\n")
  check_fit(ms, s)
  if (s == "adjusted for controls" && has_pkg("performance"))
    tryCatch(print(performance::check_collinearity(ms)), error = function(e) print(e))
  if (s == "carry-over") {
    m_base <- fit_glmer(as.formula(paste(y, "~", PRIMARY, "+ (1 | part) + (1 | case)")))
    cat("  carry-over likelihood-ratio test p =",
        signif(anova(m_base, ms)[["Pr(>Chisq)"]][2], 3), "\n")
  }
  print(planned_contrasts(ms, y)[4, ])                      # text vs heatmap
}

# ---- B11. Phase 2: participant profiles for choosing interviewees -----------
prof <- aggregate(cbind(accept_correct, reject_incorrect) ~ part, data = d,
                  FUN = mean, na.action = na.pass, na.rm = TRUE)
prof <- merge(prof, participants[, c("part", "role", "site", "training", "lextale",
                                     "experience")], by = "part")
prof[, 2:3] <- round(100 * prof[, 2:3])
write.csv(prof[order(prof$reject_incorrect), ],
          "results/interview_selection_profiles.csv", row.names = FALSE)

# ---- B12. Numbers for Table 2 and Section 5.3 --------------------------------
lab <- c(none = "No explanation", heatmap = "Heatmap", text = "Text",
         both = "Heatmap plus text")
for (f in FORMATS) {
  r <- desc[desc$format == f, ]
  cat(sprintf("%s & %.0f\\%% (%d/%d) & %.0f\\%% (%d/%d) \\\\\n", lab[f],
              r$pct_accept, r$events_accept, r$n_accept,
              r$pct_reject, r$events_reject, r$n_reject))
}
pt <- results$reject_incorrect$contrasts
pt <- pt[pt$contrast == "text vs heatmap", ]
p_txt <- if (pt$p_raw < 0.001) "p < .001" else sprintf("p = %.3f", pt$p_raw)
cat(sprintf("Primary test: OR %.2f, 95%% CI %.2f-%.2f, %s\n",
            pt$OR, pt$CI_low, pt$CI_high, p_txt))
sessionInfo()

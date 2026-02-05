
rm(list=ls())
#  specify the DGP
dgp <- "null"
J_truth <- 5000
this_date <- '_2feb2026'
file_truth <- paste0('truth_dgp_', dgp, '_J', J_truth, this_date, '.Rdata')
load(file_truth)
psi <- as.numeric(truth[1] - truth[2])


format_round <- function(x, rounder){
  format(round(as.numeric(x), rounder), nsmall = rounder)
}

# with a little help from ChatGPT


format_me <- function(J, psi, est1, pt, ci_lo, ci_hi, se, pval){
  
  data.frame(
    Stage1 = est1,
    J      = J,
    pt_ci  = paste0(
      "$",
      format_round(mean(pt * 100), 2), "\\;(",
      format_round(mean(ci_lo * 100), 2), ",\\;",
      format_round(mean(ci_hi * 100), 2),
      ")$"
    ),
    bias    = format_round(mean((pt - psi) * 100), 2),
    mc_sd   = format_round(sd(pt), 4),
    ave_se  = format_round(mean(se), 4),
    cover   = format_round(mean(ci_lo <= psi & psi <= ci_hi) * 100, 1),
    reject  = format_round(mean(pval < 0.05) * 100, 1),
    stringsAsFactors = FALSE
  )
}



Js <- c(20, 30, 50, 70)

OUT_ALL <- do.call(
  rbind,
  lapply(Js, function(J){
    
    load(paste0(dgp, "_J", J, "_rep1000_2feb2026.Rdata"))
    
    estimator_map <- list(
      gee              = list(obj = gee,              label = "GEE"),
      stmle_eligible   = list(obj = stmle_eligible,   label = "Single-Stage TMLE"),
      stmle_ratio      = list(obj = stmle_ratio,      label = "Single-Stage Ratio TMLE"),
      screened         = list(obj = screened,         label = "Screened"),
      eligible         = list(obj = eligible,         label = "Eligible"),
      unadj            = list(obj = unadj,            label = "Unadjusted"),
      tmle             = list(obj = tmle,             label = "TMLE"),
      tmle_tmle        = list(obj = tmle_tmle,        label = "TMLE–TMLE"),
      tmle_tmle_cv     = list(obj = tmle_tmle_cv,        label = "TMLE–TMLE-CV")
    )
    do.call(
      rbind,
      lapply(estimator_map, function(e){
        
        format_me(
          J = J,
          psi = psi,
          est1 = e$label,
          pt = e$obj$est,
          ci_lo = e$obj$CI.lo,
          ci_hi = e$obj$CI.hi,
          se = e$obj$se,
          pval = e$obj$pval
        )
        
      })
    )
  })
)

OUT_ALL$J <- as.integer(OUT_ALL$J)

OUT_ALL <- OUT_ALL[order(OUT_ALL$Stage1, OUT_ALL$J), ]
row.names(OUT_ALL) <- NULL

library(xtable)

xt <- xtable(
  OUT_ALL,
  align = c("l","l","l","c","c","c","c","c","c"),
  caption = "Estimator performance across sample sizes",
  label = "tab:simulation_results"
)

addtorow <- list(
  pos = list(-1),
  command = paste0(
    "\\textbf{Method} & ",
    "\\textbf{J} & ",
    "\\textbf{Pt (95\\% CI)} & ",
    "\\textbf{Bias} & ",
    "\\boldmath$\\sigma$ & ",
    "\\boldmath$\\hat{\\sigma}$ & ",
    "\\textbf{Coverage} & ",
    "\\textbf{Power} \\\\ \n"
  )
)

print(
  xt,
  include.rownames = FALSE,
  sanitize.text.function = identity,
  hline.after = c(-1, 0, nrow(OUT_ALL)),
  add.to.row = addtorow
)


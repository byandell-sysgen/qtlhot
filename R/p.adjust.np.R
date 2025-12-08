p.adjust.np <- function(tests, method = "BH")
{
  for(test.name in paste("pvals.np", c("AIC","BIC"), sep = "."))
    for (i in seq(ncol(tests[[test.name]])))
      tests[[test.name]][, i] <- stats::p.adjust(tests[[test.name]][, i],
                                                 method = method)
  tests
}

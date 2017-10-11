pdf("../Figures/lmoutput")

library(xtable)

cat("Call:
lm(formula = y ~ x1 + x2 + x3 + x4 - 1)

    Residuals:
    Min      1Q  Median      3Q     Max
    -5.7479 -1.1102 -0.0531  0.9590  5.3176

    Coefficients:
    Estimate Std. Error t value Pr(>|t|)
    x1  0.05505    0.05158   1.067    0.286
    x2  0.04689    0.05013   0.935    0.350
    x3 -0.05049    0.05030  -1.004    0.316
    x4 -0.01418    0.04958  -0.286    0.775

    Residual standard error: 1.585 on 996 degrees of freedom
    Multiple R-squared:  0.003085,	Adjusted R-squared:  -0.0009184
    F-statistic: 0.7706 on 4 and 996 DF,  p-value: 0.5444")

dev.off()

About

DARPA Drone Position Error Analysis | Math 536 Exam 1

End-to-end statistical consulting report analyzing camera-based position estimation error in a DARPA drone system. Covers the full analyst workflow: exploratory data analysis, predictor transformation, OLS regression, rigorous assumption diagnostics, and robust inference under violated assumptions.

Highlights:

- Identified heteroscedasticity and non-normality via Breusch-Pagan (BP = 216.74, p < .001) and Shapiro-Wilk tests — and correctly pivoted to nonparametric methods rather than ignoring violations
- Implemented case-resampling BCa bootstrap (B = 10,000) as a principled remedy for assumption failures
- Delivered actionable correction factors with honest uncertainty quantification at 0.1, 1, and 10 km
= Communicated technical findings to both non-technical (executive summary) and technical (appendix) audiences

Methods:
Log-linear regression · OLS estimation · Assumption testing · Nonparametric bootstrap · BCa confidence intervals
Tools: R — boot, lmtest, car, ggplot2

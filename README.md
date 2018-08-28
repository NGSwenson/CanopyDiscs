# CanopyDiscs
A collaborative project with Jenny Zambrano and Bill Fagan

# Order of analyses
1) canopy_overlap_plot.R code: to calculate crown overlaps

2) ProcessOverlapFiles.R code: go through all output files to generate a table of conspecific and heterospecific overlaps and a community data matrix

3) Jaccard.Bray.Observed.R code: calculate the similaritiy of heterospecific overlaps between all conspecifics

4) Jaccard.Randomizations.R code: shuffle names on trees in plot and recalculate #3 999 times

5) make.SES.R: calculate p-values and standardized effect sizes using output from #3 and #4

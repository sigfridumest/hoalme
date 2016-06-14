#fitting a model
WaferSubset.fit <- lmeObject(current ~ poly(voltage, 2), random = list(Wafer = pdDiag(~ 1 + poly(voltage, 1)),
		Site = pdDiag(~ 1 + poly(voltage, 1))), data = WaferSubset)
WaferSubset.fit
# hypotheses: a list per level
hyp.scales <- list(fixed = 2.7, reStruct = list( list(c(2.2, 2.2), c(1.2, .5)), list(c(2.4, 2.3), c(2.4, 2.2) ) ) )
#constrained maximization
cons.WaferSubset.fit <- conslmeObjects(WaferSubset.fit, scale = hyp.scales)
#plots
plot(cons.WaferSubset.fit) # unconstrained parameterization
plot(cons.WaferSubset.fit, FALSE) # natural paramaterization
#plotting and computing confidence intervals using r
r.WaferSubset <- signedRoot(cons.WaferSubset.fit)
#plot
plot(r.WaferSubset)
# .99 % c CI's
intervals(r.WaferSubset, .99)

#higher order analysis
#computing r*
rs.WaferSubset <- RStar(WaferSubset.fit, cons.WaferSubset.fit)
#extracting Q
q.WaferSubset <- QStatistic(rs.WaferSubset)
#plot for Q
plot(q.WaferSubset)
#plot comparing r* and r
plot(rs.WaferSubset, compare = TRUE)
#plot for r* only
plot(rs.WaferSubset, compare = FALSE)

#using hoalme method
hl.WaferSubset <- hoalme(WaferSubset.fit, cons.WaferSubset.fit)
#or extracting from an RStar object
hl.WaferSubset <- hoalme(rs.WaferSubset)
#tests of hypotheses for fixed effects coefficients and CI's for  every parameter
summary(hl.WaferSubset)

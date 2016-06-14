#fitting a model
Orth.fit <- lmeObject(distance ~ 1 + Sex + I(age - 11), random = pdSymm(~ 1 + I(age - 11)), data = OrthSub)
Orth.fit
#hypotheses
hyp.scales <- list(fixed = list(c(3, 3), c(3, 3) , c(3, 3)), reStruct = list( list(c(2.6, .75), c(.1, 2), c(1.2, .5)) ) )
#constrained maximization
cons.Orth.fit<- conslmeObjects(Orth.fit, scale = hyp.scales)
#profiles
plot(cons.Orth.fit) #plot, unconstrained parameterization
plot(cons.Orth.fit, FALSE) #plot, natural parameterization

#computing r*
rs.Orth <- RStar(Orth.fit, cons.Orth.fit)
#plot comparing r* and r
plot(rs.Orth, compare = TRUE)
#updating r* computation setting a smaller  neighborhood of the hypothesis
rs.Orth.d05 <- update(rs.Orth, .5)
#exhibits some odd points
plot(rs.Orth.d05, compare = TRUE)


intervals(rs.Orth)
sr.Orth <- signedRoot(rs.Orth)
intervals(sr.Orth)
hl.Orth <- hoalme(rs.Orth)
hl.Orth
summary(hl.Orth)


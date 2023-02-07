data=read.table("amarlx",header=TRUE)

burn_in=20

# convert vals to numeric
for (colname in colnames(data)){
  a=strsplit(data[,colname],split="D")
  data[,colname]=sapply(a, function(x) as.numeric(x[1])*10^as.numeric(x[2]))
}

plot(x = data$t,
     y = data$RA,
     xlab = "time (dimensionless)",
     ylab = "Reaction rate Aragonite",type="l",
     xlim=c(burn_in,max(data$t)),
     ylim=range(data$RA[data$t>burn_in]))

plot(x = data$t,
     y = data$RC,
     xlab = "time (dimensionless)",
     ylab = "Reaction rate Calcite",
     type="l",
     xlim=c(burn_in,max(data$t)),
     ylim=range(data$RC[data$t>burn_in]))

plot(x = data$t,
     y = data$OC,
     xlab = "time (dimensionless)",
     ylab = "Oversaturation Calcite",
     type="l",
     xlim=c(burn_in,max(data$t)),
     ylim=range(data$OC[data$t>burn_in]))

plot(x = data$t,
     y = data$OA,
     xlab = "time (dimensionless)",
     ylab = "Oversaturation Aragonite",
     type="l",
     xlim=c(burn_in,max(data$t)),
     ylim=range(data$OA[data$t>burn_in]))

SimpleExpFunc <- function(x, b0, b1, b2) {1 / (b0 - b1*exp(-x/b2))}
FitSimpleExp <- function(data, fit.range, flat.points) {
	# User inversed log linear fit to guess the start point.
	inversed.data <- data.frame(x = data[, 1], y = 1/data[, 2])
	flat.value <- median(inversed.data[flat.points, 2])
	# log transform.
	log.inversed.y <- log(flat.value - inversed.data[1:6, 2])
	log.inversed.data <- data.frame(x = inversed.data[1:6, 1], y = log.inversed.y)
	# build linear model and estimate parameters.
	linear.model <- lm(y ~ x, log.inversed.data, na.action = na.exclude)
	intercept <- coef(linear.model)[[1]]
	slope <- coef(linear.model)[[2]]
  
  if (is.na(slope) | is.na(intercept)) {
    return(NULL)
  }

	# now we got the start point.
	b1 = exp(intercept)
	b2 = -1/slope
	b0 = flat.value
  
  fit <- NULL
  try(fit <- nls(y ~ SimpleExpFunc(x, b0, b1, b2), data = data[fit.range, ],
                 start = list(b0 = b0, b1 = b1, b2 = b2)))
  return(fit)
}

# Predict 
Predict <- function(x.range, b0, b1, b2) {
	x = x.range
	y = sapply(x, SimpleExpFunc, b0 = b0, b1 = b1, b2 = b2)
	data.frame(x = x, y = y)
}

library(RJSONIO)
library(ggplot2)
library(ggthemes)
Plot <- function(data.dir, prefix, func, pos, tofit=T) {
	filename <- paste(prefix, func, paste("pos", pos, ".json", sep = ""), sep = "_")
	filepath <- file.path(data.dir, filename)
  print(filepath)
	data <- fromJSON(filepath)
	plot.x.range <- 1:50
	plot.data <- data.frame(x = data[["CtIndices"]], y = data[["Ct"]])
	plot.data <- plot.data[plot.data$x %in% plot.x.range, ]

	g <- ggplot(plot.data, aes(x = x, y = y)) + geom_point() + theme_few() +
	  xlab("Codon distance") + ylab(expression("C"["t"]))

	# Fitting
	  if (tofit) {
	    fit.range <- 1:50
	    flat.points <- 25:50
	    nlm <- FitSimpleExp(plot.data, fit.range, flat.points)
      if (!is.null(nlm)) {
        b0 = coef(nlm)[[1]]
        b1 = coef(nlm)[[2]]
        b2 = coef(nlm)[[3]]
        predict.data <- Predict(plot.x.range, b0, b1, b2)
        g <- g + geom_line(data = predict.data, aes(x = x, y = y))
        
        ks <- data$Ks
        ratio <- b1*(1+ks*4/(4-1))/(2*ks*(b0+b1))
        
        print(summary(nlm))
        print(ks)
        # annotate the plot.
        g <- g + ggtitle(sprintf("%s\nr/m = %.3f", prefix, ratio))
        fitted = TRUE
      } else {
        fitted = FALSE
      }
	    
	  } else {
	    fitted = FALSE
	  }
  
  if (!fitted) {
    g <- g + ggtitle(sprintf("%s", prefix)) + geom_line()
  }
 	g
}

Wrap <- function(prefix, data.dir, func, pos) {
	gplot <- Plot(data.dir, prefix, func, pos)
	outname <- paste(prefix, ".png", sep = "")
	outfile <- file.path("plots", outname)
	ggsave(filename = outfile, plot = gplot, width = 6, height = 4)
}

data.dir <- "cov"
func <- "CovReads"
pos <- 4
pattern <- paste(func, paste("pos", pos, ".json", sep = ""), sep = "_")
files <- list.files(data.dir, pattern = pattern)
prefixes <- gsub(paste("", pattern, sep = "_"), "", files)
for (prefix in prefixes) {
	Wrap(prefix, data.dir, func, pos)
}
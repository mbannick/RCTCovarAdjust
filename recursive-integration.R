## In what follows I am assuming the effect-size is 0!
### This would have to be modified for constructing confidence intervals


## Defining the grid on which we want to evaluate the densities
lower = -10
upper = 10
l.out = 100000
grid_x <- seq(from = lower, to = upper, length.out = l.out)
delta = (upper - lower)/l.out

### This should be the cutoff for significance from the first block (on the score scale I believe)
### (It shouldn't be 2! It should be a bit bigger, depending on our design, I just threw in 2 because it seemed semi-reasonable)
cutoff <- 2
cutoff <- 2.943058

### Here we set up the truncated distribution of the statistic from the first block.
### I standardize it (to make it a dist), though that is actually not necessary
### For future blocks this will change!
gridded_base_dist <- dnorm(grid_x)
gridded_base_dist_continue <- gridded_base_dist
gridded_base_dist_continue[which(abs(grid_x) >= cutoff)] <- 0
gridded_base_dist_continue <- gridded_base_dist_continue / (sum(gridded_base_dist_continue)*delta)

### This is the distribution of the new increment
### I believe I am implicitly doing things on the score scale for equally sized blocks
### Otherwise I should probably modify the SD below
gridded_new_increment <- dnorm(grid_x)

### Here I convolve the two distributions
### The convolve function in R is quite fast and nice!
### but it flips things in a funny way (so I flip it back)
est_flipped <- convolve(gridded_base_dist_continue,
                        gridded_new_increment)  * delta
gridded_convolved_dist <- est_flipped[c(ceiling(length(grid_x)/2):length(grid_x),
                               1:(ceiling(length(grid_x)/2)-1))]


plot(grid_x, gridded_base_dist_continue) ## the dist of the statistic at the previous stage
plot(grid_x, gridded_new_increment) ## the dist of the new increment
plot(grid_x, gridded_convolved_dist) ## the dist of the sum
sum(gridded_convolved_dist)*delta ## checking to make sure this is a dist

## doing this in function form

recursive_integ <- function(gridded_base_dist,
                            gridded_new_increment,
                            delta,
                            cutoff){
  gridded_base_dist_continue <- gridded_base_dist
  gridded_base_dist_continue[which(abs(grid_x) >= cutoff)] <- 0
  gridded_base_dist_continue <- gridded_base_dist_continue / (sum(gridded_base_dist_continue)*delta)

  est_flipped <- convolve(gridded_base_dist_continue,
                          gridded_new_increment)  * delta
  gridded_convolved_dist <- est_flipped[c(ceiling(length(grid_x)/2):length(grid_x),
                                          1:(ceiling(length(grid_x)/2)-1))]
  return(gridded_convolved_dist)
}

## Testing the function version

approx_dist <-recursive_integ(gridded_base_dist,
                              gridded_new_increment,
                              delta,
                              cutoff)
plot(grid_x,approx_dist, type='l')

## We can also do this again with the output for the second stage to get the sampling distribution for the 3rd
## I suspect I'm doing the cutoff wrong (trying to convert from the score scale to the likelihood ratio scale)

cutoff_2 <- cutoff*sqrt(2)
cutoff_2 <- 1.987662*sqrt(2)

gridded_base_dist_2 <- approx_dist

approx_dist_2 <-recursive_integ(gridded_base_dist_2,
                              gridded_new_increment,
                              delta,
                              cutoff_2)
plot(grid_x,approx_dist_2, type='l')

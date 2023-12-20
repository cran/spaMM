cat(crayon::yellow("\ntest-pdep-effects.R: "))

# Easily goes wrong with factors -> fixes for the following checks in v4.4.4

## Simple -> OK
fit <- fitme(Petal.Length ~ Sepal.Length + Species, data = iris)
plot_effects(fit, "Species")

## Drop value of one level without dropping the level itself 
iris2 <- iris[iris$Species != "setosa", ]
fit2 <- fitme(Petal.Length ~ Sepal.Length + Species, data = iris2)
plot_effects(fit2, "Species")

## Non alphabetical level orders
iris3 <- iris
iris3$Species <- factor(iris3$Species, levels = c("versicolor",
                                                  "virginica", "setosa"))
fit3 <- fitme(Petal.Length ~ Sepal.Length + Species, data = iris3)
plot_effects(fit3, "Species")

## Factors as character vector (not true factors) 
iris4 <- iris
iris4$Species <- as.character(iris4$Species)
fit4 <- fitme(Petal.Length ~ Sepal.Length + Species, data = iris4)
plot_effects(fit4, "Species")
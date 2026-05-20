requiredPackages <- c(
  "dplyr",
  "splines",
  "parallel",
  "pbapply",
  "fastDummies",
  "SMUT",
  "nlshrink"
)

for (p in requiredPackages) {
  if (!require(p, character.only = TRUE))
    install.packages(p)
  library(p, character.only = TRUE)
}
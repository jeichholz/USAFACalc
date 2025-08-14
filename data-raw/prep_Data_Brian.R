#This script prepares more data that Brian gave me for R.

ThreeBladeProp_T53=readr::read_csv("data-raw/ThreeBladeProp_T53 1.csv")
TurboNormalizer_T53=readr::read_csv("data-raw/TurboNormalizer_T53 1.csv")
Explosion1=readr::read_csv("data-raw/Explosion1.csv")

#Save as .rda in the proper directory. Add depends information to the project.
usethis::use_data(ThreeBladeProp_T53,overwrite=TRUE)
usethis::use_data(TurboNormalizer_T53,overwrite=TRUE)
usethis::use_data(Explosion1,overwrite=TRUE)

#Add documentation file that I'll fill out later.
usethis::use_r("ThreeBladeProp_T53_doc")
usethis::use_r("TurboNormalizerT53_doc")
usethis::use_r("Explosion1_doc")

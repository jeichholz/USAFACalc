#This script prepares the T53 climb data for inclusion in USAFACalc.

USFoodExpenditure=readr::read_csv("data-raw/USFoodExpenditure.csv")
SunSpotNumber2=readr::read_csv("data-raw/SunSpotNumber2.csv")
MedianHomePrice=readr::read_csv("data-raw/MedianHomePrice.csv")

#Save as .rda in the proper directory. Add depends information to the project.
usethis::use_data(USFoodExpenditure,overwrite=TRUE)
usethis::use_data(SunSpotNumber2,overwrite=TRUE)
usethis::use_data(MedianHomePrice,overwrite=TRUE)

#Add documentation file that I'll fill out later.
usethis::use_r("USFoodExpenditure_doc")
usethis::use_r("SunSpotNumber2_doc")
usethis::use_r("MedianHomePrice_doc")

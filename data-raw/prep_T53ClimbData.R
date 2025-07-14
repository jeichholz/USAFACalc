#This script prepares the T53 climb data for inclusion in USAFACalc.

T53ClimbData1=readr::read_csv("data-raw/T53ClimbData1.csv")
T53ClimbData2=readr::read_csv("data-raw/T53ClimbData2.csv")

#Save as .rda in the proper directory. Add depends information to the project.
usethis::use_data(T53ClimbData1,overwrite=TRUE)
usethis::use_data(T53ClimbData2)

#Add documentation file that I'll fill out later.
usethis::use_r("T53ClimbData1_doc")
usethis::use_r("T53ClimbData2_doc")

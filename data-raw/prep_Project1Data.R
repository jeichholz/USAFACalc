#This script prepares the T53 climb data for inclusion in USAFACalc.

pwt_africa=readr::read_csv("data-raw/pwt_africa.csv")
pwt_asia=readr::read_csv("data-raw/pwt_asia.csv")
pwt_euro=readr::read_csv("data-raw/pwt_euro.csv")

#Save as .rda in the proper directory. Add depends information to the project.
usethis::use_data(pwt_africa,overwrite=TRUE)
usethis::use_data(pwt_asia,overwrite=TRUE)
usethis::use_data(pwt_euro,overwrite=TRUE)

#Add documentation file that I'll fill out later.
usethis::use_r("pwt_africa_doc")
usethis::use_r("pwt_asia_doc")
usethis::use_r("pwt_euro_doc")

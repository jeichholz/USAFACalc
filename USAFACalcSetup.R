#Installation of required packages for Calc 1/2.

#mosaic.
install.packages("mosaic",quiet = TRUE)

#mosiacCore
install.packages("mosaicCore",quiet = TRUE)

#mosaicData
install.packages("mosaicData",quiet = TRUE)

#mosaicCalc
install.packages("mosaicCalc",quiet = TRUE)

#Manipulate
install.packages("manipulate",quiet=TRUE)

#MMAC
install.packages("MMAC",quiet=TRUE)

#caracas
install.packages("caracas",quiet = TRUE)

#set up sympy
caracas::install_sympy()

#devtools
install.packages("devtools",quiet = TRUE)

#usafacalc
devtools::install_github("jeichholz/USAFACalc",quiet = TRUE)

#Write the .Rprofile
outFile=file.path(Sys.getenv("HOME"),".Rprofile")

writeLines(c("#We are running .First.sys() here, because that loads all the default packages.  R runs .Rprofile before it loads the default packages.",
             "#So if you don't run this first, you'll load mosaic, but then you'll load other packages which mask some of what you want from mosaic.",
             "#This way we load the defaults first, load mosaic second, and let mosaic mask the defaults.",
             ".First.sys()",
             "",
             "if (!suppressPackageStartupMessages(library(mosaic,warn.conflicts=FALSE,logical.return=TRUE))){",
             "  warning(\"Error loading package mosaic.\")",
             "}",
             "",
             "if (!suppressPackageStartupMessages(library(mosaicCore,warn.conflicts=FALSE,logical.return=TRUE))){",
             "  warning(\"Error loading package mosaicCore.\")",
             "}",
             "",
             "if (!suppressPackageStartupMessages(library(mosaicData,warn.conflicts=FALSE,logical.return=TRUE))){",
             "  warning(\"Error loading package mosaicData.\")",
             "}",
             "",
             "if (!suppressPackageStartupMessages(library(manipulate,warn.conflicts=FALSE,logical.return=TRUE))){",
             "  warning(\"Error loading package manipulate.\")",
             "}",
              "",
              "if (!suppressPackageStartupMessages(library(MMAC,warn.conflicts=FALSE,logical.return=TRUE))){",
              "  warning(\"Error loading package MMAC.\")",
              "}",
              "",
              "if (!suppressPackageStartupMessages(library(mosaicCalc,warn.conflicts=FALSE,logical.return=TRUE))){",
              "  warning(\"Error loading package mosaicCalc.\")",
              "}",
              "",
              "if (!suppressPackageStartupMessages(library(USAFACalc,warn.conflicts=FALSE,logical.return=TRUE))){",
              "  warning(\"Error loading package USAFACalc.\")",
              "}"),outFile)



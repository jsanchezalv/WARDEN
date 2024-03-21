# script to update the documentation and vignettes hosted at https://https://github.com/javiersanchezalvarezevidera/RDICE/
# for each package once a change is made update the following files then run the appropriate section of code
# this will pre-compile all the vignettes and documentation to the website

##############################################
# Updates to add a new package
# when adding a new package is necessary to include links/make updates in the following files as well as in this script:
# all existing _pkgdown.yml files (update the navbar section for other packages)

# 1) make any updates to the package 
# 2) update documentation using devtools
devtools::document()

# 3) reinstall the updated package 
devtools::install()

# 4) update the file _pkgdown.yml
# 5) rebuild documentation using pkgdown
pkgdown::build_site()

# Note as the vignettes for some packages can take a long time to run it is
# also possible just to update partial sections by just running the 
# below functions without regenerating the vignettes
pkgdown::build_home()
pkgdown::build_reference()
pkgdown::build_articles_index()



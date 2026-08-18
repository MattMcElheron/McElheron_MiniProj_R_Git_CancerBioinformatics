
## Module 1: Project Setup & Version Control

In this module, we will set up a structured computational environment and establish version control using Git and GitHub. Organizing your code, data, and outputs into standard subdirectories from the outset ensures your project remains structured, while linking your RProject to GitHub ensures all updates are tracked, backed up, and fully reproducible. We want our code to work time and time again on our machine, AND to be useable in the exact same way on someone else's computer (a colleague, a reviewer, a client). Similar to lab work, we want a strict lab notebook with a very specific layout and set of instructions, and for it to actually be readable for others! Think of recipe books - they usually follow the same structure - we do not want to be on the final step of baking our cake to realize we are missing a vital ingredient, coding follows the same logic.

### Step 1: Create local directory structure
Open **RStudio**, create a new R Project linked to your cloned GitHub repository, and run the following in your Console:

```r
# Create standard project folders
dir.create("data")
dir.create("scripts")
dir.create("output")

```

Now, when someone else opens our script, they will know how our data/input/output should be structured.

# Script: 01_setup.R
# Purpose: Check and install required packages

When we start using RStudio, we will need to install packages, using `install.packages("package_name")`. We then switch on the package, using `library(package_name)`. Note the use of quotations for install but not for library. The next time we open up RStudio, we do not need to reinstall the package, but we do need to switch it on, with the `library()` function.  

The below code gives R the list of packages we need and installs them if they are missing. 

```r
required_packages <- c("dplyr", "ggplot2", "tidyverse", "ggrepel")

new_packages <- required_packages[!(required_packages %in% installed_packages)]
if(length(new_packages)) install.packages(new_packages)


```


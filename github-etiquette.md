### GitHub Etiquette for BioTIME users

Keep your scripts in your personal scripts folder. Don't look in other people's folders.

When working on group projects, move the script out of your personal folder to the relevant project folder, so that everyone can work on it.

Keep file paths short and sensible.

Upload your data on GitHub and use relative file paths - if the data are on your computer, and you have e.g. `data <- read.csv("C:/user/Documents/bird_abundance.csv")` in your code, only you will be able to read in that data, since it's only on your computer. But if you load the data from the lab's repo, every member of the lab can use it, e.g. `data <- read.csv("data/bird_abundance.csv")`, where `data` is a folder within the lab's repo.

Don't use funky characters and spaces in your file names, these cause trouble because of differences in Mac/Windows systems.

Always pull before you push in case someone has done any work since the last time you pulled - you wouldn't want anyone's work to get lost or to have to resolve many coding conflicts.

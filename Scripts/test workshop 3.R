#create annotations

#set working dir
setwd("C:/NRSC_510/WorkshopR-master")

#identify working dir
getwd()

#windows, right click, copy path

# ?setwd() gives information on the function if you aren't sure what to do

covid_data_newcol <- mutate(covid_data, ###)
                            
covid_data_newcol <- mutate(covid_data, ###)
                            
                            set.seed(2020)
                            student_grades <- data.frame(ID = 1:1000,
                                                         grade = rnorm(1000, 80, 5))
                            authors <- data.frame(
                              surname = I(c("Tukey", "Venables", "Tierney", "Ripley", "McNeil")),
                              nationality = c("US", "Australia", "US", "UK", "Australia"),
                              deceased = c("yes", rep("no", 4)))
                            
                            authorN <- within(authors, { name <- surname; rm(surname) })
                            
                            books <- data.frame(
                              name = I(c("Tukey", "Venables", "Tierney",
                                         "Ripley", "Ripley", "McNeil", "R Core")),
                              title = c("Exploratory Data Analysis",
                                        "Modern Applied Statistics ...",
                                        "LISP-STAT",
                                        "Spatial Statistics", "Stochastic Simulation",
                                        "Interactive Data Analysis",
                                        "An Introduction to R"),
                              other.author = c(NA, "Ripley", NA, NA, NA, NA,
                                               "Venables & Smith"))
                            
# test
# What happens if I add a line?

# Database connection
install.packages("RMySQL")
library(RMySQL)

db <- dbConnect(RMySQL::MySQL(),
                dbname = "company",
                host="courses.csrrinzqubik.us-east-1.rds.amazonaws.com",
                post=3306,
                user="student",
                password="datacamp"
)

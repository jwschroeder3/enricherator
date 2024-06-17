# Load the required package
library(testthat)

# Use the function test_dir to run all the tests inside 'tests' directory
result <- test_dir("tests/") 

# view the result
print(result)

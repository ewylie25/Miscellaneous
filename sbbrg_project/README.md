## SBBRG Correlation Project

###Description:
1. Verify a high level of correlation between AAPL and GOOG by calculating the correlation coefficient using daily prices over a 30-day period in January 2012.  
2. Make the assumption that GOOG should match the daily percent change of AAPL.  Create formatted output over a 29-day period in February 2012 with columns for date, stock prices, and expected stock price for GOOG based on the daily movement for AAPL.  Over this period, determine when the difference between the expected stock price and actual stock price is greater than the 1% multiplied by the correlation coefficient.  When the price is below its expected value highlight the value in green and when the price is above its expected value highlight the value in red.  
3. The data can come from a MySQL table that you create, a csv that you generate, or scraped from a website, such as Yahoo Finance. Use a public source, such as Yahoo Finance or Google to gather pricing data.  Provide a summary of the methodology used.  


####Files:

- main.php
- Database.php
- CorrelationCalculator.php
- DailyPercentChangeCalculator.php
- helper_functions.php
- parameters.json
- Calculator.php
- CorrelationCalculatiorTest.php
- style.css
- output.html

####Main

**Usage**: From CLI just call without arguments:`php main.php`  
**Dependencies**: Requires package sbbrg_project, parameters.json, style.css, PECL Stats package, and PDO_MYSQL driver.  
**Implementation Notes**: Uses the percent daily change to calculate the correlation coefficient between two companies. I was slightly unclear what was meant by 1% of the correlation coefficent since 0.50 * 0.01 is a very small number and the difference is never less than that value.  

####Database

**Contains**: Class Database  
**Dependencies**: PDO_MYSQL driver  
**Implementation Notes**: Creates persistent connection using PHP data object, see doc strings for specific implementation notes. Member of sbbrg_project package.  
**Unit Test**: TODO  

####CorrelationCalculator

**Contains**: Class CorrelationCalculator  
**Dependencies**: Class Calculator, PECL Stats package  
**Implementation Notes**: Simplified interface to calculate correlation coefficients from tabular data, see doc strings for specifics. Member of sbbrg_project.  
**Unit Test**: CorrelationCalculatorTest.php (phpunit)  

####DailyPercentChangeCalculator

**Contains**: Class DailyPercentChangeCalculator  
**Dependencies**: Class Calculator  
**Implementation Notes**: Simplified interface to calculate daily percent change from tabular data, see doc strings for specifics. Member of sbbrg_project.  
**Unit Test**: TODO  

####Calculator

**Contains**: Class Calculator  
**Dependencies**: Class Database  
**Implementation Notes**: Class to simplify and consolidate interface to My SQL database, see doc strings for specifics. Member of sbbrg_project.  
**Unit Test**: TODO  

####helper_functions

**Contains**: generate_DB_table_from_URL  
**Dependencies**: Class Database  
**Implementation Notes**: Generates My SQL table from URL. See doc strings and comments for more info. Member of sbbrg_project.  
**Unit Test**: N/A  

####parameters.json

JSON file containing parameters for project. Two major sets of parameters:  
- My SQL: username, password, database, and host ip
- Data: URL, line delimiter, field delimiter, and field definitions - name and My SQL type.

The user for My SQL must have read/write access to specified database and the database must already exist. The data URLs are links to CSV files containing price data from Jan-1-2012 to Feb-29-2012 for GOOG and AAPL, respectively, from Google Finance Historical Prices. I believe/hope these are persistent links.  

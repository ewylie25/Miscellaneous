## SBBRG Correlation Project

####Description:  
1. Verify a high level of correlation between AAPL and GOOG by calculating the correlation coefficient using daily prices over a 30-day period in January 2012.  
2. Make the assumption that GOOG should match the daily percent change of AAPL.  Create formatted output over a 29-day period in February 2012 with columns for date, stock prices, and expected stock price for GOOG based on the daily movement for AAPL.  Over this period, determine when the difference between the expected stock price and actual stock price is greater than the 1% multiplied by the correlation coefficient.  When the price is below its expected value highlight the value in green and when the price is above its expected value highlight the value in red.  
3. The data can come from a MySQL table that you create, a csv that you generate, or scraped from a website, such as Yahoo Finance. Use a public source, such as Yahoo Finance or Google to gather pricing data.  Provide a summary of the methodology used.  


####Files:

-main.php
-Database.php
-CorrelationCalculator.php
-DailyPercentChangeCalculator.php
-helper_functions.php  
-parameters.json
-Calculator.php
-CorrelationCalculatiorTest.php  

#####main

Usage:  
Dependencies:  
Implementation Notes:  

#####Database

Contains:
Dependencies:
Implementation Notes:
Unit Test:

#####CorrelationCalculator

Contains:  
Dependencies:  
Implementation Notes:  
Unit Test: CorrelationCalculatorTest.php (phpunit)  

#####DailyPercentChangeCalculator

Contains:  
Dependencies:  
Implementation Notes:  
Unit Test:  

#####Calculator

Contains:  
Dependencies:  
Implementation Notes:  
Unit Test:  

#####helper_functions

Contains:  
Dependencies:  
Implementation Notes:  
Unit Test: N/A  

#####parameters.json

My SQL: user@'%' with password, granted all privileges on database.*  make database first...

Data: CSV files containing price data from Jan-1-2012 to Feb-29-2012 for GOOG and AAPL, respectively, from Google Finance Historical Prices. I believe/hope these are persistent links.  

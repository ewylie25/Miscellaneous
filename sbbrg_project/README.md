## SBBRG Correlation Project

####Description:  
1. Verify a high level of correlation between AAPL and GOOG by calculating the correlation coefficient using daily prices over a 30-day period in January 2012.  
2. Make the assumption that GOOG should match the daily percent change of AAPL.  Create formatted output over a 29-day period in February 2012 with columns for date, stock prices, and expected stock price for GOOG based on the daily movement for AAPL.  Over this period, determine when the difference between the expected stock price and actual stock price is greater than the 1% multiplied by the correlation coefficient.  When the price is below its expected value highlight the value in green and when the price is above its expected value highlight the value in red.  
3. The data can come from a MySQL table that you create, a csv that you generate, or scraped from a website, such as Yahoo Finance. Use a public source, such as Yahoo Finance or Google to gather pricing data.  Provide a summary of the methodology used.  


####Files:

-main.php
-CorrelationCalculator.php
-helper_functions.php
-parameters.json
-CorrelationCalculatiorTest.php




My SQL $user@'%' with $pass, granted all privileges on test.* used for this script.

containing price data from Jan-1-2012 to Feb-29-2012 for GOOG and AAPL, respectively, from Google Finance Historical Prices.
Returns file contents as string - I believe/hope these are persistent links

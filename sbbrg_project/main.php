<?php
/**
 * User: liz
 * Date: 5/20/14
 *
 * Usage: CLI>php main.php
 * Input: parameters.json (implicit)
 * Output:
 */

require 'helper_functions.php';
require 'AbstractCalculator.php';
require 'CorrelationCalculator.php';
require 'DailyPercentChangeCalculator.php';
require 'Database.php';

# Loading parameters from json file
#   $raw -> String of parameters.json contents
#   $parameters -> Associative Array of parameters.json contents (assuming proper json in file, otherwise null)
$raw = file_get_contents('parameters.json');
$parameters = json_decode($raw, true);

# Connect to specified database
#   $database -> Database set using My SQL parameters in parameters.json
$database = new \sbbrg_project\Database($parameters["My SQL"]);
$database->connect();

# Generating My SQL database tables from specified data if table does not already exist.
#   $company -> String of data (company) name specified in parameters.json
#   $exists -> Boolean denoting if table named for $company already exists in database
#   $success -> Boolean denoting if generate_database_from_CSV completed without Exception
foreach(array_keys($parameters["Data"]) as $company){
    $exists = $database->isTable($company);
    if (!$exists){
        $success = \sbbrg_project\generate_database_from_CSV($database, $parameters["Data"][$company], $company);
        if (!$success){
            die("Could not load data specified in parameters.json for: ".$company."\n");
        }
    }
}

# Calculate correlation of prices of Apple and Google Stocks
#   $correlation_calc -> CorrelationCalculator using $database
$correlation_calc = new \sbbrg_project\CorrelationCalculator($database);

# Set companies to use in calculating correlation. This will accept more than two company names, but will only use the first two.
$correlation_calc->setCompanies(array('google', 'apple'));

# Calculate correlation for January and February using opening stock prices
$correlation_calc->setCorrelation('open', "2012-01-01", "2012-01-31");
$jan_open_corr = $correlation_calc ->getCorrelation();
$correlation_calc->setCorrelation('open', "2012-02-01", "2012-02-29");
$feb_open_corr = $correlation_calc ->getCorrelation();

# Calculate correlation for January and February using closing stock prices
$correlation_calc->setCorrelation('close', "2012-01-01", "2012-01-31");
$jan_close_corr = $correlation_calc ->getCorrelation();
$correlation_calc->setCorrelation('close', "2012-02-01", "2012-02-29");
$feb_close_corr = $correlation_calc ->getCorrelation();

# Calculate correlation for January and February using daily low stock prices
$correlation_calc->setCorrelation('low', "2012-01-01", "2012-01-31");
$jan_low_corr = $correlation_calc->getCorrelation();
$correlation_calc->setCorrelation('low', "2012-02-01", "2012-02-29");
$feb_low_corr = $correlation_calc ->getCorrelation();

# Calculate correlation for January and February using daily high stock prices
$correlation_calc->setCorrelation('high', "2012-01-01", "2012-01-31");
$jan_high_corr = $correlation_calc->getCorrelation();
$correlation_calc->setCorrelation('high', "2012-02-01", "2012-02-29");
$feb_high_corr = $correlation_calc ->getCorrelation();

# Calculate daily percent change in stock prices for Google and Apple
#   $percent_change_calc -> DailyPercentChangeCalculator using $database
$percent_change_calc = new \sbbrg_project\DailyPercentChangeCalculator($database);

# Set companies to use in calculating daily percent change. This will accept and use any number of companies (assuming data in $database).
$percent_change_calc->setCompanies(array('google', 'apple'));

# Get dates of interest
#   $date -> Array of Strings of business days in date range
$dates = $percent_change_calc->setDateRange("2012-02-01","2012-02-29");

# Build data for table
#   $data -> Associative Array containing table content
#   $day -> String date
#   $apple_percent_change -> Float of the percent change in stock price that day
#   $google_closing_price
$data = array();
foreach($dates as $day){
    $percent_change_calc->setDailyPercentChange('apple', $day);
    $apple_percent_change =  $percent_change_calc->getDailyPercentChange();
    $googleclose = $database->getValue('google', $day, 'close');
    $googleopen = $database->getValue('google', $day, 'open');
    $appleclose = $database->getValue('apple', $day, 'close');
    $expectedgoog = (1+$apple_percent_change)*$googleopen;
    $diff = abs($expectedgoog - $googleclose);
    $onepercentcorr =  0.01*$feb_close_corr;
    if ($diff > $onepercentcorr){
        $diffovercorr = true;
    } else{
        $diffovercorr = false;
    }
    if ($googleclose>$expectedgoog){
        $color = "#FF0000";
    } else{
        $color = "#00FF00";
    }
    $data["$day"] = array("google close" => $googleclose,
                          "apple close" => $appleclose,
                         "expected google" => $expectedgoog,
                          "color" => $color,
                          "Difference Greater than 1% of correlation" =>$diffovercorr);
}

$database->disconnect();

# TODO: Decide on output - static html page for viewing and dynamic page for server...
print "January\n";
print "Opening Prices correlation: ".$jan_open_corr."\n";
print "Closing Prices correlation: ".$jan_close_corr."\n";
print "Low Price correlation: ".$jan_low_corr."\n";
print "High Price correlation: ".$jan_high_corr."\n";
print "February\n";
print "Opening Prices correlation: ".$feb_open_corr."\n";
print "Closing Prices correlation: ".$feb_close_corr."\n";
print "Low Price correlation: ".$feb_low_corr."\n";
print "High Price correlation: ".$feb_high_corr."\n";

print_r($data);
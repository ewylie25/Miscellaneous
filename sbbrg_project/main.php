<?php
/**
 * Liz Wylie
 * May 2014
 *
 * Contains main method for sbbrg correlation project. Procedure is broken up into seven steps and variables
 *  are documented in each section of code.
 *
 * TODO: Implement Logging
 * TODO: Stricter Error handling?
 * TODO: Unit tests...
 *
 * Usage: CLI>php main.php
 * Input: parameters.json (implicit)
 * Output: output.html
 */

/*
 * Supporting functions and classes in namespace sbbrg_project
 * - helper_functions.php -> function generate_DB_table_from_URL
 * - CorrelationCalculator.php -> class CorrelationCalculator
 * - DailyPercentChangeCalculator.php -> class DailyPercentChangeCalculator
 * - Database.php -> class Database
 */
require 'helper_functions.php';
require 'CorrelationCalculator.php';
require 'DailyPercentChangeCalculator.php';
require 'Database.php';

/*
 * 1. Loading parameters from json file
 * - $raw -> String of parameters.json contents
 * - $parameters -> Associative Array of parameters.json contents (assuming proper json in file, otherwise null)
 */

$raw = file_get_contents('parameters.json');
$parameters = json_decode($raw, true);

/*
 * 2. Connect to specified database
 * - $database -> Database set using My SQL parameters in parameters.json
 */

$database = new \sbbrg_project\Database($parameters["My SQL"]);
$database->connect();

/*
 * 3. Generating My SQL database tables from specified data if table does not already exist.
 * - $company -> String, data (company) name specified in parameters.json
 * - $exists -> Boolean, denotes if table named for $company already exists in database
 * - $success -> Boolean, denotes if generate_database_from_CSV completed without Exception
 */

foreach(array_keys($parameters["Data"]) as $company){
    // Check if table exists in database
    $exists = $database->isTable($company);
    // Generate table if it doesn't work
    if (!$exists){
        $success = \sbbrg_project\generate_DB_table_from_URL($database, $parameters["Data"][$company], $company);
        if (!$success){
            // Quit if cannot load data
            die("Could not load data specified in parameters.json for: ".$company."\n");
        }
    }
}

/*
 * 4. Calculate correlation of prices of Apple and Google Stocks
 * - $correlation_calc -> CorrelationCalculator using $database
 * - $correlation -> Float, correlation coefficient calculated using percent daily change between Apple and Google stocks in January 2012
 */

$correlation_calc = new \sbbrg_project\CorrelationCalculator($database);

// Set companies to use in calculating correlation. This will accept more than two company names, but will only use the first two.
$correlation_calc->setCompanies(array('google', 'apple'));

// Set dates of interest
$correlation_calc->setDateRange("2012-01-01", "2012-01-31");

// Calculate correlation using percent daily change
$correlation_calc->setCorrelation('pdc');
$correlation = $correlation_calc ->getCorrelation();

/*
 * 5. Build data for table
 * - $percent_change_calc -> DailyPercentChangeCalculator, uses $database
 * - $dates -> Array of Strings of business days in date range
 * - $data -> Associative Array, contains table content
 * - $day -> String, date
 * - $apple_percent_change -> Float, the percent change in Apple stock price that day
 * - $google_closing_price ->05/20/14 Float, the closing price for Google stock that day
 * - $google_opening_price -> Float, the opening price for Google stock that day
 * - $apple_closing_price -> Float, the closing price for Apple stock that day
 * - $expected_google_price -> Float, the expected closing price of Google stock assuming matched with Apple stock's daily percent change that day
 * - $diff -> Float, absolute difference between actual and expected closing price of Google stock that day
 * - $one_percent_correlation -> Float, ...
 * - $is_diff_gt_opc -> Boolean, ...
 */

// Calculate daily percent change in stock prices for Google and Apple
$percent_change_calc = new \sbbrg_project\DailyPercentChangeCalculator($database);

// Set companies to use in calculating daily percent change. This will accept and use any number of companies (assuming data in $database).
$percent_change_calc->setCompanies(array('google', 'apple'));

// Get dates of interest
$dates = $percent_change_calc->setDateRange("2012-02-01","2012-02-29");

$data = array();
// Add data for each business day in February 2012
foreach($dates as $day){
    // Get Apple's percent change that day
    $percent_change_calc->setDailyPercentChange('apple', $day);
    $apple_percent_change =  $percent_change_calc->getDailyPercentChange();

    // Calculate the expected price given Apple's percent change
    $google_opening_price = $database->getValue('google', $day, 'open');
    $expected_google_price = (1+$apple_percent_change)*$google_opening_price;

    // Data for the table
    $google_closing_price = $database->getValue('google', $day, 'close');
    $apple_closing_price = $database->getValue('apple', $day, 'close');

    // Determine if the difference was greater than one percent of the correlation coefficient
    $diff = abs($expected_google_price - $google_closing_price);
    $one_percent_correlation =  0.01*$correlation;
    if ($diff > $one_percent_correlation){
        $is_diff_gt_opc = true;
    } else{
        $is_diff_gt_opc = false;
    }

    // Color the cells green if less than expected and red if greater than expected
    if ($google_closing_price>$expected_google_price){
        $color = "#FF0000";
    } else{
        $color = "#00FF00";
    }

    // Build data structure
    $data[$day] = array("date" => $day,
                        "google close" => $google_closing_price,
                        "apple close" => $apple_closing_price,
                        "expected google" => $expected_google_price,
                        "color" => $color,
                        "difference" => $diff,
                        "Difference Greater than 1% of correlation" =>($is_diff_gt_opc)? 'True' : 'False');
}

/*
 * 6. Disconnect from the database
 */
$database->disconnect();

/*
 * 7. Write output file
 * - $content -> string, HTML content
 * - $filename -> string, HTML file name
 * - $success -> bool, denoting if write to file worked
 * TODO: Integrate code into writing content, it would be more efficient for this application. Maybe use DOM.
 */
// HTML - uses a default style sheet I've used before
$content = '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <meta http-equiv="content-type" content="text/html; charset=utf-8" />
        <meta name="description" content="Output page from main.php in SBBRG Correlation Project" />
        <meta name="keywords" content="" />
        <meta name="author" content="Liz Wylie" />
        <link rel="stylesheet" type="text/css" href="style.css" media="screen" />
        <title>Correlation Project</title>
    </head>
    <div id="wrapper">
        <div id="header">
            <h2>SBB Research Group Project 1</h2>
        </div>
        <body>
            <div id="content">
                <h3>Part 1</h3>
                <p>
                    January correlation coefficient calculated using percent daily change is '.round($correlation, 4).'
                </p>
                <h3>Part 2</h3>
                <p>
                    I added all of the values asked to be calculated because I thought the prompt was ambiguous. Red denotes Actual GOOG price was greater than expected and Green denotes it was less than expected.
                </p>
                <table>
                    <tr>
                        <td>Date</td>
                        <td>APPL Close</td>
                        <td>GOOG Close</td>
                        <td>Expected GOOG Close</td>
                        <td>GOOG Difference</td>
                        <td>Greater than 1% Corr</td>
                    </tr>';

foreach($data as $d){
    $content .= '
                       <tr>';
    $content .='
                            <td>'.$d["date"].'</td>';
    $content .='
                            <td>'.$d["apple close"].'</td>';
    $content .='
                            <td>'.$d["google close"].'</td>';
    $content .='
                            <td style="background-color:'.$d['color'].'">'.round($d["expected google"], 3).'</td>';
    $content .='
                            <td>'.round($d["difference"], 3).'</td>';
    $content .='
                            <td>'.$d["Difference Greater than 1% of correlation"].'</td>';
    $content .= '
                       </tr>';
}
$content .= '
                </table>
                <h3>Part 3</h3>
                <p>
                    The data URLs are links to CSV files containing price data from Jan-1-2012 to Feb-29-2012 for GOOG and AAPL from Google Finance Historical Prices. I believe/hope these are persistent links. I worked on a 64-bit Ubuntu 12.04 LTS desktop using a php 5.5 interpreter from ppa:ondrej/php5 and My SQL 5.5. As setup, I created the test database. I created a auser@\'%\' with password, and  granted all privileges on test.*. I can provide more information if necessary.
                </p>
            </div>
            <div id="footer">
                <p>
                    <a href="https://github.com/ewylie25/Miscellaneous/tree/master/sbbrg_project">See it on GitHub</a>
                </p>
            </div>
        </body>
    </div>
</html>';

// Set file name
$filename ="output.html";

// Write to file
$success = file_put_contents($filename, $content);
if (!$success){
    die("Failed to write HTML file.");
}
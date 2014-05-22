<?php
/**
 * User: liz
 * Date: 5/20/14
 */

require 'helper_functions.php';
require 'AbstractCalculator.php';
require 'CorrelationCalculator.php';
require 'DailyPercentChangeCalculator.php';

# Loading parameters from json file
$raw = file_get_contents('parameters.json');
$parameters = json_decode($raw, true);

# Generating My SQL database tables from csv data specified in json file
$success = \sbbrg_project\generate_database_from_CSV($parameters);
if (!$success){
    die("Could not generate databases with given parameters in parameters.json");
}

# Part 1 - Verifying high correlation for January prices of Apple and Google Stocks
$calc = new \sbbrg_project\CorrelationCalculator($parameters);
$calc->setCorrelation('open', "2012-01-01", "2012-01-31");
$open_corr = $calc ->getCorrelation();
$calc->setCorrelation('close', "2012-01-01", "2012-01-31");
$close_corr = $calc ->getCorrelation();
$calc->setCorrelation('low', "2012-01-01", "2012-01-31");
$low_corr = $calc->getCorrelation();
$calc->setCorrelation('high', "2012-01-01", "2012-01-31");
$high_corr = $calc->getCorrelation();

# Part 2
$perc = new \sbbrg_project\DailyPercentChangeCalculator($parameters);
$dates = $perc->setDateRange("2012-02-01","2012-02-29");

$calc->setCorrelation('open', "2012-02-01", "2012-02-29");
$feb_open = $calc ->getCorrelation();
$calc->setCorrelation('close', "2012-02-01", "2012-02-29");
$feb_close = $calc ->getCorrelation();
$calc->setCorrelation('high', "2012-02-01", "2012-02-29");
$feb_high = $calc ->getCorrelation();
$calc->setCorrelation('low', "2012-02-01", "2012-02-29");
$feb_low = $calc ->getCorrelation();

$data = array();
foreach($dates as $day){
    $perc->setDailyPercentChange('apple', $day);
    $applechange =  $perc->getDailyPercentChange();
    $googleclose = $perc->getValue('google', $day, 'close');
    $googleopen = $perc->getValue('google', $day, 'open');
    $appleclose = $perc->getValue('apple', $day, 'close');
    $expectedgoog = (1+$applechange)*$googleopen;
    $diff = abs($expectedgoog - $googleclose);
    $onepercentcorr =  0.01*$feb_close;
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

# TODO: Decide on output - static html page for viewing and dynamic page for server...
//print "January\n";
//print "Opening Prices correlation: ".$open_corr."\n";
//print "Closing Prices correlation: ".$close_corr."\n";
//print "Low Price correlation: ".$low_corr."\n";
//print "High Price correlation: ".$high_corr."\n";
//print "February\n";
//print "Opening Prices correlation: ".$feb_open."\n";
//print "Closing Prices correlation: ".$feb_close."\n";
//print "Low Price correlation: ".$feb_low."\n";
//print "High Price correlation: ".$feb_high."\n";
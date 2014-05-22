<?php
/**
 * User: liz
 * Date: 5/20/14
 */

require 'helper_functions.php';
require 'CorrelationCalculator.php';

$raw = file_get_contents('parameters.json');
$parameters = json_decode($raw, true);

$success = \sbbrg_project\generate_database_from_CSV($parameters);
if (!$success){
    die("Could not generate databases with given parameters in parameters.json");
}

print "Opening Prices correlation: ";
$calc = new \sbbrg_project\CorrelationCalculator($parameters);
$open_corr = $calc ->getCorrelation('open');
print $open_corr."\n";
print "Closing Prices correlation: ";
$close_corr = $calc ->getCorrelation('close');
print $close_corr."\n";
print "Low Price correlation: ";
$low_corr = $calc->getCorrelation('low');
print $low_corr."\n";
print "High Price correlation: ";
$high_corr = $calc->getCorrelation('high');
print $high_corr."\n";
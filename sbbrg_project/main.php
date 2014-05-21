<?php
/**
 * Created by PhpStorm.
 * User: liz
 * Date: 5/20/14
 * Time: 12:32 PM
 */

include 'helper_functions.php';
include 'correlationcalculation.php';

$success = \sbbrg_project\generate_database_from_CSV();
if (!$success){
    die("Could not generate databases");
}

print "Open"."\n";
$open = new \sbbrg_project\CorrelationCalculation('open');
$open_corr = $open ->getCorrelation();
print $open_corr."\n";

print "Close"."\n";
$close = new \sbbrg_project\CorrelationCalculation('close');
$close_corr = $close ->getCorrelation();
print $close_corr."\n";

print "Low"."\n";
$low = new \sbbrg_project\CorrelationCalculation('low');
$low_corr = $low->getCorrelation();
print $low_corr."\n";

print "High"."\n";
$high = new \sbbrg_project\CorrelationCalculation('high');
$high_corr = $high->getCorrelation();
print $high_corr."\n";
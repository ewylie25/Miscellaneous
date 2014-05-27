<?php
/**
 * Liz Wylie
 * May 2014
 *
 * Class for calculating correlation coefficients. See doc strings.
 */

namespace sbbrg_project;
require_once('Calculator.php');

/**
 * Class CorrelationCalculator
 *
 * Class to  calculate correlation coefficients.
 *
 * Simple interface to calculate the correlation between two companies using any value in the data set or percent
 * daily change. Company names and values available are described in parameters.json.
 *
 * @package sbbrg_project
 * @see Calculator
 *
 * @param float $correlation The last calculated correlation coefficient.
 * @param string $type The type of the last calculated correlation coefficient.
 */
class CorrelationCalculator extends Calculator {
    private $correlation;
    private $type;

    /**
     * Returns last calculated correlation coefficient.
     *
     * @return float    Calculated correlation coefficient.
     */
    public function getCorrelation(){
        return $this->correlation;
    }

    /**
     * Calculates and stores correlation coefficient.
     *
     * @param string $type Select value to use in calculating correlation coefficient.
     */
    public function setCorrelation($type){
        // Store type
        $this->type = $type;
        // Instantiate data structure
        $processed_data = array();

        // Deal with type -> price of percent daily change, "pdc" or column name from table
        switch($type){
            case "pdc":
                // Iterate through first two companies specified and creates array of percent daily change values.
                for($i = 0; $i<2; $i++){
                    $processed_data[$i] = array_map(function($x){return ($x["close"]-$x["open"])/$x["open"];}, $this->data[$this->companies[$i]]);
                }
                break;
            default:
                // Iterate through first two companies specified and creates array of values from specified column.
                for($i = 0; $i<2; $i++){
                    $processed_data[$i] = array_map(function($x){return $x[$this->type];}, $this->data[$this->companies[$i]]);
                }
        }

        // Calculate correlation using stats extension
        $this->correlation = stats_stat_correlation($processed_data[0], $processed_data[1]);
    }

} 
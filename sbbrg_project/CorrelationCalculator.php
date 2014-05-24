<?php
/**
 * User: liz
 * Date: 5/20/14
 */

namespace sbbrg_project;
require_once('Calculator.php');

class CorrelationCalculator extends Calculator {

    private $correlation;

    public function getCorrelation(){
        return $this->correlation;
    }
    public function setCorrelation($type, $startdate, $enddate){
        $this->type = $type;
        $this->start = $startdate;
        $this->end = $enddate;

        # retrieve data from My SQL
        $this->setData();

        $processed_data = array();

        # Deal with type -> price of percent daily change
        switch($type){
            case "pdc":
                for($i = 0; $i<2; $i++){
                    $processed_data[$i] = array_map(function($x){return ($x["close"]-$x["open"])/$x["open"];}, $this->data[$this->companies[$i]]);
                }
                break;
            default:
                for($i = 0; $i<3; $i++){
                    $processed_data[$i] = array_map(function($x){return $x[$this->type];}, $this->data[$this->companies[$i]]);
                }
        }

        # Calculate correlation using stats extension
        $this->correlation = stats_stat_correlation($processed_data[0], $processed_data[1]);
    }

} 
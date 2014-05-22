<?php
/**
 * User: liz
 * Date: 5/20/14
 */

namespace sbbrg_project;


class CorrelationCalculator extends AbstractCalculator {

    private $correlation;

    public function getCorrelation(){
        return $this->correlation;
    }
    public function setCorrelation($price_type, $startdate, $enddate){
        $this->price_type = $price_type;
        $this->start = $startdate;
        $this->end = $enddate;
        # retrieve data from My SQL
        $this->setData();
        # Calculate correlation using stats extension
        $this->correlation = stats_stat_correlation($this->data[0], $this->data[1]);
    }
    protected function setData(){
        for($i=0; $i < count($this->datasets); $i++){
            $this->data[$i] = array_map(function($x){return $x[$this->price_type];} ,$this->fromDB($this->datasets[$i]));
        }
    }

} 
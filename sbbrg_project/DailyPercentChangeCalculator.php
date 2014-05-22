<?php
/**
 * User: liz
 * Date: 5/22/14
 */

namespace sbbrg_project;


class DailyPercentChangeCalculator extends AbstractCalculator {
    private $percentchange;

    public function getDailyPercentChange(){
        return $this->percentchange;
    }
    public function setDailyPercentChange($company, $date){
        #TODO: Check Date in Date Range
        #TODO: Catch double dates or missing dates
        foreach($this->data[$company] as $temp){
            if ($temp['date'] == $date){
                $this->percentchange = ($temp["close"] - $temp["open"])/$temp["open"];
                break;
            }
        }
    }
    public function setDateRange($startdate, $enddate){
        $this->start = $startdate;
        $this->end = $enddate;

        $this->setData();

        $begin = new \DateTime($startdate);
        $end = new \DateTime($enddate);
        $end = $end->modify( '+1 day' );
        $interval = new \DateInterval('P1D');
        $daterange = new \DatePeriod($begin, $interval ,$end);

        $temp = array();
        foreach ($daterange as $date) {
            $day_num = $date->format("N");
            $day = $date->format('Y-m-d');
            if ($day_num < 6 and $day !== "2012-02-20"){
                array_push($temp, $day);
            }
        }
        return $temp;
    }
    protected function setData(){
        for($i=0; $i < count($this->datasets); $i++){
            $this->data[$this->datasets[$i]] = $this->fromDB($this->datasets[$i]);
        }
    }


} 
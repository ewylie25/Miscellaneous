<?php
/**
 * Liz Wylie
 * May 2014
 *
 * Class for calculating daily percent change. See doc strings.
 */

namespace sbbrg_project;
require_once('Calculator.php');

/**
 * Class DailyPercentChangeCalculator
 *
 * Class to calculate the percent daily change.
 *
 * Simple interface to calculate the percent daily change for a company on a given day in a specified date range. Company
 * names and data sets are described in parameters.json.
 *
 * @package sbbrg_project
 * @see Calculator
 *
 * @param float $percent_change The last calculated daily percent change.
 */
class DailyPercentChangeCalculator extends Calculator {
    private $percent_change;

    /**
     * Returns last calculated daily percent change.
     *
     * @return float    Calculated daily percent change.
     */
    public function getDailyPercentChange(){
        return $this->percent_change;
    }

    /**
     * Calculates and stores percent daily change.
     *
     * @param string $company Name of company as specified in parameters.json
     * @param string $date Date of interest formatted as yyyy-mm-dd.
     */
    public function setDailyPercentChange($company, $date){
        // TODO: Check Date in Date Range
        // TODO: Catch double dates or missing dates
        // Return data for date of interest
        $temp = $this->data[$company][$date];

        // Calculate percent change.
        $this->percent_change = ($temp["close"] - $temp["open"])/$temp["open"];
    }

    /**
     * Sets date range and returns array of dates of business days in given range
     *
     * Overrides method in parent class.
     *
     * @param string $start_date Starting date of interest in yyyy-mm-dd format.
     * @param string $end_date Last date of interest in yyyy-mm-dd format.
     * @return array|void Array of business days (string format) in date range
     */
    public function setDateRange($start_date, $end_date){
        // Assign dates
        $this->start = $start_date;
        $this->end = $end_date;

        // retrieve data from My SQL
        $this->setData();

        // use built in date functions to return all days in range
        $begin = new \DateTime($start_date);
        $end = new \DateTime($end_date);
        $end = $end->modify( '+1 day' );
        $interval = new \DateInterval('P1D');
        $date_range = new \DatePeriod($begin, $interval ,$end);

        // filter out holidays and weekends
        $temp = array();
        foreach ($date_range as $date) {
            $day_num = $date->format("N");
            $day = $date->format('Y-m-d');
            if ($day_num < 6 and $day !== "2012-02-20"){
                array_push($temp, $day);
            }
        }
        // return business days
        return $temp;
    }
} 
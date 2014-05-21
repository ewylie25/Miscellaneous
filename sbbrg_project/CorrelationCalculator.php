<?php
/**
 * User: liz
 * Date: 5/20/14
 */

namespace sbbrg_project;


class CorrelationCalculation {
    public $price_type;
    private $correlation;
    private $google_data;
    private $apple_data;

    public function __construct($price_type){
        # Set column of data to look at
        $this->price_type = $price_type;
    }
    public function getCorrelation(){
        $this->setCorrelation();
        return $this->correlation;
    }
    private function setCorrelation(){
        # retrieve data from My SQL
        $this->setData();

        # Calculate correlation using stats extension
        $this->correlation = stats_stat_correlation($this->google_data, $this->apple_data);
    }
    private function setData(){
        $this->google_data = $this->fromDB('google');
        $this->apple_data = $this->fromDB('apple');
    }
    private function fromDB($table){
        #TODO:Unify db parameters in settings file
        $user = 'auser';
        $pass = 'apass';
        $temp = array();
        try {
            $dbh = new \PDO('mysql:host=localhost;dbname=test', $user, $pass);
            foreach($dbh->query("SELECT $this->price_type FROM $table WHERE `DATE` BETWEEN '2012-01-01' AND '2012-01-31'") as $value) {
                array_push($temp, $value[0]);
            }
            $dbh = null;
        } catch (\PDOException $e) {
            print "Error!: " . $e->getMessage() . "<br/>";
        }
        return $temp;
    }
} 
<?php
/**
 * User: liz
 * Date: 5/20/14
 *
 * TODO: PHP Statistics seems poorly documented/developed: proper way to integrate R/python with exec()?
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

        # I was a little concerned about the stats extension because it isn't very documented... totally works fine
        //$google_sd = $this->getStdDev($this->google_data);
        //$apple_sd = $this->getStdDev($this->apple_data);
        //$covariance = $this->getCovariance($this->google_data, $this->apple_data);
        //$this->correlation = $covariance/($google_sd*$apple_sd);

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
            foreach($dbh->query("SELECT $this->price_type from $table") as $value) {
                array_push($temp, $value[0]);
            }
            $dbh = null;
        } catch (\PDOException $e) {
            print "Error!: " . $e->getMessage() . "<br/>";
        }
        return $temp;
    }
    # I was a little concerned about the stats extension because it isn't very documented... totally works fine/ at least comparable to R
   //private function getStdDev($data){
   //     return stats_standard_deviation($data);
   // }
   // private function getCovariance($data1, $data2){
   //     return stats_covariance($data1, $data2);
   // }
} 
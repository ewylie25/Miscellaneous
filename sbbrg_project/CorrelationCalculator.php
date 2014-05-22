<?php
/**
 * User: liz
 * Date: 5/20/14
 */

namespace sbbrg_project;


class CorrelationCalculator {
    public $price_type;
    public $host;
    public $database;
    public $user;
    public $pass;
    private $correlation;
    private $datasets;
    private $data;

    public function __construct($parameters){
        $this->host = $parameters["My SQL"]["host"];
        $this->database = $parameters["My SQL"]["database"];
        $this->user = $parameters["My SQL"]["user"];
        $this->password = $parameters["My SQL"]["password"];
        $this->datasets = array_keys($parameters["Data"]);
    }
    public function getCorrelation($price_type){
        $this->price_type = $price_type;
        $this->setCorrelation();
        return $this->correlation;
    }
    private function setCorrelation(){
        # retrieve data from My SQL
        $this->setData();

        # Calculate correlation using stats extension
        $this->correlation = stats_stat_correlation($this->data[0], $this->data[1]);
    }
    private function setData(){
        for($i=0; $i < count($this->datasets); $i++){
            $this->data[$i] = $this->fromDB($this->datasets[$i]);
        }
    }
    private function fromDB($table){
        $temp = array();
        try {
            $dbh = new \PDO("mysql:host=$this->host;dbname=$this->database", $this->user, $this->password);
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
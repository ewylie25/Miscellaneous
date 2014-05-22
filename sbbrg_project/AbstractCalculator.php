<?php
/**
 * Created by PhpStorm.
 * User: liz
 * Date: 5/22/14
 * Time: 11:48 AM
 */

namespace sbbrg_project;


abstract class AbstractCalculator {
    public $host;
    public $database;
    public $user;
    public $pass;
    public $datasets;
    protected $data;
    protected $price_type;
    protected $start;
    protected $end;

    public function __construct($parameters){
        $this->host = $parameters["My SQL"]["host"];
        $this->database = $parameters["My SQL"]["database"];
        $this->user = $parameters["My SQL"]["user"];
        $this->password = $parameters["My SQL"]["password"];
        $this->datasets = array_keys($parameters["Data"]);
    }

    public function getValue($company, $date, $column){
        foreach($this->data[$company] as $temp){
            if ($temp['date'] == $date){
                return $temp[$column];
            }
        }
        return null;
    }

    abstract protected  function setData();

    protected function fromDB($table){
        $temp = array();
        try {
            $dbh = new \PDO("mysql:host=$this->host;dbname=$this->database", $this->user, $this->password);
            foreach($dbh->query("SELECT * FROM $table WHERE `DATE` BETWEEN \"$this->start\" AND \"$this->end\"") as $value) {
               array_push($temp, $value);
            }
            $dbh = null;
        } catch (\PDOException $e) {
            print "Error!: " . $e->getMessage() . "\n";
        }
        return $temp;
    }
} 
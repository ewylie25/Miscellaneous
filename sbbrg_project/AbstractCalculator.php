<?php
/**
 * Created by PhpStorm.
 * User: liz
 * Date: 5/22/14
 * Time: 11:48 AM
 */

namespace sbbrg_project;


abstract class AbstractCalculator {
    public $database;
    protected $companies;
    protected $data;
    protected $price_type;
    protected $start;
    protected $end;

    public function __construct($database){
        $this->database = $database;
    }

    abstract protected  function setData();

    public function setCompanies($companies){
        $this->companies = array();
        foreach ($companies as $company){
            $exists = $this->database->isTable($company);
            if (!$exists){
                die("Company not in database: ".$company."\n");
            } else{
                array_push($this->companies,$company);
            }
        }
    }

    protected function fromDB($table){
        return $this->database->getTable($table, $this->start, $this->end);
    }
} 
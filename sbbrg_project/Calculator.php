<?php
/**
 * Liz Wylie
 * May 2014
 *
 */

namespace sbbrg_project;


class Calculator {
    public $database;
    protected $companies;
    protected $data;
    protected $type;
    protected $start;
    protected $end;

    public function __construct($database){
        $this->database = $database;
    }

    protected  function setData(){
        for($i=0; $i < count($this->companies); $i++){
            $this->data[$this->companies[$i]] = $this->fromDB($this->companies[$i]);
        }
    }

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
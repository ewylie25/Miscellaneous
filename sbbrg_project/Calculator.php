<?php
/**
 * Liz Wylie
 * May 2014
 *
 * Helper parent class - see doc strings. Inherited by CorrelationCalculator and DailyPercentChangeCalculator.
 */

namespace sbbrg_project;

/**
 * Class Calculator
 *
 * Parent class to reduce code redundancy in calculations from data sets.
 *
 * @package sbbrg_project
 *
 * @var Database $database My SQL database used to retrieve data.
 * @var array $companies The names of companies(string) of interest that have tables in specified database.
 * @var array $data Map of company names(string) to the data retrieved from the database(multidimensional array).
 * @var string $start Starting date of interest in data sets in yyyy-mm-dd format.
 * @var string $end Last date of interest in data sets in yyyy-mm-dd format.
 */
class Calculator {
    protected $database;
    protected $companies;
    protected $data;
    protected $start;
    protected $end;

    /**
     * Class constructor.
     *
     * Sets database for instance.
     *
     * @param Database $database
     */
    public function __construct($database){
        $this->database = $database;
    }

    /**
     * Creates list of companies and verifies corresponding tables exist in Database.
     *
     * @param array $companies Array of company names as strings. Naming should correspond to names set in parameters.json.
     */
    public function setCompanies($companies){
        // Instantiate data structure
        $this->companies = array();

        // Iterate through company names
        foreach ($companies as $company){
            // Retrieve boolean denoting if table exists in database named as $company
            $exists = $this->database->isTable($company);
            // If it doesn't exist then something is wrong - exit code, otherwise add to class parameter.
            if (!$exists){
                die("Company not in database: ".$company."\n");
            } else{
                array_push($this->companies,$company);
            }
        }
    }

    /**
     * Sets date range of interest to calculator.
     *
     * @param string $start_date Starting date of interest in data sets in yyyy-mm-dd format.
     * @param string $end_date Last date of interest in data sets in yyyy-mm-dd format
     */
    public function setDateRange($start_date, $end_date){
        // set dates
        $this->start = $start_date;
        $this->end = $end_date;

        // retrieve data from My SQL
        $this->setData();
    }
    /**
     * Internal function.
     *
     * Iterates through set companies building an associative array mapping company names to
     * multidimensional arrays of data retrieved from Database instance.
     */
    protected  function setData(){
        // Iterate through all the companies specified when calling setCompanies.
        for($i=0; $i < count($this->companies); $i++){
            // C
            $this->data[$this->companies[$i]] = $this->fromDB($this->companies[$i]);
        }
    }

    /**
     * Internal function.
     *
     * Called to retrieve data table for specified company from Database instance.
     *
     * @param string $table company name (table name)
     * @return array    Multidimensional (associative) array reflective of database table structure.
     */
    protected function fromDB($table){
        return $this->database->getTable($table, $this->start, $this->end);
    }
} 
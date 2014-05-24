<?php
/**
 * Liz Wylie
 * May 2014
 *
 * Class for interfacing with My SQL database.
 *
 *
 *
 */

namespace sbbrg_project;


class Database {
    public $host;
    public $db;
    public $user;
    public $pass;
    private $dbh;
    private $tables;
    private $stmt;


    public function __construct($db_parameters){
        $this->host = $db_parameters["host"];
        $this->db = $db_parameters["database"];
        $this->user = $db_parameters["user"];
        $this->pass = $db_parameters["password"];
        $this->tables = array();
    }

    public function connect(){
        try{
            $this->dbh = new \PDO("mysql:host=$this->host;dbname=$this->db", $this->user, $this->pass, array(\PDO::ATTR_PERSISTENT => true));
            print "Connected\n";
        } catch (\PDOException $e){
            print "Unable to connect: " . $e->getMessage(). "\n";
        }
        # DB settings and connection settings and style adapted from PHP Data Objects Manual - https://php.net/manual/en/book.pdo.php
        $this->dbh->setAttribute(\PDO::ATTR_ERRMODE, \PDO::ERRMODE_EXCEPTION);

        foreach($this->dbh->query("SELECT * FROM information_schema.tables") as $value) {
            array_push($this->tables, $value[2]);
        }
    }

    public function disconnect(){
        $this->dbh = null;
    }

    public function createTable($table, $definition){
        # Check if table already exists
        if (in_array($table, $this->tables)){
            print "Table(".$table.") Already exists!\n";
        } else {
            # Create My SQL table
            $this->dbh->beginTransaction();
            $this->dbh->exec("CREATE TABLE $table ($definition)");
            $this->dbh->commit();
        }
    }
    public function startAddToTable($table, $field_string, $place_holders){
        # Prepare My SQL statement
        # Following example syntax from http://php.net/manual/en/pdostatement.execute.php
        # Probably should bind parameters, but I don't see how it could be utilized maliciously as is.
        $this->stmt = $this->dbh->prepare("INSERT INTO $table ($field_string) VALUES ($place_holders)");
    }
    public function addToTable($params){

        $this->stmt->execute($params);
    }
    public function endAddToTable(){
        $this->stmt= null;
    }
    public function getTable($table, $start, $end){
        if (!in_array($table, $this->tables)){
            print "No table corresponding to ".$table."\n";
            return null;
        }
        $temp = array();
        foreach($this->dbh->query("SELECT * FROM $table WHERE `DATE` BETWEEN \"$start\" AND \"$end\"") as $value) {
            array_push($temp, $value);
        }
        return $temp;
    }

    public function isTable($table){
        return in_array($table, $this->tables);
    }

    public function getValue($table, $date, $column){
        $temp = array();
        foreach($this->dbh->query("SELECT $column FROM $table WHERE `DATE` = \"$date\"") as $value){
            array_push($temp, $value[0]);
        }
        return $temp[0];
    }
} 
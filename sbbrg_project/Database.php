<?php
/**
 * Liz Wylie
 * May 2014
 *
 * Class for interfacing with My SQL database.
 *
 * DB settings and connection settings and style adapted from PHP Data Objects Manual - https://php.net/manual/en/book.pdo.php
 */

namespace sbbrg_project;

/**
 * Class Database
 *
 * Class to simplify and consolidate interface to My SQL database.
 *
 * Creates a persistent connection to a specified My SQL database using a PHP Data Object.
 *
 * @package sbbrg_project
 *
 * @param string $host My SQL host name. Specified in parameters.json.
 * @param string $db My SQL database name. Specified in parameters.json.
 * @param string $user My SQL user name with appropriate privileges. Specified in parameters.json.
 * @param string $pass My SQL password for corresponding username. Specified in parameters.json.
 * @param \PDO $dbh PHP Data Object interface for My SQL database
 * @param array $tables Array of database's tables.
 * @param \PDOStatement $stmt Prepared My SQL statement
 */
class Database {
    public $host;
    public $db;
    public $user;
    public $pass;
    private $dbh;
    private $tables;
    private $stmt;

    /**
     * Class constructor
     *
     * Sets parameters necessary for connecting to database.
     *
     * @param Array $db_parameters Associative array containing database parameters
     */
    public function __construct($db_parameters){
        $this->host = $db_parameters["host"];
        $this->db = $db_parameters["database"];
        $this->user = $db_parameters["user"];
        $this->pass = $db_parameters["password"];
        $this->tables = array();
    }

    /**
     * Connects to database specified in constructor.
     */
    public function connect(){
        // Establish a persistent connection, exit code if cannot
        try{
            $this->dbh = new \PDO("mysql:host=$this->host;dbname=$this->db", $this->user, $this->pass, array(\PDO::ATTR_PERSISTENT => true));
            print "Connected\n";
        } catch (\PDOException $e){
            die("Unable to connect: " . $e->getMessage());
        }
        $this->dbh->setAttribute(\PDO::ATTR_ERRMODE, \PDO::ERRMODE_EXCEPTION);
        // Add pre-existing tables to helper structure
        foreach($this->dbh->query("SELECT * FROM information_schema.tables") as $value) {
            array_push($this->tables, $value[2]);
        }
    }

    /**
     * Disconnect from database
     */
    public function disconnect(){
        $this->dbh = null;
    }

    /**
     * Create a new table
     *
     * @param string $table desired table name
     * @param string $definition definition string for My SQL, follows format "name_1 TYPE, name_2 TYPE,..."
     */
    public function createTable($table, $definition){
        # Check if table already exists; if it doesn't exist, create it
        if ($this->isTable($table)){
            print "Table(".$table.") Already exists!\n";
        } else {
            # Create My SQL table
            $this->dbh->beginTransaction();
            $this->dbh->exec("CREATE TABLE $table ($definition)");
            $this->dbh->commit();
        }
    }

    /**
     * Prepare statement for adding data to My SQL table.
     *
     * Must be proceeded by a call to endAddToTable to ensure database disconnects properly.
     * Following example syntax from http://php.net/manual/en/pdostatement.execute.php
     * TODO: Probably should bind parameters, but I don't see how it could be utilized maliciously as is.
     *
     * @param string $table Name of table to add to
     * @param string $field_string Comma-separated string of variable names being added to table
     * @param string $place_holders Comma-separated string of placeholders (?) for values being added to table
     */
    public function startAddToTable($table, $field_string, $place_holders){
        // Prepare My SQL statement
        $this->stmt = $this->dbh->prepare("INSERT INTO $table ($field_string) VALUES ($place_holders)");
    }

    /**
     * Add values to database table.
     *
     * Must be preceded by a call to startAddToTable
     *
     * @param array $params ordered list of value to add using prepared statement
     */
    public function addToTable($params){

        $this->stmt->execute($params);
    }

    /**
     * Terminate use of prepared statement.
     *
     * Must be called to ensure database connection closes at end of program if startAddToTable is called.
     */
    public function endAddToTable(){
        $this->stmt= null;
    }

    /**
     * Returns data in specified date range from table or null if table does not exist.
     *
     * @param string $table table name
     * @param string $start start date in yyyy-mm-dd format
     * @param string $end end date in yyyy-mm-dd format
     * @return array|null
     */
    public function getTable($table, $start, $end){
        // Check for table's existence
        if (!$this->isTable($table)){
            print "No table corresponding to ".$table."\n";
            return null;
        }
        // Set up associative array with date as key.
        $temp = array();
        foreach($this->dbh->query("SELECT * FROM $table WHERE `DATE` BETWEEN \"$start\" AND \"$end\"") as $value) {
            $temp[$value['date']]=$value;
        }
        return $temp;
    }

    /**
     * Check if table exists in database.
     *
     * @param string $table table name
     * @return bool     true if table exists, false if it does not
     */
    public function isTable($table){
        return in_array($table, $this->tables);
    }

    /**
     * Returns single, specified, value from table.
     *
     * @param string $table Table name.
     * @param string $date Date in yyyy-mm-dd format.
     * @param string $column Column (variable) name.
     * @return mixed    Value
     */
    public function getValue($table, $date, $column){
        $temp = array();
        foreach($this->dbh->query("SELECT $column FROM $table WHERE `DATE` = \"$date\"") as $value){
            array_push($temp, $value[0]);
        }
        return $temp[0];
    }
} 
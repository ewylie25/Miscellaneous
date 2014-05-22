<?php
/**
 * Liz Wylie
 * 05/20/14
 *
 * TODO: Implement Logging
 */

namespace sbbrg_project;


/***
 * Function to generate MySQL database for analysis from CSV files online
 * @param associative array $parameters
 * @return boolean $success
 *
 * TODO: Troubleshoot function failure - possible refactor into multiple steps
 * TODO: Rethink Design - only listed as separate function to modularize code - not really necessary
 */
function generate_database_from_CSV($parameters){
    #Initialize variable to keep track of status
    $success = false;

    # Retrieve database parameters as specified in parameters.json
    $host = $parameters["My SQL"]["host"];
    $db = $parameters["My SQL"]["database"];
    $user = $parameters["My SQL"]["user"];
    $pass = $parameters["My SQL"]["password"];

    # Establish persistent connection to database, test on localhost. Prints error and returns false if fails.
    try{
        $dbh = new \PDO("mysql:host=$host;dbname=$db", $user, $pass, array(\PDO::ATTR_PERSISTENT => true));
        echo "Connected\n";
    } catch (\PDOException $e){
        echo "Unable to connect: " . $e->getMessage(). "\n";
        return $success;
    }

    # Add data to database - close connection at end and if errors
    try{
        # DB settings and connection settings and style adapted from PHP Data Objects Manual - https://php.net/manual/en/book.pdo.php
        $dbh->setAttribute(\PDO::ATTR_ERRMODE, \PDO::ERRMODE_EXCEPTION);

        # Iterate through each data set and set table name based on which data set being used
        foreach($parameters["Data"] as $databasetable => $dataset){
            # Returns contents of data files as string as specified in parameters.json
            $csvdata = file_get_contents($dataset["url"]);

            # Set up parameters for data processing
            $fieldseparator = $dataset["field delimiter"];
            $lineseparator = $dataset["line delimiter"];

            # Create My SQL column definitions for given data set
            # Following example syntax from http://php.net/manual/en/pdostatement.execute.php
            # Probably should bind parameters, but I don't see how it could be utilized maliciously as is.
            $definition = "";
            $field_string = "";

            foreach($dataset["fields"] as $field){
                $definition = $definition." ".$field["name"]." ".$field["type"].",";
                $field_string = $field_string." ".$field["name"].",";
            }
            $definition = trim($definition, ',');
            $field_string = trim($field_string, ',');
            $place_holders = implode(',', array_fill(0, count($dataset["fields"]), '?'));

            # Create My SQL table
            $dbh->beginTransaction();
            $dbh->exec("CREATE TABLE $databasetable ($definition)");
            $dbh->commit();

            # Prepare My SQL statement
            $stmt = $dbh->prepare("INSERT INTO $databasetable ($field_string) VALUES ($place_holders)");

            # Break up text into lines (string -> array of strings)
            $lines = explode($lineseparator, $csvdata);

            # Remove the header - beware indexing is off using this method - there is no longer an element in $lines[0]
            unset($lines[0]);

            # Iterate through lines of csv file contents
            foreach($lines as $line) {
                # Make sure the string isn't empty
                if (!empty($line)){
                    # remove any artifacts from file formatting
                    $line = trim($line," \t");
                    $line = str_replace("\r","",$line);

                    $params = explode($fieldseparator, $line);

                    # Fix date formatting to be accepted by My SQL
                    $params[0] = strftime("%Y-%m-%d", strtotime($params[0]));

                    # Execute My SQL statement
                    $stmt->execute($params);
                } else {
                    # If empty sting - move on
                    continue;
                }
            }
        }
        $success = true;
        # close database connection, according to user comments both statement as necessary.
        # TODO: test behavior of disconnection
        $stmt = null;
        $dbh = null;
        return $success;

    } catch (\Exception $e){
        # if error occurs, print, close connection and exit
        echo "Error: ". $e->getMessage();
        $stmt = null;
        $dbh = null;
        return $success;
    }
}
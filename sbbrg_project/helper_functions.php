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
 * @param
 * @return boolean $success
 *
 * TODO: Troubleshoot function failure - possible refactor into multiple steps
 * TODO: Rethink Design - only listed as separate function to modularize code - not really necessary
 */
function generate_database_from_CSV(){
    #Initialize variable to keep track of status
    $success = false;

    # My SQL $user@'%' with $pass, granted all privileges on test.* used for this script.
    $user = 'auser';
    $pass = 'apass';

    # Links download separate csv files containing price data from Jan-1-2012 to Jan-31-2012 for GOOG and AAPL, respectively, from Google Finance Historical Prices.
    # Returns file contents as string - I believe/hope these are persistent links
    $googleData = file_get_contents("http://www.google.com/finance/historical?cid=304466804484872&startdate=Jan+1%2C+2012&enddate=Jan+31%2C+2012&num=30&ei=uG9xU7G_OdLeqAH3_IBY&output=csv");
    $yahooData = file_get_contents("http://www.google.com/finance/historical?cid=22144&startdate=Jan+1%2C+2012&enddate=Jan+31%2C+2012&num=30&ei=bnFxU5DIAo3cqQHk7oDoBA&output=csv");

    # Establish persistent connection to database, test on localhost. Prints error and returns false if fails.
    try{
        $dbh = new \PDO('mysql:host=localhost;dbname=test', $user, $pass, array(\PDO::ATTR_PERSISTENT => true));
        echo "Connected\n";
    } catch (\PDOException $e){
        echo "Unable to connect: " . $e->getMessage(). "\n";
        return $success;
    }

    # Add data to database tables - close connection at end and if errors
    try{
        # DB settings and connection settings and style adapted from PHP Data Objects Manual - https://php.net/manual/en/book.pdo.php
        $dbh->setAttribute(\PDO::ATTR_ERRMODE, \PDO::ERRMODE_EXCEPTION);

        # Set up parameters for data processing
        $fieldseparator = ",";
        $lineseparator = "\n";

        # Iterate through each data set
        foreach(array($googleData,$yahooData) as $csvdata){

            # Set table name based on which data set being used
            switch($csvdata){
                case($googleData):
                    $databasetable = 'google';
                    break;
                case($yahooData):
                    $databasetable = 'yahoo';
            }

            # Create My SQL table
            $dbh->beginTransaction();
            $dbh->exec("CREATE TABLE $databasetable (date DATE, open FLOAT, high FLOAT, low FLOAT, close FLOAT, volume INT)");
            $dbh->commit();

            # Prepare My SQL statement
            $stmt = $dbh->prepare("INSERT INTO $databasetable (date, open, high, low, close, volume) VALUES (:date, :open, :high, :low, :close, :volume)");

            # Bind parameters - there doesn't seem to be a binding for float or date - utilizing default of string.
            # TODO: Figure out proper way to bind float and date parameters...assuming this isn't it.
            $stmt->bindParam(':date', $date);
            $stmt->bindParam(':open', $open);
            $stmt->bindParam(':close', $close);
            $stmt->bindParam(':high', $high);
            $stmt->bindParam(':low', $low);
            $stmt->bindParam(':volume', $volume, \PDO::PARAM_INT);

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

                    # Break up line by delimiter and assign variable names to each element
                    list($raw_date, $open, $high, $low, $close, $volume) = explode($fieldseparator, $line);

                    # Fix date formatting to be accepted by My SQL
                    # TODO: This should be done properly such that csv could contain more months/years.
                    $parts = explode('-', $raw_date);
                    $date = "2012-01-$parts[0]";

                    # Execute My SQL statement
                    $stmt->execute();
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
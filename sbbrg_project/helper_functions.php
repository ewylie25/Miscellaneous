<?php
/**
 * Liz Wylie
 * May 2014
 *
 * Contains functions for sbbrg correlation project.
 *
 * Contents:
 * - generate_DB_table_from_URL($database, $data_set_info, $company_name)
 *
 *
 */

namespace sbbrg_project;

/***
 * Generate My SQL table from URL.
 *
 * Function takes a database instance, information about the data set from parameters.json and the company name(data set).
 *  Generates a new table in the given database from url. Assumes url contents are formatted with line and field delimiters
 *  like a CSV or TSV file. Broken up into four steps, variables documented in code.
 *
 * @param Database $database
 * @param associative array $data_set_info
 * @param string $company_name
 *
 * @return boolean $success, true if no Exception occurs during execution of function
 * */
function generate_DB_table_from_URL($database, $data_set_info, $company_name){
    // Return boolean denoting if code executed error free
    $success = true;
    try{

        /*
         * 1.GENERATING PARAMETERS
         * - $definition -> String, My SQL variable definitions(name and type) generated from field info specified in parameters.json
         * - $field_separator -> String, field delimiter as specified for data set in parameters.json
         * - $line_separator -> String, line delimiter as specified for data set in parameters.json
         * - $field_string -> String, My SQL variable names generated from field info specified in parameters.json
         * - $place_holders -> String, '?' placeholders for prepared My SQL statement, number of '?' = number of variables specified in parameters.json
         */

        // Generate parameters for creating table
        $definition = "";
        foreach($data_set_info["fields"] as $field){
            $definition = $definition." ".$field["name"]." ".$field["type"].",";
        }
        $definition = trim($definition, ',');

        // Generate parameters for writing to table
        $field_separator = $data_set_info["field delimiter"];
        $line_separator = $data_set_info["line delimiter"];
        $field_string = "";
        foreach($data_set_info["fields"] as $field){
            $field_string = $field_string." ".$field["name"].",";
        }
        $field_string = trim($field_string, ',');
        $place_holders = implode(',', array_fill(0, count($data_set_info["fields"]), '?'));

        /*
         * 2.RETRIEVING DATA:
         * - $csv_data -> String of file contents
         */

        // Returns contents of data files from url specified in parameters.json
        $csv_data = file_get_contents($data_set_info["url"]);

        /*
         * 3.CREATING TABLE:
         */
        $database->createTable($company_name, $definition);

        /*
         * 4. WRITING TO TABLE:
         * - $lines -> Array of Strings
         * - $line -> String
         * - $params -> Array of Strings
         */

        // Prepare My SQL statement
        $database->startAddToTable($company_name,$field_string,$place_holders);

        // Break up text into lines
        $lines = explode($line_separator, $csv_data);

        // Remove the header - beware indexing is off using this method - there is no longer an element in $lines[0]
        unset($lines[0]);

        // Iterate through lines of csv file contents
        foreach($lines as $line) {
            // Make sure the string isn't empty
            if (!empty($line)){
                // Remove any artifacts from file formatting
                $line = trim($line," \t");
                $line = str_replace("\r","",$line);

                // Break up line into values
                $params = explode($field_separator, $line);

                // Fix date formatting to be accepted by My SQL
                $params[0] = strftime("%Y-%m-%d", strtotime($params[0]));

                // Execute My SQL statement
                $database->addToTable($params);
            } else {
                // If empty sting - move on
                continue;
            }
        }
        // Close My SQL Transaction
        $database->endAddToTable();
        return $success;

    } catch(\Exception $e){
        // If an Exception occurs, print error message and return false
        $success=false;
        print "Unable to process data: " . $e->getMessage(). "\n";
        return $success;
    }
}
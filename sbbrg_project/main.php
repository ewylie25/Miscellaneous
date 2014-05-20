<?php
/**
 * Created by PhpStorm.
 * User: liz
 * Date: 5/20/14
 * Time: 12:32 PM
 */

include 'helper_functions.php';

$success = \sbbrg_project\generate_database_from_CSV();
if (!$success){
    die("Could not generate databases");
}

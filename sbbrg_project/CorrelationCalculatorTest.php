<?php
/**
 * User: liz
 * Date: 5/21/14
 * Created this test in order to ensure I did not break existing results in continued work
 * TODO: Better unit test
 * TODO: Some of the higher level unit test should either integrate into parameters.json or run off of a dummy data set.
 */

namespace sbbrg_project\Test;

require_once('CorrelationCalculator.php');
require_once('Database.php');

use sbbrg_project\CorrelationCalculator;
use sbbrg_project\Database;

class CorrelationCalculationTest extends \PHPUnit_Framework_TestCase {
    public $parameters;
    public $db;

    public function setUp(){
        $raw = file_get_contents('parameters.json');
        $this->parameters = json_decode($raw, true);
        $this->db = new Database($this->parameters["My SQL"]);
        $this->db->connect();
    }
    public function tearDown(){
        $this->db->disconnect();
    }

    public function testGetCorrelationAll(){
        $temp= new CorrelationCalculator($this->db);
        $temp->setCompanies(array('google', 'apple'));
        $temp->setDateRange("2012-01-01", "2012-01-31");

        $expectedResults = array('open' => -0.80021353044614,
                                'close' => -0.74911580407734,
                                'high' => -0.8003163877059,
                                'low' => -0.789410585227);
        foreach($expectedResults as $column =>$expectedResult){
            $temp->setCorrelation($column);
            $result = $temp->getCorrelation();
            $this->assertEquals($expectedResult, $result);
        }

    }
}
 
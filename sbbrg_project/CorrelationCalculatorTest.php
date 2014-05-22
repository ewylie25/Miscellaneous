<?php
/**
 * User: liz
 * Date: 5/21/14
 * Created this test in order to ensure I didn't break existing results in continued work
 * TODO: Better unit test
 */

namespace sbbrg_project\Test;

require_once('CorrelationCalculator.php');

use sbbrg_project\CorrelationCalculator;

class CorrelationCalculationTest extends \PHPUnit_Framework_TestCase {
    public $parameters;

    public function setUp(){
        $raw = file_get_contents('parameters.json');
        $this->parameters = json_decode($raw, true);
    }
    public function testGetCorrelationAll(){
        $expectedResults = array('open' => -0.80021353044614,
                                'close' => -0.74911580407734,
                                'high' => -0.8003163877059,
                                'low' => -0.789410585227);
        foreach($expectedResults as $column =>$expectedResult){
            $temp= new CorrelationCalculator($this->parameters);
            $temp->setCorrelation($column, "2012-01-01", "2012-01-31");
            $result = $temp->getCorrelation();
            $this->assertEquals($expectedResult, $result);
        }
    }
}
 
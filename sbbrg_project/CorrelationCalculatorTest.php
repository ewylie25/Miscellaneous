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

    public function testStats(){
        //
    }
    public function testGetCorrelationAll(){
        $expectedResults = array('open' => -0.80021353044614,
                                'close' => -0.74911580407734,
                                'high' => -0.8003163877059,
                                'low' => -0.789410585227);
        foreach($expectedResults as $column =>$expectedResult){
            $temp= new CorrelationCalculator($column);
            $result = $temp->getCorrelation();
            $this->assertEquals($expectedResult, $result);
        }
    }
}
 
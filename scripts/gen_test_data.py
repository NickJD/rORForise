import random
from rORForise_Testing.src.rORForise import check_pred as cp
import unittest

CDS_start = 1000
CDS_end = 2000


class TestPreds(unittest.TestCase):
    def test_correct_left_plus_plus_plus(self):
        CDS_direction = "+"
        read_direction = '+'

        # read overhangs left
        # correct start, correct frame, correct direction
        read_start = random.randint(CDS_start-100, CDS_start-1)
        read_end = random.randint(read_start+200, read_start+400)
        pred_direction = '+'
        pred_start = CDS_start - read_start + 1
        pred_end = read_end - read_start + 1
    
        answer = cp.check_pred(CDS_start, CDS_end, CDS_direction, read_start, read_end, read_direction, pred_start, pred_end, pred_direction)
        self.assertEqual(answer, set({0,8,10})) # correct

        
    def test_alt_start_left_plus_plus_plus(self):
        CDS_direction = "+"
        read_direction = '+'

        # read overhangs left
        # alternate start, correct frame, correct direction
        read_start = random.randint(CDS_start-100, CDS_start-1)
        read_end = random.randint(read_start+200, read_start+400)
        pred_direction = '+'
        pred_start = CDS_start - read_start + 1 + 3 # alt start
        pred_end = read_end - read_start + 1
    
        answer = cp.check_pred(CDS_start, CDS_end, CDS_direction, read_start, read_end, read_direction, pred_start, pred_end, pred_direction)
        self.assertEqual(answer, set({1,8,10})) # alt start


    def test_pred_before_CDS_plus_plus_plus(self):
        CDS_direction = "+"
        read_direction = '+'

        # read overhangs left
        # alternate start, correct frame, correct direction, pred before CDS
        read_start = CDS_start-100
        read_end = read_start+400
        pred_direction = '+'
        pred_start = 11
        pred_end = pred_start + 33

        answer = cp.check_pred(CDS_start, CDS_end, CDS_direction, read_start, read_end, read_direction, pred_start, pred_end, pred_direction)
        self.assertEqual(answer, set({1,8,10,12})) 
    

if __name__ == '__main__':
    unittest.main()

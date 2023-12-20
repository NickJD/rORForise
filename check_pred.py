"""
Example category 1

CDS: 686 1828 +
Read: 596 811 +
Pred (ORF2): 91 213 +

Read end = 811- 596 + 1 = 216, read start = 1  
Pred start = 91 so  this is 596+91-1 = 686 in genome/CDS numbers, and this is the correct start codon 
Pred end = 213
"""


def check_pred(cds_open, cds_close, cds_direction,
               read_open, read_close, read_direction,
               pred_start, pred_end, pred_direction):

    category = (cds_direction, read_direction, pred_direction)

    if category == ('+','+','+'): 
        read_start = 1
        read_end   = read_close - (read_open - 1)

        pred_cds_start = read_open + (pred_start - 1)
        pred_cds_end   = read_close + (pred_end - 1)

        start_diff = pred_cds_start - cds_open
        stop_diff  = pred_cds_end - cds_close

        read_captures_cds_start = read_open < cds_open
        read_captures_cds_end   = read_close > cds_close

        if read_captures_cds_start:
            if start_diff == 0:
                print("correct start")
            elif start_diff % 3 == 0:
                print("missed start")
            else:
                print("incorrect start")

        if read_captures_cds_end:
            if stop_diff == 0:
                print("correct stop")
            elif stop_diff % 3 == 0:
                print("missed stop")
            else:
                print("incorrect stop")

        elif not (read_captures_cds_end or read_captures_cds_start):
            # check the prediction fills the read? can we check frame?
            print("spanning")

    elif category == ('+','-','-'):
        print("case 3")

    elif category == ('-','+','-'):
        print("case 6")

    elif category == ('-','-','+'):
        print("case 8")



        
        
def main():
    check_pred(686,1828,'+', 596,811,'+', 91,213,'+') # cat 1, left


main()

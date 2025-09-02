from rORForise_Testing.src.rORForise.check_pred import check_pred
from rORForise_Testing.src.rORForise.check_pred import inverse_answers

##################################### CAT 1

#### LEFT
answer_vals = check_pred(686, 1828, '+', 580, 786, '+', 107, 205, '+')  # cat 1, left
print(answer_vals)  # This should be correct
for a in answer_vals:
    print(inverse_answers[a])
print()
#
answer_vals = check_pred(686, 1828, '+', 604, 811, '+', 2, 76, '+')  # cat 1, left - incorrect
print(answer_vals)  # This should be failing because the pred ends before the cds start
for a in answer_vals:
    print(inverse_answers[a])
print()
#
answer_vals = check_pred(686, 1828, '+', 604, 811, '+', 83, 205, '+')  # cat 1, left
print(answer_vals)  # same read as above but it starts within the CDS  - I have confirmed the start is correct and the end is the run off (-2 bases)
for a in answer_vals:
    print(inverse_answers[a])
print()


#### MID
answer_vals = check_pred(686, 1828, '+', 1281, 1473, '+', 3, 191, '+')  # cat 1, mid
print(answer_vals)  # Correct - Cuts off first 2 nts and last 1 nts. Makes sense as it is frame 3/6
for a in answer_vals:
    print(inverse_answers[a])
print()


answer_vals = check_pred(686, 1828, '+', 1341, 1559, '+', 3, 215, '+')  # cat 1, mid
print(answer_vals)  # Correct
for a in answer_vals:
    print(inverse_answers[a])
print()
# This is an interesting one -
# Its a frame 3 prediction which is correct in that it get the right frame to be a middle prediction but does not include its last 3 base codon onto the end.
#The AA from FGS ends with 'NELKDALQRIQTLAQNE' but it should be 'NELKDALQRIQTLAQNER'. The last codon 'AGA' is cut off by FGS.
#The first 2 bases are not included (makes sense as pred start is 3) but the pre end should continue until the end of the read...??


answer_vals = check_pred(686, 1828, '+', 1355, 1538, '+', 1, 180, '+')  # cat 1, mid
print(answer_vals)  # correct but again FGS did not include the last 3 bases - a whole codon?? In this case it was 'CAA'
for a in answer_vals:
    print(inverse_answers[a])
print()
### I couldn't find any case 1 examples of incorrect MIDs which I think kinda makes sense.


#### RIGHT
answer_vals = check_pred(686, 1828, '+', 1678, 1885, '+', 2, 151, '+')  # cat 1, right
print(answer_vals)  # Correct
for a in answer_vals:
    print(inverse_answers[a])
print()


answer_vals = check_pred(686, 1828, '+', 1783, 1997, '+', 46, 213, '+')  # cat 1, right
print(answer_vals)  # Incorrect
for a in answer_vals:
    print(inverse_answers[a])
print()
#Correct Frame is frame 2 at the very end of the read but a longer frame 1 was chosen instead.
#Correct prediction would have been from the second base of the read (Skip the first 'T' and start with 'AAA') \\
# then until 'TAA' at position 45, there I think the incorrect prediction starts at 46.



answer_vals = check_pred(686, 1828, '+', 1739, 1930, '+', 90, 188, '+')  # cat 1, right
print(answer_vals)  # Incorrect
for a in answer_vals:
    print(inverse_answers[a])
print()
# Correct Frame is 1. Prediction should have started at position 1, the first codon AAT and ended with 'TAA' at position 90.
# The prediction started at the 'A'TG where the correct prediction would have ended with the same position at TA'A'. (TA'A'TG).
# I checked my tool and my tool fails at the same read. Because 90-188 is longer than 1-90. Makes sense.



################################# CAT 2 TODO

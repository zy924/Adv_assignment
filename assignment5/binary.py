"""
Return a set of all binary strings of length n which contains k of zero
"""

from itertools import repeat
from itertools import chain
from itertools import permutations
def zbits(n,k):
	zeros = repeat('0',k) #create k of zero
	ones = repeat('1',n-k) #create n-k of one
	binary_list = list(chain(zeros,ones)) #create a list of length n binary string
	all_binary = set(["".join(item) for item in permutations(binary_list,n)])
	return all_binary

if __name__ == "__main__":
	assert zbits(4, 3) == {'0001', '0010', '0100', '1000'}
	assert zbits(4, 1) == {'0111', '1011', '1101', '1110'}
	assert zbits(5, 4) == {'00001', '00100', '01000', '10000', '00010'}
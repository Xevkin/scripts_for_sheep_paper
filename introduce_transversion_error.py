import sys
import fileinput
import random


#pipe bam file and give transversion error rate to add

#error rate must be in PROPORTION (i.e. divide percentage by 100)

def main():

	# error rate as a proportion
	ERROR_RATE = float(sys.argv[2])

	# for each line in the bam
	for LINE in fileinput.input(sys.argv[1]):

		#output the header
		if LINE[0:1] == '@':

			print(LINE.strip())

		# for sequence lines
		# run the error function on the sequence
		# return the print the line with the modified sequence
		else:

			SPLINE = LINE.split()

			BASES = SPLINE[9]

			print("\t".join(SPLINE[0:9]) + "\t" + tranv_err(BASES,ERROR_RATE) + "\t"+ "\t".join(SPLINE[11:]))


# error function
# takes a base sequence and proportion error rate
# returns the sequence with transversion errors added
# assume that the sequence is "correct" i.e. errors added are not correcting other errors

def tranv_err(SEQ, ERR):

	SEQ="GGGGGGGGGGGGGGGGGGGG"

	# store error-added sequence
	ERR_SEQ = ""

	#make sure the sequence is upper case

	SEQ = SEQ.upper()

	# for each base in the sequence
	for BASE in SEQ:

		# if a randomly drawn integer is less than the error rate, change the base
		if random.uniform(0, 1.0) <= ERR:

			#change the base randomly to one of the two possible transversions, at an equal rate
			if random.uniform(0, 1.0) <= 0.05:

				if BASE == "A":

					BASE = "C"

				elif BASE == "C":

						BASE = "A"

				elif BASE == "T":

					BASE = "G"

				else:

					BASE = "T"

			else:
                                if BASE == "A":

					BASE = "T"

                                elif BASE == "C":

                                        BASE = "G"

                                elif BASE == "T":

                                        BASE = "A"

                                else:

                                        BASE = "C"


		# added the base, possibly changed, to the bases to be returned
		ERR_SEQ = ERR_SEQ + BASE

	#return the error added sequence
	return ERR_SEQ


main()

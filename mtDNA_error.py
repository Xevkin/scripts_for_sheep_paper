#!/bin/python

import sys

import re

MATCH_COUNT = 0

MISMATCH_COUNT = 0

TRANS_MISMATCH_COUNT = 0

# input are a pileup file aligned to the mitochondrial fasta, and the length of the mitochondrial fasta
# the consensus base is inferred from the data, in order to determine the rate of minority bases

with open(sys.argv[1]) as FILE:

	for LINE in FILE:

		# dummy variable in case we don't call a base at a site

		CALL="_"

		SPLINE = LINE.rstrip("\t").split("\t")

		# check in case there is an issue with the positional information - go no further than length of mtDNA

		if int(SPLINE[1]) >= int(sys.argv[2]):

			continue

		BASES = SPLINE[4]

		#clean up the base string
		#we filter out insertions and deletions
		BASES = re.sub("\^[ABCDEFGHIJKLMNO]","", BASES)

		BASES = re.sub("[+-]1[ACGTcgat]", "", BASES)

		BASES = re.sub("[+-]2[ACGTcgat][ACGTcgat]", "", BASES)

		BASES = re.sub("[+-]3[ACGTcgat][ACGTcgat][ACGTcgat]", "", BASES)

		BASES =	BASES.replace("$","").replace('*',"").replace(",",".").replace(".",SPLINE[2]).upper()

		#filter to keep sites with at max 2 nucleotides

		if  len(BASES) < 3:

			continue

		#count the number of each nucleotide
		A_count = BASES.count("A")

		T_count = BASES.count("T")

		C_count = BASES.count("C")

		G_count = BASES.count("G")

		#for each nucleotide, if that nucleotide is the most commonly observed nucelotide at that site, get a count of the other observed bases

		if all(i < A_count for i in [T_count, C_count, G_count]):

			# in this instance we assume A is the true nucleotide which would be called in our mtDNA fasta using a consensus approach
			CALL = "A"

			# track the cumulative match count of bases
			MATCH_COUNT = MATCH_COUNT + A_count

			#tracking the cumulative mismatch count, and cumulative transversion count
			MISMATCH_COUNT = MISMATCH_COUNT + G_count + C_count + T_count

			TRANS_MISMATCH_COUNT = TRANS_MISMATCH_COUNT + G_count

			continue

		# repeat for T being the most common nucleotide, etc
                if all(i < T_count for i in [A_count, C_count, G_count]):

                        CALL = "T"

                        MATCH_COUNT = MATCH_COUNT + T_count

                        MISMATCH_COUNT = MISMATCH_COUNT +  G_count + C_count + A_count

			TRANS_MISMATCH_COUNT = TRANS_MISMATCH_COUNT + A_count

			continue

                if all(i < C_count for i in [T_count, A_count, G_count]):

                        CALl = "C"

                        MATCH_COUNT = MATCH_COUNT + C_count

                        MISMATCH_COUNT = MISMATCH_COUNT +  G_count + A_count + T_count

			TRANS_MISMATCH_COUNT = TRANS_MISMATCH_COUNT + T_count

			continue

                if all(i < G_count for i in [T_count, C_count, A_count]):

                        CALL = "G"

                        MATCH_COUNT = MATCH_COUNT + G_count

                        MISMATCH_COUNT = MISMATCH_COUNT +  A_count + C_count + T_count

			TRANS_MISMATCH_COUNT = TRANS_MISMATCH_COUNT + A_count

			continue

		# now for sites with more than one most common base e.g. 2 A, 2T
		# remove our "true" base from our observed bases
		BASES.replace(CALL, "")

		# get a list of the unique nucleotides remaining
		BASES_SET=set(list(BASES))

		# if we see more than one bases after removing our consensus call
		if len(BASES_SET) > 1:

			if CALL != "_":

				NUCS = ["A",  "G", "T", "C"].remove(CALL)

			if len(BASES_SET) == 2:

				#if we have two bases which are not the consensus base with the same count, choose one to add to the match count and use the same value in the mismatch count
				MATCH_COUNT = MATCH_COUNT + BASES.count(max(BASES))

				MISMATCH_COUNT = MISMATCH_COUNT + BASES.count(max(BASES))


				# if the two bases are transition, add to the transition couunt
				if BASES_SET in [set(['A', 'G']), set(['C', 'T'])]:

					TRANS_MISMATCH_COUNT = TRANS_MISMATCH_COUNT + BASES.count(max(BASES))

			#if we have 3 or 4 possible bases, do the same but add all other base counts to the mismatch count. these sites are not included in the transition count
			elif len(BASES_SET) >= 3:

				MATCH_COUNT = MATCH_COUNT + BASES.count(max(BASES))

				BASES_REMOVED = BASES.replace(max(BASES),"")

				MISMATCH_COUNT = MISMATCH_COUNT + sum([BASES_REMOVED.count(i) for i in ["A",  "G", "T", "C"]])


	#return the matchs, mismatches, mismatches minus transitions i.e. transversions, and the total number of sites
	print(sys.argv[1].split("/")[-1].split("_")[0] + "," + str(MATCH_COUNT) + "," + str(MISMATCH_COUNT) + "," +  str(MISMATCH_COUNT - TRANS_MISMATCH_COUNT)  + "," +  str(MATCH_COUNT + MISMATCH_COUNT))


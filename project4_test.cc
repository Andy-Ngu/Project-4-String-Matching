///////////////////////////////////////////////////////////////////////////////
// maxprotein_test.cc
//
// Unit tests for maxprotein.hh
//
///////////////////////////////////////////////////////////////////////////////


#include <cassert>
#include <sstream>

#include "project4.hh"
#include "rubrictest.hh"
#include "timer.hh"

int main() {
	Rubric rubric;

	ProteinVector proteins;
	auto load_successful = load_proteins(proteins, "proteins_large.txt");
	assert( load_successful );

    rubric.criterion("load_protein still works", 1, [&]() {
		     TEST_EQUAL("size", 6639, proteins.size());
		   });

    BlosumPenaltyArray bpa;
	load_successful = load_blosum_file(bpa, "blosum62.txt");
	assert( load_successful );

	std::string align_string1;
	std::string align_string2;

  rubric.criterion("local_alignment trivial equal cases", 2,
		   [&]() {
		     auto soln = local_alignment("", "", bpa, align_string1, align_string2);
		     TEST_EQUAL("same", 0, soln);
		     soln = local_alignment("A", "A", bpa, align_string1, align_string2);
		     TEST_EQUAL("same 1", 4, soln);
		     soln = local_alignment("AB", "AB", bpa, align_string1, align_string2);
		     TEST_EQUAL("same 2", 8, soln);
		     soln = local_alignment("ABC", "ABC", bpa, align_string1, align_string2);
		     TEST_EQUAL("same 3", 17, soln);
		   });
  
  rubric.criterion("local_alignment trivial deletion/insertion cases", 2,
		   [&]() {
		     auto soln = local_alignment("A", "", bpa, align_string1, align_string2);
		     TEST_EQUAL("deletion A", 0, soln);
		     soln = local_alignment("", "A", bpa, align_string1, align_string2);
		     TEST_EQUAL("insertion A", 0, soln);
		     soln = local_alignment("ABC", "AC", bpa, align_string1, align_string2);
		     TEST_EQUAL("deletion B", 9, soln);
		     soln = local_alignment("AC", "ABC", bpa, align_string1, align_string2);
		     TEST_EQUAL("insertion B", 9, soln);
		   });
  
  rubric.criterion("local_alignment trivial substitution cases", 2,
		   [&]() {
		     auto soln = local_alignment("ABC", "XBC", bpa, align_string1, align_string2);
		     TEST_EQUAL("substitution 1", 13, soln);
		     soln = local_alignment("XBC", "ABC", bpa, align_string1, align_string2);
		     TEST_EQUAL("substitution 1", 13, soln);
		     soln = local_alignment("AXC", "ABC", bpa, align_string1, align_string2);
		     TEST_EQUAL("substitution 2", 12, soln);
		     soln = local_alignment("ABC", "AXC", bpa, align_string1, align_string2);
		     TEST_EQUAL("substitution 2", 12, soln);
		   });
  
  rubric.criterion("local_alignment simple local cases", 2,
		   [&]() {
		     auto soln = local_alignment("DEFG", "ABCDEFGHIJ", bpa, align_string1, align_string2);
		     TEST_EQUAL("alignment string 1a", "DEFG", align_string1);
		     TEST_EQUAL("alignment string 1b", "DEFG", align_string2);
		     soln = local_alignment("DEGH", "ABCDEFGHIJ", bpa, align_string1, align_string2);
		     TEST_EQUAL("alignment string 1a", "DE*GH", align_string1);
		     TEST_EQUAL("alignment string 1b", "DEFGH", align_string2);
		     soln = local_alignment("DEXFG", "ABCDEFGHIJ", bpa, align_string1, align_string2);
		     TEST_EQUAL("alignment string 1a", "DEXFG", align_string1);
		     TEST_EQUAL("alignment string 1b", "DE*FG", align_string2);
		   });
  
  return rubric.run();
}


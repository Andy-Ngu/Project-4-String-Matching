
#include "project4.hh"
#include "timer.hh"

int main() {
	ProteinVector proteins;
	auto load_successful = load_proteins(proteins, "proteins_large.txt");
	assert( load_successful );

    	BlosumPenaltyArray bpa;
	load_successful = load_blosum_file(bpa, "blosum62.txt");
	assert( load_successful );

	std::vector<std::string> testProteins;
	testProteins.push_back("PIEPCMGA");
	testProteins.push_back("TQGASNIGE");
	testProteins.push_back("ALAKLIRYGG");
	testProteins.push_back("CSNPNLSDFGR");
	testProteins.push_back("MYPEPTIDE");

	std::cout << "------------------- Dynamic Programming ------------------" << std::endl;
	for (int i = 0; i < testProteins.size(); i++) {
		std::string align_string1;
		std::string align_string2;
		std::string searchString = 	testProteins[i];
		Timer timer;
		std::cout << "String to Match = " << testProteins[i] << std::endl;
		std::shared_ptr<Protein> best_protein = local_alignment_best_match(proteins, searchString, 
		bpa, align_string1, align_string2);
		std::cout << best_protein->description << std::endl;
		std::cout << align_string1 << std::endl;
		std::cout << align_string2 << std::endl;
		std::cout << timer.elapsed() << std::endl;
	}

  return 0;
}

